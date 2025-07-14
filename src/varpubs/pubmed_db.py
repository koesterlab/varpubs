from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Iterable
import logging
import sqlalchemy
from sqlmodel import Field, SQLModel, Session, Integer
from varpubs.hgvs_extractor import extract_hgvsp_from_vcf
from Bio import Entrez
import re

# Set up logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PubmedArticle(SQLModel, table=True):
    """
    Represents a PubMed article record in the database.
    Stores metadata such as title, abstract, authors, journal, publication date, and DOI.
    """

    pmid: int = Field(Integer, primary_key=True, nullable=False)
    title: str
    abstract: str
    authors: str
    journal: str
    pub_date: str
    doi: str


class TermToPMID(SQLModel, table=True):
    """
    Maps a searched variant term (HGVS.p) to PubMed article PMIDs.

    This table captures which search term (e.g., "KRAS p.G12D") led to which PubMed article.
    It is designed with a composite primary key (term, pmid, gene) to prevent duplicate entries
    and to reflect the fact that a single term may map to multiple articles and vice versa.

    Attributes:
        term (str): The HGVS.p variant term used in the PubMed query.
        pmid (int): PubMed ID of the retrieved article.
        gene (str): Gene name parsed from the variant (e.g., "KRAS").
    """

    term: str = Field(nullable=False)
    pmid: int = Field(Integer, nullable=False)
    gene: str = Field(nullable=False)

    # Composite primary key to prevent duplicates
    __table_args__ = (sqlalchemy.PrimaryKeyConstraint("term", "pmid", "gene"),)


@dataclass
class PubmedDB:
    """
    Class responsible for building and managing a local DuckDB database
    that maps HGVS.p variant annotations to PubMed articles.
    """

    path: Path
    vcf_paths: Iterable[Path]
    email: str
    _engine: Optional[sqlalchemy.engine.base.Engine] = field(init=False, default=None)

    def deploy(self) -> None:
        """
        Main pipeline entry point.
        - Extracts terms from VCFs
        - Queries PubMed for those terms
        - Saves article metadata and mapping info to the database
        """
        logger.info(f"Deploying database at: {self.path}")
        logger.info(f"Input VCF paths: {[str(p) for p in self.vcf_paths]}")
        logger.info(f"Entrez email: {self.email}")

        self.create_tables()
        terms = self.extract_terms()
        Entrez.email = self.email

        with Session(self.engine) as session:
            for term in terms:
                self.process_term(term, session)
            session.commit()

        logger.info("Data committed to the database.")

    def create_tables(self) -> None:
        """
        Creates required tables in the DuckDB database if they don't already exist.
        """
        if self.path.parent != Path():
            self.path.parent.mkdir(parents=True, exist_ok=True)

        logger.info("Creating or updating tables...")
        SQLModel.metadata.create_all(self.engine)

        if not self.path.exists():
            logger.info("New database created.")
        else:
            logger.info("Existing database found. Tables checked or updated.")

    def extract_terms(self) -> set[str]:
        """
        Extracts all unique variant terms (HGVS.p) from provided VCF files.

        Returns:
            set[str]: Unique set of variant terms for querying PubMed.
        """
        terms = set()
        for vcf_path in self.vcf_paths:
            terms.update(extract_hgvsp_from_vcf(str(vcf_path)))
        return terms

    def process_term(self, term: str, session: Session) -> None:
        """
        For a given variant term, searches PubMed and stores the article metadata
        and the mapping (term → pmid) into the database.

        This method performs the following:
        - Queries PubMed using the variant term.
        - Parses the gene name from the term (e.g., "KRAS" from "KRAS p.G12D").
        - Checks if the mapping (term, pmid, gene) already exists in the TermToPMID table.
        - If the article corresponding to a PMID is not already stored in PubmedArticle,
        fetches metadata and stores it.
        - Saves the term-to-article mapping in TermToPMID for traceability.
        """
        logger.info(f"Searching PubMed for: {term}")
        try:
            handle = Entrez.esearch(db="pubmed", term=term, retmax=3)
            record = Entrez.read(handle)
            ids = record.get("IdList", [])
            if not ids:
                logger.info(f"No PubMed results for: {term}")
                return

            match = re.match(r"(\w+)\sp\.", term)
            gene = match.group(1) if match else "UNKNOWN"

            for pmid in map(int, ids):
                # Skip if this mapping already exists
                mapping_exists = session.exec(
                    sqlalchemy.select(TermToPMID).where(
                        TermToPMID.term == term,
                        TermToPMID.pmid == pmid,
                        TermToPMID.gene == gene,
                    )
                ).first()
                if mapping_exists:
                    logger.info(
                        f"Already mapped: (term={term}, PMID={pmid}, gene={gene})"
                    )
                    continue
                # Add article metadata if it's not already stored
                article_exists = session.exec(
                    sqlalchemy.select(PubmedArticle).where(PubmedArticle.pmid == pmid)
                ).first()
                if not article_exists:
                    try:
                        fetch = Entrez.efetch(
                            db="pubmed",
                            id=str(pmid),
                            rettype="medline",
                            retmode="xml",
                        )
                        article_data = Entrez.read(fetch)
                        pubmed_articles = article_data.get("PubmedArticle", [])
                        if not pubmed_articles:
                            logger.warning(f"No PubmedArticle found for PMID {pmid}")
                            continue
                        article = (
                            pubmed_articles[0]
                            .get("MedlineCitation", {})
                            .get("Article", {})
                        )
                        parsed = self.parse_article(article, pmid)
                        if parsed:
                            session.add(parsed)
                    except Exception as e:
                        logger.error(f"Failed to fetch article for PMID {pmid}: {e}")
                        continue
                # Store the mapping term → pmid with the exact searched term
                session.add(TermToPMID(term=term, pmid=pmid, gene=gene))

        except Exception as e:
            logger.error(f"Skipping term {term} due to error: {e}")

    def parse_article(self, article: dict, pmid: int) -> Optional[PubmedArticle]:
        """
        Converts raw article data from Entrez into a PubmedArticle DB record.

        Args:
            article (dict): Article data parsed from PubMed.
            pmid (int): PubMed ID for the article.

        Returns:
            Optional[PubmedArticle]: A PubmedArticle object or None if parsing failed.
        """
        try:
            title = article.get("ArticleTitle", "")

            # Extract abstract text
            abstract = ""
            try:
                abstract_obj = article.get("Abstract", {})
                abstract_text = abstract_obj.get("AbstractText", "")
                if isinstance(abstract_text, list):
                    abstract = " ".join(str(p) for p in abstract_text)
                elif isinstance(abstract_text, str):
                    abstract = abstract_text
                else:
                    abstract = str(abstract_text)
            except Exception as e:
                logger.warning(f"Abstract parse failed for PMID {pmid}: {e}")
            # Extract author names
            authors = ""
            try:
                author_list = article.get("AuthorList", [])
                if isinstance(author_list, list):
                    authors = ", ".join(
                        f"{a.get('LastName', '')} {a.get('ForeName', '')}".strip()
                        for a in author_list
                        if isinstance(a, dict) and "LastName" in a
                    )
            except Exception as e:
                logger.warning(f"Author parse failed for PMID {pmid}: {e}")

            journal = article.get("Journal", {}).get("Title", "")
            pub_date = (
                article.get("Journal", {})
                .get("JournalIssue", {})
                .get("PubDate", {})
                .get("Year", "")
            )
            # Extract DOI if available
            doi = ""
            try:
                eloc = article.get("ELocationID", [{}])
                if isinstance(eloc, list) and eloc and isinstance(eloc[0], dict):
                    doi = eloc[0].get("#text", "")
            except Exception as e:
                logger.warning(f"DOI parse failed for PMID {pmid}: {e}")

            return PubmedArticle(
                pmid=pmid,
                title=title,
                abstract=abstract,
                authors=authors,
                journal=journal,
                pub_date=pub_date,
                doi=doi,
            )
        except Exception as e:
            logger.error(f"Failed to parse article for PMID {pmid}: {e}")
            return None

    @property
    def engine(self) -> sqlalchemy.engine.base.Engine:
        """
        Lazily initializes and returns the SQLAlchemy engine for DuckDB.

        Returns:
            sqlalchemy.engine.base.Engine: Connection engine to the DuckDB database.
        """
        if self._engine is None:
            self._engine = sqlalchemy.create_engine(f"duckdb:///{self.path}")
        return self._engine
