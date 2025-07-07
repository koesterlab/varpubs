from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Iterable
import logging
import sqlalchemy
from sqlmodel import Field, SQLModel, Session, Integer
from sqlmodel import Field, SQLModel, Session, Integer
from varpubs.hgvs_extractor import extract_hgvsp_from_vcf
from Bio import Entrez

# Set up logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PubmedArticle(SQLModel, table=True):
    """Table to store metadata of PubMed articles."""
    pmid: int = Field(Integer, primary_key=True, nullable=False)
    title: str
    abstract: str
    authors: str
    journal: str
    pub_date: str
    doi: str


class TermToPMID(SQLModel, table=True):
    """Mapping table between HGVS.p terms and PubMed PMIDs."""
    term: str = Field(nullable=False)
    pmid: int = Field(Integer, primary_key=True, nullable=False)
    __table_args__ = (sqlalchemy.PrimaryKeyConstraint("term", "pmid"),)


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
        Main pipeline entry point. Extracts HGVS.p terms from VCFs,
        queries PubMed for article metadata, and populates the database.
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
        Creates required tables if they do not already exist.
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
        Extracts unique HGVS.p terms from the provided VCF files.

        Returns:
            set[str]: A set of extracted HGVS.p terms.
        """
        terms = set()
        for vcf_path in self.vcf_paths:
            terms.update(extract_hgvsp_from_vcf(str(vcf_path)))
        return terms

    def process_term(self, term: str, session: Session) -> None:
        """
        Processes a single HGVS.p term:
        - Queries PubMed for article PMIDs.
        - Parses and stores article metadata if not already stored.
        - Maps term to corresponding PMIDs.

        Args:
            term (str): The HGVS.p term to process.
            session (Session): The current DB session.
        """
        logger.info(f"Searching PubMed for: {term}")
        try:
            handle = Entrez.esearch(db="pubmed", term=term, retmax=3)
            record = Entrez.read(handle)
            ids = record.get("IdList", [])
            if not ids:
                logger.info(f"No PubMed results found for HGVS.p term: {term}")
                return

            for pmid in map(int, ids):
                mapping_exists = session.exec(
                    sqlalchemy.select(TermToPMID).where(
                        TermToPMID.term == term,
                        TermToPMID.pmid == pmid
                    )
                ).first()
                if mapping_exists:
                    logger.info(f"Skipping (term={term}, PMID={pmid}): already mapped.")
                    continue

                article_exists = session.exec(
                    sqlalchemy.select(PubmedArticle).where(
                        PubmedArticle.pmid == pmid
                    )
                ).first()

                if not article_exists:
                    fetch = Entrez.efetch(db="pubmed", id=str(pmid), rettype="medline", retmode="xml")
                    article_data = Entrez.read(fetch)
                    pubmed_articles = article_data.get("PubmedArticle", [])
                    if not pubmed_articles:
                        logger.warning(f"No PubmedArticle found for PMID {pmid}")
                        continue
                    article = pubmed_articles[0].get("MedlineCitation", {}).get("Article", {})
                    parsed = self.parse_article(article, pmid)
                    if parsed:
                        session.add(parsed)

                session.add(TermToPMID(term=term, pmid=pmid))
        except Exception as e:
            logger.error(f"Skipping term {term} due to error: {e}")

    def parse_article(self, article: dict, pmid: int) -> Optional[PubmedArticle]:
        """
        Parses article metadata into a PubmedArticle object.

        Args:
            article (dict): The raw article dictionary from Entrez.
            pmid (int): The PubMed ID of the article.

        Returns:
            Optional[PubmedArticle]: The parsed article, or None if parsing failed.
        """
        try:
            title = article.get("ArticleTitle", "")

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

            authors = ""
            try:
                author_list = article.get("AuthorList", [])
                if isinstance(author_list, list):
                    authors = ", ".join(
                        f"{a.get('LastName', '')} {a.get('ForeName', '')}".strip()
                        for a in author_list if isinstance(a, dict) and "LastName" in a
                    )
            except Exception as e:
                logger.warning(f"Author parse failed for PMID {pmid}: {e}")

            journal = article.get("Journal", {}).get("Title", "")
            pub_date = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {}).get("Year", "")
            doi = ""
            try:
                eloc = article.get("ELocationID", [{}])
                if isinstance(eloc, list) and eloc and isinstance(eloc[0], dict):
                    doi = eloc[0].get("#text", "")
            except:
                pass

            return PubmedArticle(
                pmid=pmid,
                title=title,
                abstract=abstract,
                authors=authors,
                journal=journal,
                pub_date=pub_date,
                doi=doi
            )
        except Exception as e:
            logger.error(f"Failed to parse article for PMID {pmid}: {e}")
            return None

    @property
    def engine(self) -> sqlalchemy.engine.base.Engine:
        """
        Returns the SQLAlchemy engine for the DuckDB database.

        Returns:
            sqlalchemy.engine.base.Engine: The database engine.
        """
        if self._engine is None:
            self._engine = sqlalchemy.create_engine(f"duckdb:///{self.path}")
        return self._engine
