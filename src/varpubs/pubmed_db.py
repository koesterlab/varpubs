from dataclasses import dataclass, field
from pathlib import Path
from typing import cast, Optional, Iterable
import logging
import sqlalchemy
from sqlmodel import Field, SQLModel, Session, Integer
from varpubs.hgvs_extractor import extract_hgvsp_from_vcf
from Bio import Entrez
import re
from sqlmodel import select

# Set up logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PubmedArticle(SQLModel, table=True):
    pmid: int = Field(Integer, primary_key=True, nullable=False)
    title: str
    abstract: str
    authors: str
    journal: str
    pub_date: str
    doi: str


class TermToPMID(SQLModel, table=True):
    term: str = Field(nullable=False)
    pmid: int = Field(Integer, nullable=False)
    gene: str = Field(nullable=False)

    __table_args__ = (sqlalchemy.PrimaryKeyConstraint("term", "pmid", "gene"),)


@dataclass
class PubmedDB:
    path: Path
    vcf_paths: Iterable[Path]
    email: str
    _engine: Optional[sqlalchemy.engine.base.Engine] = field(init=False, default=None)

    def deploy(self) -> None:
        logger.info(f"Deploying database at: {self.path}")
        logger.info(f"Input VCF paths: {[str(p) for p in self.vcf_paths]}")
        logger.info(f"Entrez email: {self.email}")

        self.create_tables()
        terms = self.extract_terms()
        Entrez.email = self.email

        with Session(self.engine) as session:
            all_pmids = set()
            term_to_pmids = {}

            for term in terms:
                try:
                    handle = Entrez.esearch(db="pubmed", term=term, retmax=3)
                    record = cast(dict, Entrez.read(handle))
                    pmids = list(map(int, record.get("IdList", [])))
                    if pmids:
                        match = re.match(r"(\w+)\sp\.", term)
                        gene = match.group(1) if match else "UNKNOWN"
                        term_to_pmids[term] = [(pmid, gene) for pmid in pmids]
                        all_pmids.update(pmids)
                except Exception as e:
                    logger.error(f"Error fetching PMIDs for term '{term}': {e}")

            # Fetch metadata in batch
            self.fetch_articles_metadata(session, list(all_pmids))

            # Store term→pmid mappings
            for term, pmid_gene_list in term_to_pmids.items():
                for pmid, gene in pmid_gene_list:
                    exists = session.exec(
                        select(TermToPMID).where(
                            TermToPMID.term == term,
                            TermToPMID.pmid == pmid,
                            TermToPMID.gene == gene,
                        )
                    ).first()
                    if not exists:
                        session.add(TermToPMID(term=term, pmid=pmid, gene=gene))

            session.commit()
            logger.info("Data committed to the database.")

    def create_tables(self) -> None:
        if self.path.parent != Path():
            self.path.parent.mkdir(parents=True, exist_ok=True)

        logger.info("Creating or updating tables...")
        SQLModel.metadata.create_all(self.engine)

        if not self.path.exists():
            logger.info("New database created.")
        else:
            logger.info("Existing database found. Tables checked or updated.")

    def extract_terms(self) -> set[str]:
        terms = set()
        for vcf_path in self.vcf_paths:
            terms.update(extract_hgvsp_from_vcf(str(vcf_path)))
        return terms

    def fetch_articles_metadata(self, session: Session, pmids: list[int]) -> None:
        BATCH_SIZE = 50
        for i in range(0, len(pmids), BATCH_SIZE):
            batch = pmids[i : i + BATCH_SIZE]
            try:
                fetch = Entrez.efetch(
                    db="pubmed",
                    id=",".join(map(str, batch)),
                    rettype="medline",
                    retmode="xml",
                )
                data = cast(dict, Entrez.read(fetch))
                for item in data.get("PubmedArticle", []):
                    medline = cast(dict, item.get("MedlineCitation", {}))
                    article = cast(dict, medline.get("Article", {}))
                    pmid = int(cast(str, medline.get("PMID", "0")))
                    exists = session.exec(
                        select(PubmedArticle).where(PubmedArticle.pmid == pmid)
                    ).first()
                    if not exists:
                        parsed = self.parse_article(article, pmid)
                        if parsed:
                            session.add(parsed)
            except Exception as e:
                logger.error(f"Batch fetch failed for PMIDs {batch}: {e}")

    def parse_article(self, article: dict, pmid: int) -> Optional[PubmedArticle]:
        try:
            title = article.get("ArticleTitle", "")
            abstract = ""
            try:
                abstract_obj = (
                    article.get("Abstract", {}) if isinstance(article, dict) else {}
                )
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
                        for a in author_list
                        if isinstance(a, dict) and "LastName" in a
                    )
            except Exception as e:
                logger.warning(f"Author parse failed for PMID {pmid}: {e}")

            journal = (
                article.get("Journal", {}).get("Title", "")
                if isinstance(article, dict)
                else ""
            )
            pub_date = (
                article.get("Journal", {})
                .get("JournalIssue", {})
                .get("PubDate", {})
                .get("Year", "")
            )
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
        if self._engine is None:
            self._engine = sqlalchemy.create_engine(f"duckdb:///{self.path}")
        return self._engine
