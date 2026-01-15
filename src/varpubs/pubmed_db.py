from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Iterable
import logging
import sqlalchemy
from sqlmodel import Field, SQLModel, Session, Integer
from varpubs.hgvs_extractor import extract_hgvsp_from_vcf
from sqlmodel import select
from pubgator import PubGator

logger = logging.getLogger(__name__)


class PubmedArticle(SQLModel, table=True):
    """Normalized PubMed article metadata stored locally."""

    pmid: int = Field(Integer, primary_key=True, nullable=False)
    title: str
    abstract: str
    authors: str
    journal: str
    doi: str


class BioconceptToPMID(SQLModel, table=True):
    """Mapping between extracted bioconcept (e.g. @VARIANT_p.V560del_KIT_human) and PMID."""

    bioconcept: str = Field(nullable=False, primary_key=True)
    pmid: int = Field(Integer, nullable=False, primary_key=True)


@dataclass
class PubmedDB:
    """End-to-end pipeline: extract terms → search PubMed → store articles → map terms."""

    path: Path
    vcf_paths: Iterable[Path]
    species: str
    _engine: Optional[sqlalchemy.engine.base.Engine] = field(init=False, default=None)

    def deploy(self) -> None:
        """Create/update DB, fetch articles, and populate term-to-PMID links."""
        logger.info(f"Deploying database at: {self.path}")
        logger.info(f"Input VCF paths: {[str(p) for p in self.vcf_paths]}")

        pg = PubGator(max_requests_per_second=1)
        self.create_tables()

        relations = dict()

        bioconcepts = set()
        for vcf_path in self.vcf_paths:
            bioconcepts.update(extract_hgvsp_from_vcf(str(vcf_path), self.species))

        for bioconcept in bioconcepts:
            # Consider adding parameter to manually set max_ret
            retries = 5
            publications = list()
            for i in range(1, retries + 1):
                try:
                    publications = pg.search(bioconcept, max_ret=25)
                    break
                except ValueError:
                    logger.warn(
                        f"Failed search for: {bioconcept}. Try {i} of {retries}"
                    )
                    if i == retries:
                        raise ValueError("Max retries reached.")
            pmids = [publication.pmid for publication in publications]
            relations[bioconcept] = pmids

            logger.debug(f"Found {len(pmids)} PMIDs for {bioconcept}")

        with Session(self.engine) as session:
            pmids = [pmid for _, pmids in relations.items() for pmid in pmids]
            to_be_fetched = [
                pmid
                for pmid in pmids
                if not session.exec(
                    select(PubmedArticle).where(PubmedArticle.pmid == pmid)
                ).first()
            ]

            BATCH_SIZE = 100
            for i in range(0, len(to_be_fetched), BATCH_SIZE):
                batch = to_be_fetched[i : i + BATCH_SIZE]
                annotations = pg.export_publications(pmids=batch)
                for annotation in annotations.documents:
                    article = self.parse_article(annotation)
                    session.add(article)

            for bioconcept, pmids in relations.items():
                for pmid in pmids:
                    exists = session.exec(
                        select(BioconceptToPMID).where(
                            BioconceptToPMID.bioconcept == bioconcept,
                            BioconceptToPMID.pmid == pmid,
                        )
                    ).first()
                    if not exists:
                        session.add(BioconceptToPMID(bioconcept=bioconcept, pmid=pmid))

            session.commit()
            logger.info("Data committed to the database.")

    def create_tables(self) -> None:
        """Create tables if missing; safe to call multiple times."""
        self.path.parent.mkdir(parents=True, exist_ok=True)

        logger.info("Creating or updating tables...")
        SQLModel.metadata.create_all(self.engine)

        if not self.path.exists():
            logger.info("New database created.")
        else:
            logger.info("Existing database found. Tables checked or updated.")

    def parse_article(self, article) -> Optional[PubmedArticle]:
        title, authors, journal, doi, abstract = "", "", "", "", ""

        for passage in article.passages:
            if passage.infons.get("type") == "title":
                title = passage.text
                authors = passage.infons.get("authors", "")
                meta = passage.infons.get("journal")
                if meta:
                    journal, _, rest = meta.partition(";")
                    doi = rest.split("doi:", 1)[1] if rest and "doi:" in rest else ""
            if passage.infons.get("type") == "abstract":
                abstract = passage.text

        return PubmedArticle(
            pmid=article.id,
            title=title,
            authors=authors,
            journal=journal,
            doi=doi,
            abstract=abstract,
        )

    @property
    def engine(self) -> sqlalchemy.engine.base.Engine:
        """Lazy-create SQLAlchemy engine for DuckDB."""
        if self._engine is None:
            self._engine = sqlalchemy.create_engine(f"duckdb:///{self.path}")
        return self._engine
