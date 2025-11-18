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
from varpubs.hgvs_extractor import AA3_TO_1

logger = logging.getLogger(__name__)


def _variant_synonyms(term: str):
    """Generate equivalent variant spellings (e.g., p.Ser339Leu → Ser339Leu, S339L)."""
    m = re.match(r"^([A-Za-z0-9_-]+)\s+(.*)$", term.strip())
    if not m:
        return None, []
    gene, rest = m.group(1), m.group(2)
    syns = set()

    m3 = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", rest)
    if m3:
        ref3, pos, alt3 = m3.groups()
        syns.add(f"p.{ref3}{pos}{alt3}")
        syns.add(f"{ref3}{pos}{alt3}")
        ref1, alt1 = AA3_TO_1.get(ref3), AA3_TO_1.get(alt3)
        if ref1 and alt1:
            syns.add(f"{ref1}{pos}{alt1}")  # S339L
    else:
        m1 = re.match(r"([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])", rest)
        if m1:
            syns.add("".join(m1.groups()))
        else:
            syns.add(rest.strip())

    return gene, list(syns)


def _normalize_tokens(text: str):
    """Lowercase and split text into alphanumeric tokens (punctuation-insensitive)."""
    return re.sub(r"[^a-z0-9]+", " ", text.lower()).split()


def matches_variant_loose(term: str, text: str, token_window: int = 40) -> bool:
    """Loose match: gene and any variant synonym must appear within a token window."""
    gene, syns = _variant_synonyms(term)
    if not gene:
        return False
    tokens = _normalize_tokens(text)
    if not tokens:
        return False

    gene_tok = gene.lower()
    pos_gene = [i for i, t in enumerate(tokens) if t == gene_tok]
    if not pos_gene:
        return False

    def seq_positions(seq_tokens):
        """Find start positions where seq_tokens appear consecutively in tokens."""
        L = len(seq_tokens)
        if L == 0:
            return []
        return [
            i for i in range(len(tokens) - L + 1) if tokens[i : i + L] == seq_tokens
        ]

    for s in syns:
        seq = _normalize_tokens(s)
        if not seq:
            continue
        starts = seq_positions(seq)
        if not starts:
            continue
        for g in pos_gene:
            for sidx in starts:
                if abs(g - sidx) <= token_window:
                    return True
    return False


class PubmedArticle(SQLModel, table=True):
    """Normalized PubMed article metadata stored locally."""

    pmid: int = Field(Integer, primary_key=True, nullable=False)
    title: str
    abstract: str
    authors: str
    journal: str
    pub_date: str
    doi: str


class TermToPMID(SQLModel, table=True):
    """Mapping between extracted VCF term and matched PMID (composite PK)."""

    term: str = Field(nullable=False, primary_key=True)
    pmid: int = Field(Integer, nullable=False, primary_key=True)
    gene: str = Field(nullable=False, primary_key=True)


@dataclass
class PubmedDB:
    """End-to-end pipeline: extract terms → search PubMed → store articles → map terms."""

    path: Path
    vcf_paths: Iterable[Path]
    email: str
    batch_size: int = 5
    _engine: Optional[sqlalchemy.engine.base.Engine] = field(init=False, default=None)

    def deploy(self) -> None:
        """Create/update DB, fetch articles, and populate term-to-PMID links."""
        logger.info(f"Deploying database at: {self.path}")
        logger.info(f"Input VCF paths: {[str(p) for p in self.vcf_paths]}")
        logger.info(f"Entrez email: {self.email}")

        self.create_tables()
        terms = self.extract_terms()
        Entrez.email = self.email

        pmids = []
        batches = [
            list(terms)[i : i + self.batch_size]
            for i in range(0, len(terms), self.batch_size)
        ]
        for batch in batches:
            # Concatenetante all variants in one batch via OR query to reduce API calls.
            logger.info(f"Querying PubMed for batch {batch}")
            query_terms = set()
            for term in batch:
                query_terms.add(term)
                gene, syns = _variant_synonyms(term)
                if gene:
                    for s in syns:
                        query_terms.add(f"{gene} {s}")
            query = " OR ".join(f'"{qterm}"' for qterm in query_terms)
            try:
                retmax = 10000
                handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax)
                record = cast(dict, Entrez.read(handle))
                fetched_pmids = list(map(int, record.get("IdList", [])))
                if len(fetched_pmids) == retmax:
                    logger.warning(
                        f"PubMed API limit reached for batch {batch}. Consider decreasing batch size parameter."
                    )
                pmids.extend(fetched_pmids)
            except Exception as e:
                logger.error(f"Error querying PubMed: {e}")
                return

        logger.info(f"Found {len(pmids)} PMIDs")

        with Session(self.engine) as session:
            """B) Batch-fetch article metadata for all PMIDs."""
            self.fetch_articles_metadata(session, pmids)

            """C) Build Term→PMID using loose matching on title+abstract."""
            term_to_pmids = {}
            for pmid in pmids:
                article = session.exec(
                    select(PubmedArticle).where(PubmedArticle.pmid == pmid)
                ).first()
                if not article:
                    continue
                title_abstract = article.title + " " + article.abstract

                for term in terms:
                    if matches_variant_loose(term, title_abstract):
                        gene = term.split()[0]
                        term_to_pmids.setdefault(term, []).append((pmid, gene))

            """D) Upsert mappings (composite PK prevents duplicates)."""
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
        """Create tables if missing; safe to call multiple times."""
        self.path.parent.mkdir(parents=True, exist_ok=True)

        logger.info("Creating or updating tables...")
        SQLModel.metadata.create_all(self.engine)

        if not self.path.exists():
            logger.info("New database created.")
        else:
            logger.info("Existing database found. Tables checked or updated.")

    def extract_terms(self) -> set[str]:
        """Extract HGVS protein terms and short forms from annotated VCFs."""
        terms = set()
        for vcf_path in self.vcf_paths:
            terms.update(extract_hgvsp_from_vcf(str(vcf_path)))
        return terms

    def fetch_articles_metadata(self, session: Session, pmids: list[int]) -> None:
        """Fetch PubMed metadata in batches and store any missing articles."""
        BATCH_SIZE = 250
        for idx, i in enumerate(range(0, len(pmids), BATCH_SIZE), start=1):
            batch = pmids[i : i + BATCH_SIZE]
            logger.info(
                f"Fetching metadata for batch {idx}/{(len(pmids) + BATCH_SIZE - 1) // BATCH_SIZE}"
            )
            try:
                fetch = Entrez.efetch(
                    db="pubmed",
                    id=",".join(map(str, batch)),
                    rettype="medline",
                    retmode="xml",
                )
                data = cast(dict, Entrez.read(fetch, ignore_errors=True))
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
        """Parse MEDLINE XML into PubmedArticle; tolerate missing fields."""
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
                if isinstance(eloc, list):
                    for e in eloc:
                        if getattr(e, "attributes", {}).get("EIdType") == "doi":
                            doi = str(e)
                            break
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
        """Lazy-create SQLAlchemy engine for DuckDB."""
        if self._engine is None:
            self._engine = sqlalchemy.create_engine(f"duckdb:///{self.path}")
        return self._engine
