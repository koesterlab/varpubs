# -----------------------------------------------------------------------------
# This module defines a PubmedDB class that:
# - Parses a VCF file to extract HGVS.p terms.
# - Queries PubMed via the Entrez API for abstracts mentioning each term.
# - Stores the retrieved article metadata into a DuckDB database using SQLModel.
# The PubmedArticle table uses a composite primary key: (term, pmid).
# -----------------------------------------------------------------------------

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import sqlalchemy
from sqlmodel import Field, SQLModel, Session
from varpubs.hgvs_extractor import extract_hgvsp_from_vcf
from Bio import Entrez

# Table schema for storing PubMed articles related to each HGVS.p term.
# Each (term, pmid) pair is unique.
class PubmedArticle(SQLModel, table=True):
    term: str = Field(nullable=False)
    pmid: int = Field(nullable=False)
    title: str
    abstract: str
    authors: str
    journal: str
    pub_date: str
    doi: str

    __table_args__ = (sqlalchemy.PrimaryKeyConstraint("term", "pmid"),)

# Class responsible for deploying and populating the PubMed literature database.
# Uses DuckDB for storage and Entrez API for fetching article metadata.
@dataclass
class PubmedDB:
    path: Path
    vcf_path: Path
    email: str
    _engine: Optional[sqlalchemy.engine.base.Engine] = field(init=False, default=None)

    # Main deployment method:
    # - Initializes the DuckDB database.
    # - Extracts HGVS.p terms from the input VCF file.
    # - Queries PubMed for each term.
    # - Inserts new article entries into the database (if not already present).
    def deploy(self) -> None:
        print(f"[INFO] Deploying database at: {self.path}")
        print(f"[INFO] VCF path: {self.vcf_path}")
        print(f"[INFO] Entrez email: {self.email}")

        if self.path.exists():
            print("[INFO] Removing existing database file.")
            self.path.unlink()

        if self.path.parent != Path():
            self.path.parent.mkdir(parents=True, exist_ok=True)

        print("[INFO] Creating tables...")
        SQLModel.metadata.create_all(self.engine)
        print("[SUCCESS] Database deployed.")

        hgvsp_terms = extract_hgvsp_from_vcf(str(self.vcf_path))
        Entrez.email = self.email

        with Session(self.engine) as session:
            for term in set(hgvsp_terms):
                try:
                    print(f"[INFO] Searching PubMed for: {term}")
                    handle = Entrez.esearch(db="pubmed", term=term, retmax=3)
                    record = Entrez.read(handle)
                    ids = record.get("IdList", [])
                    if not ids:
                        print(f"[INFO] No PubMed results found for HGVS.p term: {term}")

                    for pmid in ids:
                        fetch = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="xml")
                        article_data = Entrez.read(fetch)

                        try:
                            pubmed_articles = article_data.get("PubmedArticle", [])
                            if not pubmed_articles:
                                print(f"[WARN] No PubmedArticle found for PMID {pmid}")
                                continue

                            article = pubmed_articles[0].get("MedlineCitation", {}).get("Article", {})
                        except Exception as e:
                            print(f"[ERROR] Failed to extract article for PMID {pmid}: {e}")
                            continue

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
                            print(f"[WARN] Abstract parse failed for PMID {pmid}: {e}")

                        authors = ""
                        try:
                            author_list = article.get("AuthorList", [])
                            if isinstance(author_list, list):
                                authors = ", ".join(
                                    f"{a.get('LastName', '')} {a.get('ForeName', '')}".strip()
                                    for a in author_list if isinstance(a, dict) and "LastName" in a
                                )
                        except Exception as e:
                            print(f"[WARN] Author parse failed for PMID {pmid}: {e}")

                        journal = ""
                        try:
                            journal = article.get("Journal", {}).get("Title", "")
                        except:
                            pass

                        pub_date = ""
                        try:
                            pub_date = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {}).get("Year", "")
                        except:
                            pass

                        doi = ""
                        try:
                            eloc = article.get("ELocationID", [{}])
                            if isinstance(eloc, list) and eloc and isinstance(eloc[0], dict):
                                doi = eloc[0].get("#text", "")
                        except:
                            pass

                        existing = session.exec(
                            sqlalchemy.select(PubmedArticle).where(
                                PubmedArticle.term == term,
                                PubmedArticle.pmid == int(pmid)
                            )
                        ).first()

                        if existing:
                            print(f"[INFO] Skipping (term={term}, PMID={pmid}): already in database.")
                            continue

                        session.add(PubmedArticle(
                            term=term,
                            pmid=int(pmid),
                            title=title,
                            abstract=abstract,
                            authors=authors,
                            journal=journal,
                            pub_date=pub_date,
                            doi=doi
                        ))
                except Exception as e:
                    print(f"[ERROR] Skipping term {term} or PMID {pmid} due to error: {e}")
            session.commit()

    # Lazy initialization of the DuckDB engine using the provided database path.
    @property
    def engine(self) -> sqlalchemy.engine.base.Engine:
        if self._engine is None:
            self._engine = sqlalchemy.create_engine(f"duckdb:///{self.path}")
        return self._engine
