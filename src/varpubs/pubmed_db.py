from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Iterable
import sqlalchemy
from sqlmodel import Field, SQLModel, Session, Integer
from varpubs.hgvs_extractor import extract_hgvsp_from_vcf
from Bio import Entrez

# New table: PMIDs -> Article metadata
class PubmedArticle(SQLModel, table=True):
    pmid: int = Field(Integer, primary_key=True, nullable=False)
    title: str
    abstract: str
    authors: str
    journal: str
    pub_date: str
    doi: str

# New table: term + pmid mapping
class TermToPMID(SQLModel, table=True):
    term: str = Field(nullable=False)
    pmid: int = Field(foreign_key="pubmedarticle.pmid", nullable=False)

    __table_args__ = (sqlalchemy.PrimaryKeyConstraint("term", "pmid"),)


@dataclass
class PubmedDB:
    path: Path
    vcf_paths: Iterable[Path]
    email: str
    _engine: Optional[sqlalchemy.engine.base.Engine] = field(init=False, default=None)

    def deploy(self) -> None:
        print(f"[INFO] Deploying database at: {self.path}")
        print(f"[INFO] Input VCF paths: {[str(p) for p in self.vcf_paths]}")
        print(f"[INFO] Entrez email: {self.email}")

        if not self.path.exists():
            if self.path.parent != Path():
                self.path.parent.mkdir(parents=True, exist_ok=True)
            print("[INFO] Creating tables...")
            SQLModel.metadata.create_all(self.engine)
            print("[SUCCESS] New database created.")
        else:
            print("[INFO] Database already exists. Will append to it.")

        hgvsp_terms = set()
        for vcf_path in self.vcf_paths:
            hgvsp_terms.update(extract_hgvsp_from_vcf(str(vcf_path)))

        Entrez.email = self.email

        with Session(self.engine) as session:
            for term in hgvsp_terms:
                try:
                    print(f"[INFO] Searching PubMed for: {term}")
                    handle = Entrez.esearch(db="pubmed", term=term, retmax=3)
                    record = Entrez.read(handle)
                    ids = record.get("IdList", [])
                    if not ids:
                        print(f"[INFO] No PubMed results found for HGVS.p term: {term}")
                        continue

                    for pmid in ids:
                        pmid = int(pmid)

                        # Check if (term, pmid) already exists
                        mapping_exists = session.exec(
                            sqlalchemy.select(TermToPMID).where(
                                TermToPMID.term == term,
                                TermToPMID.pmid == pmid
                            )
                        ).first()

                        if mapping_exists:
                            print(f"[INFO] Skipping (term={term}, PMID={pmid}): already mapped.")
                            continue

                        # Check if article already exists
                        article_exists = session.exec(
                            sqlalchemy.select(PubmedArticle).where(
                                PubmedArticle.pmid == pmid
                            )
                        ).first()

                        if not article_exists:
                            fetch = Entrez.efetch(db="pubmed", id=str(pmid), rettype="medline", retmode="xml")
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

                            journal = article.get("Journal", {}).get("Title", "")
                            pub_date = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {}).get("Year", "")
                            doi = ""
                            try:
                                eloc = article.get("ELocationID", [{}])
                                if isinstance(eloc, list) and eloc and isinstance(eloc[0], dict):
                                    doi = eloc[0].get("#text", "")
                            except:
                                pass

                            session.add(PubmedArticle(
                                pmid=pmid,
                                title=title,
                                abstract=abstract,
                                authors=authors,
                                journal=journal,
                                pub_date=pub_date,
                                doi=doi
                            ))

                        session.add(TermToPMID(term=term, pmid=pmid))

                except Exception as e:
                    print(f"[ERROR] Skipping term {term} due to error: {e}")

            session.commit()

    @property
    def engine(self) -> sqlalchemy.engine.base.Engine:
        if self._engine is None:
            self._engine = sqlalchemy.create_engine(f"duckdb:///{self.path}")
        return self._engine
