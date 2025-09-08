import logging
import csv
from pathlib import Path
from typing import Optional
from sqlmodel import Session, select
from varpubs.hgvs_extractor import extract_hgvsp_from_vcf
from varpubs.pubmed_db import PubmedArticle, TermToPMID, PubmedDB
from varpubs.summarize import PubmedSummarizer


def summarize_variants(
    db_path: Path,
    vcf_path: Path,
    summarizer: PubmedSummarizer,
    out_path: Optional[Path] = None,
):
    """
    Extracts variant terms from a VCF file, finds related PubMed articles from the database,
    summarizes them using the given summarizer, and optionally saves the summaries to a CSV file.
    """

    terms = extract_hgvsp_from_vcf(str(vcf_path))
    db = PubmedDB(path=db_path, vcf_paths=[], email="")
    engine = db.engine

    with Session(engine) as session:
        rows = []

        for term in terms:
            logging.info(f"Summarizing abstracts for variant term: {term}")
            mappings = session.exec(
                select(TermToPMID).where(TermToPMID.term == term)
            ).all()
            pmids = set(m.pmid for m in mappings)
            summaries = []
            for pmid in pmids:
                article = session.exec(
                    select(PubmedArticle).where(PubmedArticle.pmid == pmid)
                ).first()
                if not article:
                    continue

                summary_text = summarizer.summarize_article(article, term)
                summaries.append(tuple([article, summary_text]))
            if summaries:
                summary = summarizer.summarize(summaries, term)
                rows.append(
                    tuple([term, ",".join(f"{pmid}" for pmid in pmids), summary])
                )

        if out_path:
            with open(out_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(["Variant", "PMIDs", "Summary"])
                writer.writerows(rows)
