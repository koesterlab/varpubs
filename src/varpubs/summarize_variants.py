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
    judges: Optional[list[str]] = None,
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
            summaries = {}
            for pmid in pmids:
                article = session.exec(
                    select(PubmedArticle).where(PubmedArticle.pmid == pmid)
                ).first()
                if not article:
                    continue

                summary_text = summarizer.summarize_article(article, term)
                scores = {}
                if not judges:
                    judges = []
                for judge in judges:
                    scores[judge] = summarizer.judge(article, judge)
                summaries[pmid] = {"summary": summary_text, "scores": scores}
            top_summaries = sorted(
                summaries.items(),
                key=lambda x: sum(x[1]["scores"].values()) if x[1]["scores"] else 0,
                reverse=True,
            )[:50]
            top_summaries = [(pmid, data["summary"]) for pmid, data in top_summaries]

            if top_summaries:
                summary = summarizer.summarize(top_summaries, term)
                gene, symbol = term.split(" ")
                rows.append(
                    tuple(
                        [gene, symbol, summary, ",".join(f"{pmid}" for pmid in pmids)]
                    )
                )

        if out_path:
            with open(out_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(["symbol", "hgvsp", "summary", "pmids"])
                writer.writerows(rows)
