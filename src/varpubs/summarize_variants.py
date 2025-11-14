import csv
import logging
from pathlib import Path
from typing import List, Optional

from sqlmodel import Session, select

from varpubs.cache import Cache, Judge, Summary
from varpubs.hgvs_extractor import extract_hgvsp_from_vcf
from varpubs.pubmed_db import PubmedArticle, PubmedDB, TermToPMID
from varpubs.summarize import PubmedSummarizer


def summarize_variants(
    db_path: Path,
    vcf_path: Path,
    summarizer: PubmedSummarizer,
    out_path: Optional[Path] = None,
    judges: Optional[list[str]] = None,
    output_cache: Optional[Path] = None,
):
    """
    Extracts variant terms from a VCF file, finds related PubMed articles from the database,
    summarizes them using the given summarizer, and optionally saves the summaries to a CSV file.
    """

    terms = extract_hgvsp_from_vcf(str(vcf_path))
    db = PubmedDB(path=db_path, vcf_paths=[], email="")
    engine = db.engine
    cache = summarizer.settings.cache
    judgements: List[Judge] = []

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

                summary_text = (
                    cache.lookup_summary(
                        term,
                        pmid,
                        summarizer.settings.model,
                        summarizer.summary_prompt_hash(),
                    )
                    if cache
                    else None
                )

                if not summary_text:
                    summary_text = summarizer.summarize_article(article, term)

                scores = {}
                if not judges:
                    judges = []
                for judge in judges:
                    score = (
                        cache.lookup_judge(
                            term,
                            pmid,
                            summarizer.settings.model,
                            judge,
                            summarizer.judge_prompt_hash(),
                        )
                        if cache
                        else None
                    )
                    if not score:
                        score = summarizer.judge(article, judge)
                        judgements.append(
                            Judge(
                                term=term,
                                pmid=pmid,
                                model=summarizer.settings.model,
                                judge=judge,
                                score=score,
                                prompt_hash=summarizer.judge_prompt_hash(),
                            )
                        )
                    scores[judge] = score
                summaries[pmid] = {
                    "article": article,
                    "summary": summary_text,
                    "scores": scores,
                    "term": term,
                }
            sorted_summaries = sorted(
                summaries.items(),
                key=lambda x: sum(x[1]["scores"].values()) if x[1]["scores"] else 0,
                reverse=True,
            )[:50]
            top_summaries: list[tuple[PubmedArticle, str]] = [
                (data["article"], data["summary"]) for pmid, data in sorted_summaries
            ]

            if top_summaries:
                summary = summarizer.summarize(top_summaries, term)
                gene, symbol = term.split(" ")
                rows.append(
                    tuple(
                        [gene, symbol, summary, ",".join(f"{pmid}" for pmid in pmids)]
                    )
                )

            if output_cache:
                ocache = Cache(output_cache)
                ocache.deploy()
                s: List[Summary] = [
                    Summary(
                        term=data["term"],
                        pmid=pmid,
                        model=summarizer.settings.model,
                        summary=data["summary"],
                        prompt_hash=summarizer.summary_prompt_hash(),
                    )
                    for pmid, data in summaries.items()
                ]
                ocache.write_summaries(s)
                ocache.write_judges(judgements)

        if out_path:
            with open(out_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(["symbol", "hgvsp", "summary", "pmids"])
                writer.writerows(rows)
