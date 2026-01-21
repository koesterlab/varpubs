import logging
from pathlib import Path
from typing import List, Optional

from sqlmodel import Session, select
from cyvcf2 import VCF, Writer

from varpubs.cache import Cache, Judge, Summary
from varpubs.hgvs_extractor import extract_hgvsp_from_vcf
from varpubs.pubmed_db import PubmedArticle, PubmedDB, BioconceptToPMID
from varpubs.summarize import PubmedSummarizer
from varpubs.hgvs_extractor import (
    bioconcept_to_hgvsp_gene,
    get_annotation_field_index,
    extract_bioconcept_from_record,
)


def summarize_variants(
    db_path: Path,
    vcf_path: Path,
    summarizer: PubmedSummarizer,
    species: str,
    out_path: Optional[Path] = None,
    judges: Optional[list[str]] = None,
    output_cache: Optional[Path] = None,
):
    """
    Extracts variant terms from a VCF file, finds related PubMed articles from the database,
    summarizes them using the given summarizer, and optionally saves the summaries to a CSV file.
    """

    bioconcepts = extract_hgvsp_from_vcf(str(vcf_path), species)
    db = PubmedDB(path=db_path, vcf_paths=[], species=species, max_publications=50)
    engine = db.engine
    cache = summarizer.settings.cache
    judgements: List[Judge] = []

    with Session(engine) as session:
        vcf = VCF(vcf_path)
        vcf.add_info_to_header(
            {
                "ID": "publication_summaries",
                "Description": "Summary of related PubMed articles for each transcript.",
                "Type": "String",
                "Number": ".",
            }
        )
        vcf.add_info_to_header(
            {
                "ID": "PMIDs",
                "Description": "PubMed IDs of related articles for each transcript.",
                "Type": "String",
                "Number": ".",
            }
        )
        if judges:
            for judge in judges:
                vcf.add_info_to_header(
                    {
                        "ID": f"{judge}_score",
                        "Description": f"Varpubs judgement score for {judge}.",
                        "Type": "Number",
                        "Number": ".",
                    }
                )
        vcf_out = Writer(out_path, vcf)
        hgvsp_index = get_annotation_field_index(vcf, "HGVSp")
        gene_index = get_annotation_field_index(vcf, "SYMBOL")
        for record in vcf:
            bioconcepts = extract_bioconcept_from_record(
                record, hgvsp_index, gene_index, species
            )
            rec_summaries = []
            rec_pmids = []
            rec_judgements = []
            for bioconcept in bioconcepts:
                logging.info(f"Summarizing abstracts for: {bioconcept}")
                mappings = session.exec(
                    select(BioconceptToPMID).where(
                        BioconceptToPMID.bioconcept == bioconcept
                    )
                ).all()
                pmids = set(m.pmid for m in mappings)
                summaries = {}
                for pmid in pmids:
                    article = session.exec(
                        select(PubmedArticle).where(PubmedArticle.pmid == pmid)
                    ).first()
                    if not article:
                        continue

                    cached_summary = (
                        cache.lookup_summary(
                            bioconcept,
                            pmid,
                            summarizer.settings.model,
                            summarizer.summary_prompt_hash(),
                        )
                        if cache
                        else None
                    )
                    hgvsp, gene = bioconcept_to_hgvsp_gene(bioconcept)
                    summary_text = (
                        cached_summary.summary
                        if cached_summary
                        else summarizer.summarize_article(article, f"{gene} {hgvsp}")
                    )

                    scores: dict[str, int] = {}
                    if not judges:
                        judges = []
                    for judge in judges:
                        score = (
                            cache.lookup_judge(
                                bioconcept,
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
                                    term=bioconcept,
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
                        "term": bioconcept,
                    }
                final_summaries: list[tuple[PubmedArticle, str]] = [
                    (data["article"], data["summary"]) for data in summaries.values()
                ]
                judge_scores: List[dict[str, int]] = [
                    data["scores"] for data in summaries.values()
                ]

                hgvs, gene = bioconcept_to_hgvsp_gene(bioconcept)
                summary = summarizer.summarize(final_summaries, f"{gene} {hgvs}")
                rec_summaries.append(summary.replace(",", "%2C"))
                rec_pmids.append("|".join(f"{pmid}" for pmid in pmids))
                rec_judgements.append(judge_scores)

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

            record.INFO["publication_summaries"] = ",".join(rec_summaries)
            record.INFO["PMIDs"] = ",".join(rec_pmids)
            if judges:
                for judge in judges:
                    record.INFO[f"{judge}_score"] = ",".join(
                        score for score in judge_scores[judge]
                    )
            vcf_out.write_record(record)
        vcf_out.close()
        # if out_path:
        #     with open(out_path, "w", newline="", encoding="utf-8") as f:
        #         writer = csv.writer(f)
        #         writer.writerow(["symbol", "hgvsp", "summary", "pmids"])
        #         writer.writerows(rows)
