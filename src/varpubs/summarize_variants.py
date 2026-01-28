import logging
from statistics import mean
from pathlib import Path
from typing import List, Optional, Dict, Tuple

from sqlmodel import Session, select
from cyvcf2 import VCF, Writer

from varpubs.cache import Cache, Judge, Summary
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
    judges: Optional[List[str]] = None,
    output_cache: Optional[Path] = None,
):
    """
    Extracts variant terms from a VCF file, finds related PubMed articles from the database,
    summarizes them using the given summarizer, and optionally saves the summaries to a CSV file.
    """
    if judges is None:
        judges = []
    db = PubmedDB(path=db_path, vcf_paths=[], species=species, max_publications=50)
    engine = db.engine
    cache = summarizer.settings.cache

    with Session(engine) as session:
        vcf = VCF(vcf_path)
        total_record = sum(1 for _ in vcf)
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
        for judge in judges:
            vcf.add_info_to_header(
                {
                    "ID": f"{judge}_score",
                    "Description": f"Varpubs judgement score for {judge}.",
                    "Type": "Float",
                    "Number": ".",
                }
            )
        vcf_out = Writer(out_path, vcf)
        hgvsp_index = get_annotation_field_index(vcf, "HGVSp")
        gene_index = get_annotation_field_index(vcf, "SYMBOL")
        if output_cache:
            ocache = Cache(output_cache)
            ocache.deploy()
        else:
            ocache = None
        for i, record in enumerate(vcf, start=1):
            logging.info(f"Processing vcf record {i}/{total_record}")
            bioconcepts = extract_bioconcept_from_record(
                record, hgvsp_index, gene_index, species
            )
            rec_summaries = []
            rec_pmids = []
            rec_judgements: List[List[Dict[str, int]]] = []
            for bioconcept in bioconcepts:
                judgements: List[Dict] = []
                summaries = {}
                logging.info(f"Summarizing abstracts for: {bioconcept}")
                mappings = session.exec(
                    select(BioconceptToPMID).where(
                        BioconceptToPMID.bioconcept == bioconcept
                    )
                ).all()
                pmids = set(m.pmid for m in mappings)

                for pmid in pmids:
                    article = session.exec(
                        select(PubmedArticle).where(PubmedArticle.pmid == pmid)
                    ).first()
                    if not article:
                        continue

                    if cache:
                        cached_summary = cache.lookup_summary(
                            bioconcept,
                            pmid,
                            summarizer.settings.model,
                            summarizer.summary_prompt_hash(),
                        )
                    elif ocache:
                        cached_summary = ocache.lookup_summary(
                            bioconcept,
                            pmid,
                            summarizer.settings.model,
                            summarizer.summary_prompt_hash(),
                        )
                    else:
                        cached_summary = None

                    hgvsp, gene = bioconcept_to_hgvsp_gene(bioconcept)
                    if cached_summary:
                        summary_text = cached_summary.summary
                    else:
                        logging.info(
                            f"No summary cache entry found for {bioconcept} (pmid: {pmid})"
                        )
                        summary_text = summarizer.summarize_article(
                            article, f"{gene} {hgvsp}"
                        )

                    scores: Dict[str, int] = {}
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
                                ).model_dump()
                            )
                        scores[judge] = score
                    summaries[pmid] = {
                        "article": article,
                        "summary": summary_text,
                        "scores": scores,
                        "term": bioconcept,
                    }
                final_summaries: List[Tuple[PubmedArticle, str]] = [
                    (data["article"], data["summary"]) for data in summaries.values()
                ]
                judge_scores: List[Dict[str, int]] = [
                    data["scores"] for data in summaries.values()
                ]

                hgvs, gene = bioconcept_to_hgvsp_gene(bioconcept)
                summary = (
                    summarizer.summarize(final_summaries, f"{gene} {hgvs}")
                    if final_summaries
                    else "."
                )
                rec_summaries.append(summary.replace(",", "%2C"))
                rec_pmids.append(
                    "|".join(f"{pmid}" for pmid in pmids) if pmids else "."
                )
                rec_judgements.append(judge_scores)

                if ocache:
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
                    ocache.write_judges([Judge(**j) for j in judgements])
            record.INFO["publication_summaries"] = ",".join(rec_summaries)
            record.INFO["PMIDs"] = ",".join(rec_pmids)
            for judge in judges:
                mean_scores = [
                    str(mean(pmid_score[judge] for pmid_score in transcript_scores))
                    if transcript_scores
                    else "."
                    for transcript_scores in rec_judgements
                ]
                record.INFO[f"{judge}_score"] = ",".join(mean_scores)
            vcf_out.write_record(record)
        vcf_out.close()
