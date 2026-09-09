import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import Dict, List, Optional, Set, Tuple

from cyvcf2 import VCF, Writer
from sqlalchemy.engine import Engine
from sqlmodel import Session, select

from varpubs.cache import Cache, Judge, Summary
from varpubs.hgvs_extractor import (
    bioconcept_to_hgvsp_gene,
    extract_bioconcept_from_record,
    get_annotation_field_index,
    table_delimiter,
    table_row_to_bioconcept,
)
from varpubs.pubmed_db import BioconceptToPMID, PubmedArticle, PubmedDB
from varpubs.summarize import PubmedSummarizer
from varpubs.utils import extend_vep_header


@dataclass
class TranscriptRecord:
    pmids: Set[int]
    summary: str
    judges: List[Dict[str, int]]

    def mean_score(self, judge: str) -> str:
        if self.judges:
            return str(mean(pmid_score[judge] for pmid_score in self.judges))
        else:
            return ""

    def join_pmids(self) -> str:
        if self.pmids:
            return "&".join(str(pmid) for pmid in list(self.pmids))
        else:
            return ""


def process_bioconcept(
    bioconcept: str,
    session: Session,
    summarizer: PubmedSummarizer,
    judges: List[str],
    ocache: Optional[Cache],
) -> TranscriptRecord:
    """
    Finds PubMed articles related to a single bioconcept, summarizes and judges them,
    and returns the resulting per-transcript record. Independent of the input format.
    """
    cache = summarizer.settings.cache
    judgements: List[Dict] = []
    summaries = {}
    mappings = session.exec(
        select(BioconceptToPMID).where(BioconceptToPMID.bioconcept == bioconcept)
    ).all()
    # Skip synonymous variants ("=" not in bioconcept) and variants without an hgvsp annotation while still creating a TranscriptRecord per transcript
    # VARIANT__ checks whether f"@VARIANT_{hgvsp}..." actually contains an hgvsp value
    pmids = (
        set(m.pmid for m in mappings)
        if "=" not in bioconcept and "VARIANT__" not in bioconcept
        else set()
    )
    if pmids:
        logging.info(f"Summarizing abstracts for: {bioconcept}")

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
            summary_text = summarizer.summarize_article(article, f"{gene} {hgvsp}")

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
                if score:
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
            scores[judge] = score or 1
        summaries[pmid] = {
            "article": article,
            "summary": summary_text,
            "scores": scores,
            "term": bioconcept,
        }
    judge_scores: List[Dict[str, int]] = [data["scores"] for data in summaries.values()]

    hgvs, gene = bioconcept_to_hgvsp_gene(bioconcept)
    final_summary = ""
    for judge in judges:
        relevant_summaries = [
            (data["article"], data["summary"])
            for data in summaries.values()
            if data["scores"].get(judge) > 1
        ]
        if not relevant_summaries:
            continue
        judge_term_summary = summarizer.summarize(
            relevant_summaries, f"{gene} {hgvs}", judge
        )
        final_summary += f"{judge}:\n\n{judge_term_summary}\n\n"

    transcript_record = TranscriptRecord(
        pmids=pmids,
        summary=final_summary,
        judges=judge_scores,
    )
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
            if data["summary"]
        ]
        ocache.write_summaries(s)
        ocache.write_judges([Judge(**j) for j in judgements])
    return transcript_record


def _prepare(
    db_path: Path,
    species: str,
    judges: List[str],
    output_cache: Optional[Path],
) -> Tuple[Engine, Optional[Cache]]:
    """Open the article database engine and (optionally) the output cache."""
    if not judges:
        raise ValueError("At least one judge must be specified for summarization.")
    engine = PubmedDB(
        path=db_path, vcf_paths=[], species=species, max_publications=50
    ).engine
    ocache = Cache(output_cache) if output_cache else None
    if ocache:
        ocache.deploy()
    return engine, ocache


def summarize_variants(
    db_path: Path,
    vcf_path: Path,
    summarizer: PubmedSummarizer,
    species: str,
    judges: List[str],
    out_path: Optional[Path] = None,
    output_cache: Optional[Path] = None,
):
    """
    Extracts variant terms from a VCF file, finds related PubMed articles from the database,
    summarizes them using the given summarizer, and writes the annotated VCF.
    """
    engine, ocache = _prepare(db_path, species, judges, output_cache)

    with Session(engine) as session:
        vcf = VCF(vcf_path)
        total_record = sum(1 for _ in vcf)
        vcf = VCF(vcf_path)
        extend_vep_header(
            vcf,
            [
                "VARPUBS_SUMMARY",
                "VARPUBS_PMIDS",
                *[f"VARPUBS_{j}_SCORE" for j in judges],
            ],
            "ANN",
        )
        vcf_out = Writer(out_path, vcf)
        hgvsp_index = get_annotation_field_index(vcf, "HGVSp")
        gene_index = get_annotation_field_index(vcf, "SYMBOL")
        for i, record in enumerate(vcf, start=1):
            logging.info(f"Processing vcf record {i}/{total_record}")
            bioconcepts = extract_bioconcept_from_record(
                record, hgvsp_index, gene_index, species
            )
            transcript_records: Dict[str, TranscriptRecord] = {}
            for bioconcept in bioconcepts:
                if bioconcept not in transcript_records:
                    transcript_records[bioconcept] = process_bioconcept(
                        bioconcept, session, summarizer, judges, ocache
                    )

            transcript_infos: List[TranscriptRecord] = [
                transcript_records[bioconcept] for bioconcept in bioconcepts
            ]
            # get record info ANN/CSQ, split on transcripts, then zip together with our TranscriptRecord list and append in correct order according
            # to extend_vep_header call earlier
            transcript_annotations: List[str] = record.INFO["ANN"].split(",")
            ann: List[str] = []
            for transcript_info, ann_str in zip(
                transcript_infos, transcript_annotations
            ):
                # Commas separate transcripts in the ANN field, so encode them
                summary = transcript_info.summary.replace(",", "%2C")
                transcript_annotation = (
                    f"{ann_str}|{summary}|{transcript_info.join_pmids()}"
                )
                for judge in judges:
                    transcript_annotation = (
                        f"{transcript_annotation}|{transcript_info.mean_score(judge)}"
                    )
                ann.append(transcript_annotation)

            record.INFO["ANN"] = ",".join(ann)
            vcf_out.write_record(record)
        vcf_out.close()


def summarize_variants_table(
    db_path: Path,
    table_path: Path,
    summarizer: PubmedSummarizer,
    species: str,
    judges: List[str],
    gene_column: str,
    variant_column: str,
    out_path: Path,
    output_cache: Optional[Path] = None,
    delimiter: Optional[str] = None,
):
    """
    Reads a TSV/CSV of gene + variant rows, summarizes related PubMed articles per
    row, and writes the input table back with appended varpubs columns.
    """
    engine, ocache = _prepare(db_path, species, judges, output_cache)
    delimiter = delimiter or table_delimiter(table_path)
    score_columns = [f"varpubs_{judge}_score" for judge in judges]
    new_columns = ["varpubs_summary", "varpubs_pmids", *score_columns]

    records: Dict[str, TranscriptRecord] = {}
    matched = 0
    total = 0
    with (
        Session(engine) as session,
        open(table_path, newline="") as infile,
        open(out_path, "w", newline="") as outfile,
    ):
        reader = csv.DictReader(infile, delimiter=delimiter)
        columns = reader.fieldnames or []
        missing = [c for c in (gene_column, variant_column) if c not in columns]
        if missing:
            raise ValueError(
                f"Column(s) {missing} not found in {table_path}. Available: {columns}"
            )
        clashing = [c for c in new_columns if c in columns]
        if clashing:
            raise ValueError(
                f"Input table already contains output column(s) {clashing}"
            )

        writer = csv.DictWriter(
            outfile,
            fieldnames=[*columns, *new_columns],
            delimiter=delimiter,
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            total += 1
            bioconcept = table_row_to_bioconcept(
                row[gene_column], row[variant_column], species
            )
            if bioconcept not in records:
                records[bioconcept] = process_bioconcept(
                    bioconcept, session, summarizer, judges, ocache
                )
            transcript = records[bioconcept]
            if transcript.pmids:
                matched += 1
            row["varpubs_summary"] = transcript.summary
            row["varpubs_pmids"] = transcript.join_pmids()
            for judge, column in zip(judges, score_columns):
                row[column] = transcript.mean_score(judge)
            writer.writerow(row)

    logging.info(f"{matched}/{total} table rows matched the database")
    if total and not matched:
        logging.warning(
            "No table rows matched the database. Did you run deploy-db over this table?"
        )
