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
    summarizes them using the given summarizer, and optionally saves the summaries to a file.
    
    Example terminal usage:
    pixi run varpubs summarize-variants \
    --db_path mydb.duckdb \
    --vcf_path /Users/ilmaytas/Desktop/annotated.vcf \
    --hf_token hf_... \
    --output /Users/ilmaytas/Desktop/summary.txt
    """
    terms = extract_hgvsp_from_vcf(str(vcf_path))  # Extract variant terms from VCF
    db = PubmedDB(path=db_path, vcf_paths=[], email="")
    engine = db.engine

    with Session(engine) as session:
        all_summaries = []

        for term in terms:
            print(f"\n Variant term: {term}")
            mappings = session.exec(
                select(TermToPMID).where(TermToPMID.term == term)
            ).all()
            pmids = set(m.pmid for m in mappings)  # Get unique PMIDs

            for pmid in pmids:
                article = session.exec(
                    select(PubmedArticle).where(PubmedArticle.pmid == pmid)
                ).first()
                if not article:
                    continue

                summary = summarizer.summarize(article)  # Summarize article content
                output = f" Variant: {term}\n PMID: {pmid}\n Summary:\n{summary}\n🔗 DOI: {article.doi or 'N/A'}\n"
                print(output)
                all_summaries.append(output)

        if out_path:
            with open(out_path, "w") as f:
                f.write("\n".join(all_summaries))  # Save all summaries to file
