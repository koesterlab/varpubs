import csv
from pathlib import Path

from cyvcf2 import VCF
from sqlmodel import Session

from varpubs.pubmed_db import BioconceptToPMID, PubmedArticle, PubmedDB
from varpubs.summarize import PubmedSummarizer, Settings
from varpubs.summarize_variants import (
    process_bioconcept,
    summarize_variants,
    summarize_variants_table,
)

JUDGES = ["therapy related", "oncogenicity"]


class StubSummarizer(PubmedSummarizer):
    """Deterministic summarizer that returns fixed text and scores without a network call."""

    def summarize_article(self, article, term):
        return f"summary of {article.pmid}"

    def judge(self, article, term):
        return 4

    def summarize(self, texts, term, judge):
        return f"{judge} findings"


def test_summarize_variants_annotates_every_transcript(tmp_path):
    """Against an empty database no article is found, so the pipeline runs end to
    end without an LLM call and still appends the VARPUBS fields to every
    transcript of every record."""
    db_path = tmp_path / "empty.duckdb"
    PubmedDB(
        path=db_path, vcf_paths=[], species="human", max_publications=50
    ).create_tables()

    out_path = tmp_path / "out.vcf"
    summarizer = PubmedSummarizer(Settings(api_key="", role="physician"))
    summarize_variants(
        db_path=db_path,
        vcf_path=Path("tests/resources/annotated.vcf"),
        summarizer=summarizer,
        species="human",
        judges=JUDGES,
        out_path=out_path,
    )

    out = VCF(str(out_path))
    ann_desc = out.get_header_type("ANN")["Description"]
    for field in [
        "VARPUBS_SUMMARY",
        "VARPUBS_PMIDS",
        *[f"VARPUBS_{j}_SCORE" for j in JUDGES],
    ]:
        assert field in ann_desc

    records = list(out)
    assert records
    for record in records:
        for transcript in record.INFO["ANN"].split(","):
            # empty summary, empty pmids and one empty score per judge are appended
            assert transcript.endswith("|" * (2 + len(JUDGES)))


def test_process_bioconcept_summarizes_and_judges(tmp_path):
    """With a mapped article, the bioconcept is summarized and judged and the
    resulting record carries the pmids, per-judge scores and summary text."""
    db = PubmedDB(
        path=tmp_path / "db.duckdb", vcf_paths=[], species="human", max_publications=50
    )
    db.create_tables()
    bioconcept = "@VARIANT_p.R1748*_NF1_human"
    with Session(db.engine) as session:
        session.add(
            PubmedArticle(
                pmid=1, title="t", abstract="a", authors="x", journal="j", doi="d"
            )
        )
        session.add(BioconceptToPMID(bioconcept=bioconcept, pmid=1))
        session.commit()

        summarizer = StubSummarizer(Settings(api_key="", role="physician"))
        record = process_bioconcept(
            bioconcept, session, summarizer, ["therapy related"], None
        )

    assert record.pmids == {1}
    assert record.judges == [{"therapy related": 4}]
    assert record.mean_score("therapy related") == "4"
    assert record.summary.startswith("therapy related:")


def test_summarize_variants_table(tmp_path):
    """A gene + variant table is annotated in place: input columns are preserved
    and the varpubs columns are appended, with mapped rows summarized and unmapped
    rows left empty."""
    db = PubmedDB(
        path=tmp_path / "db.duckdb", vcf_paths=[], species="human", max_publications=50
    )
    db.create_tables()
    with Session(db.engine) as session:
        session.add(
            PubmedArticle(
                pmid=1, title="t", abstract="a", authors="x", journal="j", doi="d"
            )
        )
        session.add(BioconceptToPMID(bioconcept="@VARIANT_p.R1748*_NF1_human", pmid=1))
        session.commit()
    db.engine.dispose()  # release the DuckDB file before the pipeline reopens it

    table = tmp_path / "in.tsv"
    table.write_text(
        "gene\tvariant\tnote\nNF1\tp.R1748*\tkeep me\nTP53\tp.R175H\tno lit\n"
    )
    out_path = tmp_path / "out.tsv"
    summarize_variants_table(
        db_path=tmp_path / "db.duckdb",
        table_path=table,
        summarizer=StubSummarizer(Settings(api_key="", role="physician")),
        species="human",
        judges=["therapy related"],
        gene_column="gene",
        variant_column="variant",
        out_path=out_path,
    )

    with open(out_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    assert [r["note"] for r in rows] == ["keep me", "no lit"]  # passthrough preserved
    mapped, unmapped = rows
    assert mapped["varpubs_pmids"] == "1"
    assert mapped["varpubs_therapy related_score"] == "4"
    assert mapped["varpubs_summary"].startswith("therapy related:")
    assert "\n" in mapped["varpubs_summary"]  # multi-line cell round-trips
    assert unmapped["varpubs_pmids"] == ""
    assert unmapped["varpubs_summary"] == ""
