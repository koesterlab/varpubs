from pathlib import Path

import pytest
from cyvcf2 import VCF

from varpubs.hgvs_extractor import (
    extract_bioconcepts_from_table,
    extract_hgvsp_from_vcf,
    get_annotation_field_index,
    table_row_to_bioconcept,
)


def test_get_annotation_field_index():
    vcf_path = Path("tests/resources/annotated.vcf")
    vcf = VCF(str(vcf_path))
    assert get_annotation_field_index(vcf, "HGVSp") == 10
    assert get_annotation_field_index(vcf, "SYMBOL") == 3


def test_extract_hgvsp_from_vcf():
    vcf_path = Path("tests/resources/annotated.vcf")
    terms = extract_hgvsp_from_vcf(str(vcf_path), "human")
    assert "@VARIANT_p.S183L_MFSD2A_human" in terms


def test_extract_bioconcepts_from_table(tmp_path):
    table = tmp_path / "variants.tsv"
    table.write_text("gene\tvariant\nNF1\tp.R1748*\nBRAF\tp.V600E\n")
    bioconcepts = extract_bioconcepts_from_table(str(table), "gene", "variant", "human")
    assert bioconcepts == {
        "@VARIANT_p.R1748*_NF1_human",
        "@VARIANT_p.V600E_BRAF_human",
    }


def test_table_row_to_bioconcept_strips_whitespace():
    # spreadsheet exports carry trailing spaces and non-breaking spaces
    assert (
        table_row_to_bioconcept("KRAS ", "p.G12S\xa0", "human")
        == "@VARIANT_p.G12S_KRAS_human"
    )


def test_extract_bioconcepts_from_table_missing_column(tmp_path):
    table = tmp_path / "variants.tsv"
    table.write_text("gene\tprotein\nNF1\tp.R1748*\n")
    with pytest.raises(ValueError, match="not found"):
        extract_bioconcepts_from_table(str(table), "gene", "variant", "human")
