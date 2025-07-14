import pytest
from cyvcf2 import VCF
from pathlib import Path

from varpubs.hgvs_extractor import get_annotation_field_index, extract_hgvsp_from_vcf


def test_get_annotation_field_index():
    vcf_path = Path("tests/resources/annotated.vcf")
    vcf = VCF(str(vcf_path))
    assert get_annotation_field_index(vcf, "HGVS.p") == 10
    assert get_annotation_field_index(vcf, "Gene_Name") == 3


def test_extract_hgvsp_from_vcf():
    vcf_path = Path("tests/resources/annotated.vcf")
    terms = extract_hgvsp_from_vcf(str(vcf_path))
    assert "ALAD p.Lys59Asn" in terms
    assert "ALAD G59" in terms
