from pathlib import Path

from cyvcf2 import VCF

from varpubs.hgvs_extractor import extract_hgvsp_from_vcf, get_annotation_field_index


def test_get_annotation_field_index():
    vcf_path = Path("tests/resources/annotated.vcf")
    vcf = VCF(str(vcf_path))
    assert get_annotation_field_index(vcf, "HGVSp") == 10
    assert get_annotation_field_index(vcf, "SYMBOL") == 3


def test_extract_hgvsp_from_vcf():
    vcf_path = Path("tests/resources/annotated.vcf")
    terms = extract_hgvsp_from_vcf(str(vcf_path), "human")
    assert "@VARIANT_p.S183L_MFSD2A_human" in terms
