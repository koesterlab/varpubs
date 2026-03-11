from cyvcf2 import VCF

from varpubs.utils import extend_vep_header


def test_extend_vep_header_ann():
    vcf = VCF("tests/resources/annotated.vcf")
    new_fields = ["VARPUBS_SUMMARY", "VARPUBS_PMIDS", "VARPUBS_SCORE"]
    extend_vep_header(vcf, new_fields, csq_id="ANN")
    description = vcf.get_header_type("ANN")["Description"]
    for field in new_fields:
        assert field in description
