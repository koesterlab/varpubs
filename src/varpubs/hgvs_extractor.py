import logging
from cyvcf2 import VCF
from typing import Tuple, List, Any
from hgvs.parser import Parser

logger = logging.getLogger(__name__)


def get_annotation_field_index(vcf: VCF, field: str) -> int:
    for rec in vcf.header_iter():
        info = rec.info()
        if info.get("ID") == "ANN":
            desc = info.get("Description", "")
            if field in desc:
                fields_str = desc.split("': '")[-1].rstrip("'")
                fields = [f.strip() for f in fields_str.split("|")]
                try:
                    return fields.index(field)
                except ValueError as e:
                    raise RuntimeError(f"{field} not found in ANN header") from e
    raise RuntimeError("ANN field not found in VCF header")


def extract_hgvsp_from_vcf(vcf_path: str, species: str) -> set[str]:
    vcf = VCF(vcf_path)
    hgvsp_index = get_annotation_field_index(vcf, "HGVSp")
    gene_index = get_annotation_field_index(vcf, "SYMBOL")
    terms: List[str] = []

    for record in vcf:
        terms.extend(
            extract_bioconcept_from_record(record, hgvsp_index, gene_index, species)
        )
    return set(terms)


def extract_bioconcept_from_record(
    record: Any, hgvsp_index: int, gene_index: int, species: str
) -> List[str]:
    ann = record.INFO.get("ANN")
    bioconcepts = []
    hgvs_parser = Parser()
    if ann:
        for ann_entry in ann.split(","):
            fields = ann_entry.split("|")
            # TODO: Consider filtering only canonical transcripts
            if len(fields) > max(hgvsp_index, gene_index):
                if not fields[hgvsp_index]:
                    logger.warning(f"HGVSp entry is empty: {ann_entry}")
                    continue
                hgvsp = fields[hgvsp_index]
                # skip synonymous variants
                if "%3D" in hgvsp:
                    continue
                hgvsp_single = (
                    hgvs_parser.parse(hgvsp)
                    .format(conf={"p_3_letter": False})
                    .split(":")[1]
                )
                gene = fields[gene_index]

                # Create bioconcept for querying pubtator
                bioconcepts.append(
                    hgvsp_gene_to_bioconcept(hgvsp_single, gene, species)
                )
    return bioconcepts


def hgvsp_gene_to_bioconcept(hgvsp: str, gene: str, species: str) -> str:
    return f"@VARIANT_{hgvsp}_{gene}_{species}"


def bioconcept_to_hgvsp_gene(bioconcept: str) -> Tuple[str, str]:
    hgvsp, gene = bioconcept.split("_")[1:3]
    return hgvsp, gene
