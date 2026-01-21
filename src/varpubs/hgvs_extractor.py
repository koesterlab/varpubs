import logging
from cyvcf2 import VCF
from typing import Tuple, List, Any

logger = logging.getLogger(__name__)
# 3-letter to 1-letter amino acid codes
AA3_TO_1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Glu": "E",
    "Gln": "Q",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
    "*": "*",
    "del": "del",
}

AA1_TO_3 = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "E": "Glu",
    "Q": "Gln",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "*": "Ter",
}


def to_3_letter(term: str) -> str:
    return "".join(AA1_TO_3.get(aa, aa) for aa in term)


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
    term_set: set[str] = set()

    for record in vcf:
        term_set.union(set(extract_bioconcept_from_record(record, hgvsp_index, gene_index, species)))
    return term_set

def extract_bioconcept_from_record(record: Any, hgvsp_index: int, gene_index: int, species: str) -> List[str]:
    ann = record.INFO.get("ANN")
    bioconcepts = []
    if ann:
        for ann_entry in ann.split(","):
            fields = ann_entry.split("|")
            # TODO: Consider filtering only canonical transcripts
            if len(fields) > max(hgvsp_index, gene_index):
                if not fields[hgvsp_index]:
                    logger.warning(f"HGVSp entry is empty: {ann_entry}")
                    continue
                hgvsp = fields[hgvsp_index].split(":")[1]
                gene = fields[gene_index]
                if not hgvsp.startswith("p."):
                    logger.warning(
                        f"HGVSp entry does not seem to be valid: {hgvsp}"
                    )
                    continue

                # skip synonymous variants
                if "%3D" in hgvsp:
                    continue

                for long, short in AA3_TO_1.items():
                    hgvsp = hgvsp.replace(long, short)

                # Create bioconcept for querying pubtator
                bioconcepts.append(hgvsp_gene_to_bioconcept(hgvsp, gene, species))
    return bioconcepts


def hgvsp_gene_to_bioconcept(hgvsp: str, gene: str, species: str) -> str:
    return f"@VARIANT_{hgvsp}_{gene}_{species}"


def bioconcept_to_hgvsp_gene(bioconcept: str) -> Tuple[str, str]:
    hgvsp, gene = bioconcept.split("_")[1:3]
    return hgvsp, gene
