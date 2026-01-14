import logging
from cyvcf2 import VCF

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
        ann = record.INFO.get("ANN")
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
                    term_set.add(f"@VARIANT_{hgvsp}_{gene}_{species}")

    return term_set
