import re
from cyvcf2 import VCF

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


def extract_hgvsp_from_vcf(vcf_path: str) -> list[str]:
    vcf = VCF(vcf_path)
    hgvsp_index = get_annotation_field_index(vcf, "HGVS.p")
    gene_index = get_annotation_field_index(vcf, "Gene_Name")
    term_set = set()

    for record in vcf:
        ann = record.INFO.get("ANN")
        if ann:
            for ann_entry in ann.split(","):
                fields = ann_entry.split("|")
                if len(fields) > max(hgvsp_index, gene_index):
                    hgvsp = fields[hgvsp_index]
                    gene = fields[gene_index]

                    if not hgvsp.startswith("p."):
                        continue

                    # Match full HGVS.p format: e.g., p.Gly12Cys
                    match = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})?", hgvsp)
                    if not match:
                        continue

                    ref_aa_3 = match.group(1)
                    pos = match.group(2)
                    alt_aa_3 = match.group(3)

                    # Skip synonymous mutations (e.g., p.Tyr516Tyr)
                    if alt_aa_3 and ref_aa_3 == alt_aa_3:
                        continue

                    ref_aa_1 = AA3_TO_1.get(ref_aa_3)
                    alt_aa_1 = AA3_TO_1.get(alt_aa_3) if alt_aa_3 else None

                    # Add full term: e.g., 'KRAS p.Gly12Cys'
                    term_set.add(f"{gene} {hgvsp}")

                    # Add short form(s): e.g., 'KRAS G12' and 'KRAS G12C'
                    if ref_aa_1:
                        term_set.add(f"{gene} {ref_aa_1}{pos}")
                        if alt_aa_1:
                            term_set.add(f"{gene} {ref_aa_1}{pos}{alt_aa_1}")

    return sorted(term_set)
