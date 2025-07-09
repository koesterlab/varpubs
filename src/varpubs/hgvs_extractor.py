import re
from cyvcf2 import VCF


def get_hgvsp_index(vcf: VCF) -> int:
    """
    Finds the index of the 'HGVS.p' field in the ANN (annotation) column of the VCF header.

    Args:
        vcf (VCF): A parsed VCF object from cyvcf2.

    Returns:
        int: The index position of the 'HGVS.p' field in ANN entries.

    Raises:
        RuntimeError: If 'ANN' or 'HGVS.p' is not found in the VCF header.
    """
    for rec in vcf.header_iter():
        info = rec.info()
        if info.get("ID") == "ANN":
            desc = info.get("Description", "")
            if "HGVS.p" in desc:
                fields_str = desc.split("': '")[-1].rstrip("'")
                fields = [f.strip() for f in fields_str.split("|")]
                try:
                    return fields.index("HGVS.p")
                except ValueError as e:
                    raise RuntimeError("HGVS.p not found in ANN header") from e
    raise RuntimeError("ANN field not found in VCF header")


def extract_hgvsp_from_vcf(vcf_path: str) -> list[str]:
    """
    Extracts all unique HGVS.p terms from the given VCF file and converts them into two formats:
    - Full form (e.g., 'KRAS p.Gly12Cys')
    - Short form (e.g., 'KRAS G12')

    Filters out synonymous (silent) mutations.

    Args:
        vcf_path (str): Path to the annotated VCF file.

    Returns:
        list[str]: A sorted list of unique variant terms to be searched in PubMed.
    """
    vcf = VCF(vcf_path)
    hgvsp_index = get_hgvsp_index(vcf)
    gene_index = 3  # Fixed index for gene symbol in ANN field
    term_set = set()

    for record in vcf:
        ann = record.INFO.get("ANN")
        if ann:
            for ann_entry in ann.split(","):
                fields = ann_entry.split("|")
                if len(fields) > max(hgvsp_index, gene_index):
                    hgvsp = fields[hgvsp_index]
                    gene = fields[gene_index]

                    # Skip entries not starting with 'p.' (protein-level annotation)
                    if not hgvsp.startswith("p."):
                        continue

                    # Skip synonymous mutations (e.g., p.Tyr516Tyr)
                    aa_change = re.findall(r"[A-Za-z]+", hgvsp)
                    if len(aa_change) == 2 and aa_change[0] == aa_change[1]:
                        continue

                    # Add full term: e.g., 'KRAS p.Gly12Cys'
                    term_set.add(f"{gene} {hgvsp}")

                    # Add short form: e.g., 'KRAS G12'
                    match = re.match(r"p\.\D+(\d+)", hgvsp)
                    if match:
                        term_set.add(f"{gene} G{match.group(1)}")

    return sorted(term_set)
