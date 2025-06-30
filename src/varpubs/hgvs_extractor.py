# -----------------------------------------------------------------------------
# This module provides functions to extract HGVS.p terms from a VCF file.
# Specifically, it parses the VCF ANN (annotation) field to identify protein-level 
# variant annotations (e.g., "p.Val600Glu") and returns them in the form "GENE p.XXX".
# These extracted terms can be used for further analysis or querying external databases 
# such as PubMed.
# -----------------------------------------------------------------------------

from cyvcf2 import VCF

# Finds the index of the "HGVS.p" field in the VCF ANN header.
# Returns the index position to locate HGVS.p terms in annotation entries.
def get_hgvsp_index(vcf: VCF) -> int:
    for rec in vcf.header_iter():
        info = rec.info()
        if info.get("ID") == "ANN":
            desc = info.get("Description", "")
            if "HGVS.p" in desc:
                fields_str = desc.split("': '")[-1].rstrip("'")
                fields = [f.strip() for f in fields_str.split("|")]
                try:
                    return fields.index("HGVS.p")
                except ValueError:
                    raise RuntimeError("HGVS.p not found in ANN header")
    raise RuntimeError("ANN field not found in VCF header")

# Extracts protein-level variant annotations (HGVS.p) from a VCF file.
# Returns a list of strings formatted as "GENE p.XXX".
def extract_hgvsp_from_vcf(vcf_path: str) -> list[str]:
    vcf = VCF(vcf_path)
    hgvsp_index = get_hgvsp_index(vcf)
    gene_index = 3  
    hgvsp_list = []
    for record in vcf:
        ann = record.INFO.get("ANN")
        if ann:
            for ann_entry in ann.split(","):
                fields = ann_entry.split("|")
                if len(fields) > max(hgvsp_index, gene_index):
                    hgvsp = fields[hgvsp_index]
                    gene = fields[gene_index]
                    if hgvsp.startswith("p."):
                        hgvsp_list.append(f"{gene} {hgvsp}")  

    return hgvsp_list
