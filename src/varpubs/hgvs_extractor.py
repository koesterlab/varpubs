from cyvcf2 import VCF

def get_hgvsp_index(vcf: VCF) -> int:
    """
    Finds the index of the "HGVS.p" field in the VCF ANN header.

    This function parses the header of a VCF file to locate the index of the "HGVS.p"
    field in the ANN (annotation) metadata. This index is necessary to extract
    protein-level variant annotations from the ANN field.

    Args:
        vcf (VCF): A parsed cyvcf2.VCF object.

    Returns:
        int: The index of the "HGVS.p" field in the ANN description.

    Raises:
        RuntimeError: If the ANN field or HGVS.p is not found in the VCF header.
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
                except ValueError:
                    raise RuntimeError("HGVS.p not found in ANN header")
    raise RuntimeError("ANN field not found in VCF header")


def extract_hgvsp_from_vcf(vcf_path: str) -> list[str]:
    """
    Extracts HGVS.p terms from the ANN field of a VCF file.

    This function opens a VCF file, identifies the ANN annotation field, and 
    extracts all protein-level variant annotations (e.g., "p.Val600Glu") 
    along with their associated gene names. The returned values are formatted 
    as "GENE p.XXX".

    Args:
        vcf_path (str): Path to the VCF file to be parsed.

    Returns:
        list[str]: A list of strings in the format "GENE p.XXX".
    """
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
