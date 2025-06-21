from cyvcf2 import VCF

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
                        hgvsp_list.append(f"{gene} {hgvsp}")  # değişiklik burada

    return hgvsp_list
