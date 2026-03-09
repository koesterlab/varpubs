from typing import List

from cyvcf2 import VCF


def extend_vep_header(vcf: VCF, new_fields: List[str], csq_id="CSQ"):
    description = vcf.get_header_type(csq_id)["Description"]
    parts = description.split("'")
    to_add = " | ".join(new_fields)
    new_description = f"{parts[0][1:]}'{parts[1]} | {to_add}'"
    vcf.remove_header(csq_id)
    vcf.add_info_to_header(
        {
            "ID": csq_id,
            "Number": ".",
            "Type": "String",
            "Description": new_description,
        }
    )
