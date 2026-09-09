import csv
import logging
from pathlib import Path
from typing import Any, List, Optional, Tuple

from cyvcf2 import VCF
from hgvs.exceptions import HGVSParseError
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
            if len(fields) > max(hgvsp_index, gene_index):
                hgvsp = fields[hgvsp_index]
                try:
                    hgvsp_single = (
                        hgvs_parser.parse(hgvsp.replace("%3D", "="))
                        .format(conf={"p_3_letter": False})
                        .split(":")[1]
                    )
                except HGVSParseError as e:
                    logger.warning(f"Unable to parse hgvsp: '{hgvsp}'\n{e}")
                    hgvsp_single = hgvsp.replace("%3D", "=")
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


def table_row_to_bioconcept(gene: str, variant: str, species: str) -> str:
    # The variant column is expected to already hold a one-letter HGVS protein
    # change (e.g. p.N331I). Strip surrounding whitespace (including non-breaking
    # spaces from spreadsheet exports) so the bioconcept matches the database.
    return hgvsp_gene_to_bioconcept(variant.strip(), gene.strip(), species)


def table_delimiter(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def extract_bioconcepts_from_table(
    table_path: str,
    gene_column: str,
    variant_column: str,
    species: str,
    delimiter: Optional[str] = None,
) -> set[str]:
    delimiter = delimiter or table_delimiter(Path(table_path))
    with open(table_path, newline="") as infile:
        reader = csv.DictReader(infile, delimiter=delimiter)
        columns = reader.fieldnames or []
        missing = [c for c in (gene_column, variant_column) if c not in columns]
        if missing:
            raise ValueError(
                f"Column(s) {missing} not found in {table_path}. Available: {columns}"
            )
        return {
            table_row_to_bioconcept(row[gene_column], row[variant_column], species)
            for row in reader
        }
