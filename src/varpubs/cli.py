from simple_parsing import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from typing import List
from varpubs.pubmed_db import PubmedDB


@dataclass
class DeployDBArgs:
    """
    Command-line arguments for deploying the PubMed variant database.

    Attributes:
        db_path (Path): Path where the DuckDB database file will be saved.
        vcf_paths (List[Path]): List of VCF file paths to extract HGVS.p variants from.
        email (str): Email address required for querying the Entrez API.
    """

    db_path: Path
    vcf_paths: List[Path]
    email: str


def main():
    """
    Entry point for the CLI. Handles argument parsing and calls the deploy function.

    Supports the following subcommand:
        - deploy-db: Creates or updates a DuckDB database by extracting HGVS.p terms
        from VCFs and fetching related PubMed articles using Entrez.
    """
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: deploy-db
    deploy_parser = subparsers.add_parser("deploy-db", help="Deploy the database")
    deploy_parser.add_arguments(DeployDBArgs, dest="args")

    args = parser.parse_args()

    # Dispatch: if the subcommand is "deploy-db", initialize PubmedDB and deploy
    if args.command == "deploy-db":
        db = PubmedDB(
            path=args.args.db_path, vcf_paths=args.args.vcf_paths, email=args.args.email
        )
        db.deploy()
