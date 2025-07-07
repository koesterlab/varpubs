from simple_parsing import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from typing import List
from varpubs.pubmed_db import PubmedDB

"""
See README.md for CLI usage instructions.
"""

@dataclass
class DeployDBArgs:
    db_path: Path
    vcf_paths: List[Path]
    email: str

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    deploy_parser = subparsers.add_parser("deploy-db", help="Deploy the database")
    deploy_parser.add_arguments(DeployDBArgs, dest="args")

    args = parser.parse_args()

    if args.command == "deploy-db":
        db = PubmedDB(
            path=args.args.db_path,
            vcf_paths=args.args.vcf_paths,
            email=args.args.email
        )
        db.deploy()
