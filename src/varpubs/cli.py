'''
-------------------------------------------------------------
Usage:
1) Deploy the database by running:

python -m varpubs.cli deploy-db \
    --db_path PATH/TO/OUTPUT_DB.duckdb \
    --vcf_paths PATH/TO/FILE1.vcf PATH/TO/FILE2.vcf ... \
    --email your_email@example.com

Description:
- Parses one or more VCF files to extract HGVS.p terms.
- Queries PubMed using Entrez API and stores abstracts in DuckDB.
- Requires your email for NCBI Entrez API access.

2) After running, to inspect the database:

Option A - Using DuckDB CLI:
duckdb PATH/TO/OUTPUT_DB.duckdb

Then inside the CLI, run:
> SELECT * FROM pubmedarticle;

Option B - Export table(s) to CSV from inside DuckDB:
> COPY pubmedarticle TO 'pubmed_articles.csv' (HEADER, DELIMITER ',');

-------------------------------------------------------------
'''

from simple_parsing import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from typing import List
from varpubs.pubmed_db import PubmedDB

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

if __name__ == "__main__":
    print(">>> CLI module invoked")
    main()
