import logging
from simple_parsing import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional
from varpubs.pubmed_db import PubmedDB
from varpubs.summarize import Settings, PubmedSummarizer


@dataclass
class DeployDBArgs:
    """
    Command-line arguments for deploying the PubMed variant database.

    - db_path: Path to the DuckDB database file to be created or updated.
    - vcf_paths: List of VCF files containing variant information.
    - email: Email address used for Entrez API access.
    """

    db_path: Path
    vcf_paths: List[Path]
    email: str


@dataclass
class SummarizeArgs:
    """
    Command-line arguments for summarizing PubMed articles related to variants.

    - db_path: Path to the existing DuckDB database file.
    - vcf_path: A single annotated VCF file with variant terms.
    - api_key: Hugging Face API token for model access.
    - llm_url: Base URL for LLM API (Must follow the openai API format)
    - model: The LLM model used for summarization (default: teuken-7b-instruct-research).
    - output: Optional path to save the final variant summary file (CSV).
    """

    db_path: Path
    vcf_path: Path
    api_key: str
    llm_url: str
    model: str = "teuken-7b-instruct-research"
    output: Optional[Path] = None


def main():
    """
    Entry point for the varpubs CLI tool.

    Supports two commands:
    - 'deploy-db': Parses VCFs and populates the DuckDB database with PubMed entries.
    - 'summarize-variants': Summarizes articles related to variants using an LLM model.
    """
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity level (use -v, -vv, -vvv)",
    )

    deploy_parser = subparsers.add_parser("deploy-db", help="Deploy the database")
    deploy_parser.add_arguments(DeployDBArgs, dest="args")

    summarize_parser = subparsers.add_parser(
        "summarize-variants", help="Summarize variants using LLM"
    )
    summarize_parser.add_arguments(SummarizeArgs, dest="args")

    args = parser.parse_args()

    level = logging.WARNING
    if args.verbose == 1:
        level = logging.INFO
    elif args.verbose >= 2:
        level = logging.DEBUG
    logging.basicConfig(level=level)

    if args.command == "deploy-db":
        db = PubmedDB(
            path=args.args.db_path, vcf_paths=args.args.vcf_paths, email=args.args.email
        )
        db.deploy()

    elif args.command == "summarize-variants":
        from varpubs.summarize_variants import summarize_variants

        summarizer = PubmedSummarizer(
            settings=Settings(
                api_key=args.args.api_key,
                model=args.args.model,
                base_url=args.args.llm_url,
            )
        )

        summarize_variants(
            db_path=args.args.db_path,
            vcf_path=args.args.vcf_path,
            summarizer=summarizer,
            out_path=args.args.output,
        )
