from simple_parsing import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional
from varpubs.pubmed_db import PubmedDB
from varpubs.summarize import HFSettings, PubmedSummarizer


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
    - hf_token: Hugging Face API token for model access.
    - hf_model: The LLM model used for summarization (default: Mistral 7B).
    - output: Optional path to save the final variant summaries.
    """

    db_path: Path
    vcf_path: Path
    hf_token: str
    hf_model: str = "mistralai/Mistral-7B-Instruct-v0.1"
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

    # --- deploy-db command setup ---
    deploy_parser = subparsers.add_parser("deploy-db", help="Deploy the database")
    deploy_parser.add_arguments(DeployDBArgs, dest="args")

    # --- summarize-variants command setup ---
    summarize_parser = subparsers.add_parser(
        "summarize-variants", help="Summarize variants using LLM"
    )
    summarize_parser.add_arguments(SummarizeArgs, dest="args")

    args = parser.parse_args()

    if args.command == "deploy-db":
        # Build and deploy the PubMed database using VCF files
        db = PubmedDB(
            path=args.args.db_path, vcf_paths=args.args.vcf_paths, email=args.args.email
        )
        db.deploy()

    elif args.command == "summarize-variants":
        # Import summarization function and run it using provided model/token
        from varpubs.summarize_variants import summarize_variants

        summarizer = PubmedSummarizer(
            settings=HFSettings(
                token=args.args.hf_token,
                model=args.args.hf_model,
            )
        )

        summarize_variants(
            db_path=args.args.db_path,
            vcf_path=args.args.vcf_path,
            summarizer=summarizer,
            out_path=args.args.output,
        )
