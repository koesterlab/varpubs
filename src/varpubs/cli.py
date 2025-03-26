from simple_parsing import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from varpubs.pubmed_db import PubmedDB


@dataclass
class Settings:
    db_path: Path


def get_argument_parser() -> ArgumentParser:
    parser: ArgumentParser = ArgumentParser()
    parser.add_arguments(Settings, dest="options")
    subparsers = parser.add_subparsers(help="subcommand help")
    subparsers.add_parser("deploy-db", help="Deploy the database")
    return parser


def main():
    args = get_argument_parser().parse_args()
    db = PubmedDB(path=args.db_path)
    if args.subparser_name == "deploy-db":
        db.deploy()
