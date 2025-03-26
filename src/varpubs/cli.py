from pathlib import Path

from varpubs.pubmed_db import deploy
from varpubs.settings import DeployPubmedDBSettings


def main():
    settings = DeployPubmedDBSettings().parse_args()
    deploy(settings)
