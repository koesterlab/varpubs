from varpubs.settings import DeployPubmedDBSettings

# if storing in duckdb or sqlite, use the sqlmodel package


def deploy(settings: DeployPubmedDBSettings):
    # get pubmed data from their servers
    # store them in a suitable database format
    ...