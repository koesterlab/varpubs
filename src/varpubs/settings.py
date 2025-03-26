from pathlib import Path
from tap import Tap

class DeployPubmedDBSettings(Tap):
    output_path: Path