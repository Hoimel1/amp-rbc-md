from importlib import metadata

__all__ = [
    "fasta_to_pdb",
    "martinize_wrap",
    "build_membrane",
    "gen_mdp",
    "gmx_runner",
    "analyse",
    "judge",
]

try:
    __version__: str = metadata.version("amp-rbc-md")
except metadata.PackageNotFoundError:
    __version__ = "0.0.0"
