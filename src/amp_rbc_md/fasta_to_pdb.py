from __future__ import annotations

import random
from pathlib import Path
from textwrap import dedent
from typing import Tuple

from Bio.Seq import Seq  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore
from Bio import SeqIO  # type: ignore

from .utils import LOGGER, ensure_dir, set_seed


def fasta_to_pdb(sequence: str, out_dir: str | Path) -> Path:
    """Konvertiere eine Peptidsequenz in eine Dummy-PDB-Datei.

    In Produktionsumgebungen würde hier ColabFold/ESMFold aufgerufen.
    Für die CI erzeugen wir deterministische Zufallskoordinaten, damit die
    Downstream-Pipeline agieren kann, ohne die externen Modelle zu benötigen.
    """
    set_seed()
    out_dir = ensure_dir(out_dir)
    pdb_path = Path(out_dir) / "model_dummy.pdb"

    # Schreibe minimalen PDB-Header mit zufälligen CA-Positionen
    coords: Tuple[float, float, float]
    lines: list[str] = []
    for idx, aa in enumerate(sequence, start=1):
        x, y, z = (random.uniform(-1, 1) * 0.5 * idx for _ in range(3))
        lines.append(
            (
                "ATOM  {idx:5d}  CA  {aa:>3s} A{res:4d}    "
                "{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C".format(
                    idx=idx,
                    aa=aa,
                    res=idx,
                    x=x,
                    y=y,
                    z=z,
                )
            )
        )
    lines.append("END")
    pdb_path.write_text("\n".join(lines))

    # Speichere FASTA parallel, nützlich für Referenz.
    fasta_path = Path(out_dir) / "sequence.fasta"
    SeqIO.write(
        SeqRecord(Seq(sequence), id="peptide", description=""), fasta_path, "fasta"
    )

    LOGGER.info("Dummy-PDB unter %s erzeugt", pdb_path)
    return pdb_path


__all__ = ["fasta_to_pdb"]
