from __future__ import annotations

import random
from pathlib import Path
from textwrap import dedent
from typing import Tuple

from Bio.Seq import Seq  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore
from Bio import SeqIO  # type: ignore

try:
    import esm
    ESM_AVAILABLE = True
except ImportError:
    ESM_AVAILABLE = False

from .utils import LOGGER, ensure_dir, set_seed


def fasta_to_pdb(sequence: str, out_dir: str | Path) -> Path:
    """Konvertiere eine Peptidsequenz in eine PDB-Datei mit ESMFold.

    Verwendet ESMFold für State-of-the-Art Strukturvorhersage.
    """
    set_seed()
    out_dir = ensure_dir(out_dir)
    pdb_path = Path(out_dir) / "model.pdb"

    if not ESM_AVAILABLE:
        raise RuntimeError(
            "ESMFold ist nicht verfügbar. Bitte installieren Sie es mit:\n"
            "pip install fair-esm"
        )

    try:
        LOGGER.info("Verwende ESMFold für State-of-the-Art Strukturvorhersage")
        
        # Lade aktuelles ESMFold-Modell (liefert vollständige Heavy-Atom Koordinaten)
        model = esm.pretrained.esmfold_v1()
        model = model.eval()

        import torch
        with torch.no_grad():
            output = model.infer_pdb(sequence)

        # Speichere PDB-Datei mit vollständigen Atomen
        pdb_path.write_text(output)
        
        LOGGER.info("ESMFold Struktur unter %s erzeugt", pdb_path)
        
    except Exception as e:
        LOGGER.error("ESMFold fehlgeschlagen: %s", e)
        raise RuntimeError(f"ESMFold konnte keine Struktur für Sequenz '{sequence}' vorhersagen: {e}")

    # Speichere FASTA parallel, nützlich für Referenz.
    fasta_path = Path(out_dir) / "sequence.fasta"
    SeqIO.write(
        SeqRecord(Seq(sequence), id="peptide", description=""), fasta_path, "fasta"
    )

    return pdb_path


__all__ = ["fasta_to_pdb"]
