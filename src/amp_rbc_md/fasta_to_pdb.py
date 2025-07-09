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

    Verwendet ESMFold v1 für State-of-the-Art Strukturvorhersage.
    """
    set_seed()
    out_dir = ensure_dir(out_dir)
    pdb_path = Path(out_dir) / "model.pdb"

    if not ESM_AVAILABLE:
        raise RuntimeError("ESMFold ist nicht verfügbar. Bitte installieren Sie es mit: pip install fair-esm")

    try:
        LOGGER.info("Verwende ESMFold v1 für State-of-the-Art Strukturvorhersage")
        
        # Lade ESMFold v1 Modell
        model = esm.pretrained.esmfold_v1()
        
        # Vorhersage der Struktur
        output = model.infer_pdb(sequence)
        
        # Speichere PDB-Datei
        pdb_path.write_text(output)
        
        LOGGER.info("ESMFold v1 Struktur unter %s erzeugt", pdb_path)
        
    except ImportError as e:
        if "openfold" in str(e):
            raise RuntimeError(
                "ESMFold v1 benötigt openfold. Bitte installieren Sie es mit:\n"
                "pip install openfold\n"
                "oder\n"
                "conda install -c conda-forge openfold"
            )
        else:
            raise RuntimeError(f"ESMFold konnte nicht geladen werden: {e}")
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
