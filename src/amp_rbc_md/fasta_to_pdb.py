from __future__ import annotations

import random
from pathlib import Path
from textwrap import dedent
from typing import Tuple

from Bio.Seq import Seq  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore
from Bio import SeqIO  # type: ignore

try:
    import colabfold
    COLABFOLD_AVAILABLE = True
except ImportError:
    COLABFOLD_AVAILABLE = False

from .utils import LOGGER, ensure_dir, set_seed


def fasta_to_pdb(sequence: str, out_dir: str | Path) -> Path:
    """Konvertiere eine Peptidsequenz in eine PDB-Datei mit ColabFold.

    Verwendet ColabFold für State-of-the-Art Strukturvorhersage.
    """
    set_seed()
    out_dir = ensure_dir(out_dir)
    pdb_path = Path(out_dir) / "model.pdb"

    if not COLABFOLD_AVAILABLE:
        raise RuntimeError(
            "ColabFold ist nicht verfügbar. Bitte installieren Sie es mit:\n"
            "pip install colabfold"
        )

    try:
        LOGGER.info("Verwende ColabFold für State-of-the-Art Strukturvorhersage")
        
        # Erstelle temporäre FASTA-Datei
        temp_fasta = Path(out_dir) / "temp.fasta"
        temp_fasta.write_text(f">peptide\n{sequence}\n")
        
        # ColabFold Vorhersage
        from colabfold.batch import get_queries, run
        from colabfold.download import download_alphafold_params
        
        # Lade AlphaFold Parameter für Monomere
        download_alphafold_params("alphafold2_ptm")
        
        # Führe Vorhersage aus (korrigierte API)
        queries = get_queries(str(temp_fasta))
        results = run(queries, result_dir=str(out_dir), use_templates=False, custom_template_path=None)
        
        # Finde die beste PDB-Datei
        pdb_files = list(out_dir.glob("*.pdb"))
        if pdb_files:
            best_pdb = pdb_files[0]  # Nimm die erste gefundene PDB
            best_pdb.rename(pdb_path)
        else:
            raise RuntimeError("ColabFold konnte keine PDB-Datei erzeugen")
        
        # Lösche temporäre Dateien
        temp_fasta.unlink(missing_ok=True)
        
        LOGGER.info("ColabFold Struktur unter %s erzeugt", pdb_path)
        
    except Exception as e:
        LOGGER.error("ColabFold fehlgeschlagen: %s", e)
        raise RuntimeError(f"ColabFold konnte keine Struktur für Sequenz '{sequence}' vorhersagen: {e}")

    # Speichere FASTA parallel, nützlich für Referenz.
    fasta_path = Path(out_dir) / "sequence.fasta"
    SeqIO.write(
        SeqRecord(Seq(sequence), id="peptide", description=""), fasta_path, "fasta"
    )

    return pdb_path


__all__ = ["fasta_to_pdb"]
