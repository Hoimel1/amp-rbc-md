from __future__ import annotations

import random
from pathlib import Path
from textwrap import dedent
from typing import Tuple

from Bio.Seq import Seq  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore
from Bio import SeqIO  # type: ignore

try:
    from colabfold import batch
    COLABFOLD_AVAILABLE = True
except ImportError:
    COLABFOLD_AVAILABLE = False

from .utils import LOGGER, ensure_dir, set_seed


def fasta_to_pdb(sequence: str, out_dir: str | Path) -> Path:
    """Konvertiere eine Peptidsequenz in eine PDB-Datei mit ColabFold (AlphaFold).

    Verwendet ColabFold für State-of-the-Art Strukturvorhersage basierend auf AlphaFold2.
    """
    set_seed()
    out_dir = ensure_dir(out_dir)
    pdb_path = Path(out_dir) / "model.pdb"

    if not COLABFOLD_AVAILABLE:
        raise RuntimeError(
            "ColabFold ist nicht verfügbar. Bitte installieren Sie es mit:\n"
            "pip install colabfold[alphafold]"
        )

    try:
        LOGGER.info("Verwende ColabFold (AlphaFold2) für State-of-the-Art Strukturvorhersage")
        
        # ColabFold für Strukturvorhersage verwenden
        # ColabFold ist eine optimierte AlphaFold2-Implementierung
        from colabfold.utils import setup_logging
        
        # Logging konfigurieren
        log_file = Path(out_dir) / "colabfold.log"
        setup_logging(log_file)
        
        # Sequenz für ColabFold vorbereiten
        # ColabFold erwartet ein spezifisches Format
        job_name = "peptide_prediction"
        sequences = [sequence]
        
        # ColabFold-Batch-Ausführung
        result_dir = Path(out_dir) / "colabfold_results"
        result_dir.mkdir(exist_ok=True)
        
        # Führe ColabFold aus
        batch.run(
            sequences=sequences,
            job_name=job_name,
            result_dir=str(result_dir),
            use_amber=False,  # Schneller ohne AMBER-Refinement
            use_templates=False,  # Keine Templates für Peptide
            custom_msa=None,  # Keine benutzerdefinierten MSA
            pair_mode="unpaired+paired",  # Standard für Monomere
            host_url="https://api.colabfold.com",  # Standard-Server
        )
        
        # Finde die beste PDB-Datei (normalerweise die erste)
        pdb_files = list(result_dir.glob("*.pdb"))
        if not pdb_files:
            raise RuntimeError("ColabFold hat keine PDB-Datei erzeugt")
        
        # Kopiere die beste PDB-Datei
        best_pdb = pdb_files[0]
        import shutil
        shutil.copy2(best_pdb, pdb_path)
        
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
