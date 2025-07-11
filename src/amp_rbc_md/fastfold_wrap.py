from __future__ import annotations

"""Wrapper für FastFold / OpenFold Inferenz (Single-sequence).

Voraussetzungen
---------------
1. `openfold` + `deepspeed` müssen installiert sein (siehe environment.yml).
2. Die nötigen Datenbanken (MSA, Templates) und Modellgewichte müssen vorhanden sein.
   Wir erwarten, dass die Umgebungvariable ``OPENFOLD_DATA`` auf das Wurzelverzeichnis
   der Standard-AF2-Datenbanken zeigt (gleiches Layout wie im OpenFold-README).
3. Die (mitgelieferten) DeepMind-Gewichte reichen für Monomer-Vorhersagen.

Funktionsweise
--------------
Für eine übergebene Sequenz wird eine temporäre FASTA geschrieben und anschließend
``python -m openfold.run_pretrained_openfold`` als Sub-Prozess aufgerufen.
Das Skript legt die PDB im Ausgabe-Ordner ab. Wir parsen den Dateinamen und
geben den Pfad zurück.

Falls etwas fehlschlägt, wird **kein Dummy** erzeugt – stattdessen eine klare
Exception, damit Up-Stream-Code das Problem erkennt.
"""

from pathlib import Path
import os
import subprocess
import shutil
from typing import Optional

from .utils import LOGGER, ensure_dir

__all__ = ["predict"]


class FastFoldError(RuntimeError):
    """Fehler bei der FastFold-Inferenz."""


def _find_fastfold_script() -> list[str]:
    """Ermittle Aufruf für FastFold-Inferenz.

    Bevorzugt lokales FastFold inference.py, Fallback: OpenFold falls verfügbar.
    """

    # Prüfe, ob lokales FastFold inference.py existiert
    fastfold_repo = Path(__file__).parent.parent.parent / "fastfold_repo" / "inference.py"
    if fastfold_repo.exists():
        return ["python", str(fastfold_repo)]

    # Prüfe, ob fastfold als Python-Modul importierbar ist
    try:
        import importlib
        importlib.import_module("fastfold")
        # FastFold hat ein inference.py Skript
        return ["python", "-m", "fastfold.inference"]
    except ModuleNotFoundError:
        pass

    # Fallback: OpenFold (falls verfügbar)
    try:
        import importlib
        importlib.import_module("openfold")
        # Prüfe ob run_pretrained_openfold verfügbar ist
        try:
            importlib.import_module("openfold.run_pretrained_openfold")
            return ["python", "-m", "openfold.run_pretrained_openfold"]
        except ModuleNotFoundError:
            pass
    except ModuleNotFoundError:
        pass

    # Suche nach Skripten im PATH
    script = shutil.which("fastfold")
    if script:
        return [script]

    script = shutil.which("run_pretrained_openfold.py")
    if script:
        return ["python", script]

    raise FastFoldError(
        "Weder FastFold noch OpenFold verfügbar. Bitte 'pip install fastfold' oder 'pip install openfold' ausführen."
    )


def predict(
    sequence: str,
    out_dir: str | Path,
    *,
    db_root: str | Path | None = None,
    model_device: str | None = None,
) -> Path:
    """Führe FastFold/OpenFold aus und liefere Pfad zur erzeugten PDB.

    Parameters
    ----------
    sequence : str
        Aminosäure-Sequenz (1-Letter-FASTA).
    out_dir : Path | str
        Ausgabeverzeichnis (wird erstellt).
    db_root : Path | str | None, optional
        Wurzelpfad der Datenbanken. Standard: ``$OPENFOLD_DATA``.
    model_device : str | None, optional
        GPU-Device (z.B. "cuda:0"). Lässt sich ansonsten von OpenFold wählen.
    """

    out_dir = ensure_dir(out_dir)
    fasta_path = Path(out_dir) / "query.fasta"
    fasta_path.write_text(f">query\n{sequence}\n")

    # Wenn bereits PDB existiert, nichts tun
    pdb_path = Path(out_dir) / "query_unrelaxed_model_1.pdb"
    if pdb_path.exists():
        LOGGER.info("FastFold: PDB existiert bereits unter %s", pdb_path)
        return pdb_path

    # Datenbank-Root ermitteln
    db_root = Path(db_root or os.environ.get("OPENFOLD_DATA", ""))
    if not db_root.exists():
        raise FastFoldError(
            "OPENFOLD_DATA nicht gesetzt oder Verzeichnis existiert nicht: %s" % db_root
        )

    mmcif_dir = db_root / "pdb_mmcif" / "mmcif_files"
    if not mmcif_dir.exists():
        raise FastFoldError("pdb_mmcif/mmcif_files nicht gefunden unter %s" % db_root)

    cmd = _find_fastfold_script() + [
        str(fasta_path),
        str(mmcif_dir),
        "--output_dir",
        str(out_dir),
        "--gpus",
        "1",  # Single GPU für jetzt
        "--enable_workflow",
        "--inplace",
    ]

    # Häufig benötigte Datenbanken
    db_flags = {
        "--uniref90_database_path": db_root / "uniref90" / "uniref90.fasta",
        "--mgnify_database_path": db_root / "mgnify" / "mgy_clusters_2018_12.fa",
        "--pdb70_database_path": db_root / "pdb70" / "pdb70",
        "--uniclust30_database_path": db_root / "uniclust30" / "uniclust30_2018_08" / "uniclust30_2018_08",
        "--bfd_database_path": db_root / "bfd" / "bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
    }
    for flag, path in db_flags.items():
        if path.exists():
            cmd += [flag, str(path)]
        else:
            LOGGER.warning("FastFold: Datenbank fehlt (%s), Flag wird ausgelassen.", path)

    if model_device:
        cmd += ["--model_device", model_device]

    LOGGER.info("Starte FastFold/OpenFold: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        raise FastFoldError(f"FastFold/OpenFold schlug fehl: {exc}") from exc

    # Finde erzeugte PDB (OpenFold legt *_unrelaxed_model_1.pdb an)
    pdbcandidates = list(Path(out_dir).glob("*_model_*.pdb"))
    if not pdbcandidates:
        raise FastFoldError("OpenFold hat keine PDB erzeugt")

    pdb_path = pdbcandidates[0]
    LOGGER.info("FastFold: Struktur erzeugt (%s)", pdb_path)
    return pdb_path 