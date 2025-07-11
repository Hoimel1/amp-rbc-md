from __future__ import annotations

from pathlib import Path
from typing import Literal, Optional

from Bio.Seq import Seq  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore
from Bio import SeqIO  # type: ignore

from .utils import LOGGER, ensure_dir, set_seed

# Optional ColabFold import
try:
    from colabfold import batch  # type: ignore

    COLABFOLD_AVAILABLE: bool = True
except ImportError:  # pragma: no cover
    COLABFOLD_AVAILABLE = False

# FastFold/OpenFold Wrapper
try:
    from .fastfold_wrap import predict as fastfold_predict  # noqa: WPS433

    FASTFOLD_AVAILABLE: bool = True
except Exception:  # pragma: no cover
    FASTFOLD_AVAILABLE = False

Engine = Literal["fastfold", "colabfold"]

__all__ = ["fasta_to_pdb", "Engine"]


class StructurePredictionError(RuntimeError):
    """Generischer Fehler der Strukturvorhersage."""


def _run_colabfold(sequence: str, out_dir: Path) -> Path:
    if not COLABFOLD_AVAILABLE:
        raise StructurePredictionError("ColabFold nicht installiert.")

    from colabfold.utils import setup_logging  # type: ignore

    result_dir = out_dir / "colabfold_results"
    result_dir.mkdir(exist_ok=True)

    setup_logging(result_dir / "colabfold.log")

    LOGGER.info("Starte ColabFold für Strukturvorhersage …")
    queries = [("peptide", sequence, None)]

    batch.run(
        queries=queries,
        result_dir=str(result_dir),
        num_models=1,
        is_complex=False,
        num_recycles=3,
        model_type="auto",
        msa_mode="mmseqs2_uniref_env",
        use_templates=False,
    )

    pdb_files = list(result_dir.glob("*.pdb"))
    if not pdb_files:
        raise StructurePredictionError("ColabFold hat keine PDB erzeugt")

    pdb_path = out_dir / "model.pdb"
    pdb_files[0].rename(pdb_path)
    return pdb_path


def _run_fastfold(sequence: str, out_dir: Path) -> Path:
    if not FASTFOLD_AVAILABLE:
        raise StructurePredictionError("FastFold/OpenFold nicht installiert.")

    return fastfold_predict(sequence, out_dir)


def fasta_to_pdb(
    sequence: str,
    out_dir: str | Path,
    *,
    engine: Engine = "fastfold",
    db_root: Optional[str | Path] = None,
) -> Path:
    """Generiere eine PDB-Datei für *sequence* in *out_dir*.

    Parameters
    ----------
    sequence : str
        Aminosäure-Sequenz.
    out_dir : str | Path
        Verzeichnis für Ausgabedateien.
    engine : {"fastfold", "colabfold"}
        Backend für Strukturvorhersage.
    db_root : Path | str | None
        Wird an FastFold weitergereicht (s. ``OPENFOLD_DATA``).
    """

    set_seed()
    out_dir = ensure_dir(out_dir)

    # FASTA parallel ablegen (immer)
    fasta_path = Path(out_dir) / "sequence.fasta"
    SeqIO.write(
        SeqRecord(Seq(sequence), id="peptide", description=""), fasta_path, "fasta"
    )

    if engine == "fastfold":
        try:
            return _run_fastfold(sequence, out_dir)
        except Exception as err:
            LOGGER.error("FastFold fehlgeschlagen: %s", err)
            raise StructurePredictionError(err) from err
    elif engine == "colabfold":
        try:
            return _run_colabfold(sequence, out_dir)
        except Exception as err:
            LOGGER.error("ColabFold fehlgeschlagen: %s", err)
            raise StructurePredictionError(err) from err
    else:  # pragma: no cover
        raise ValueError(f"Unbekanntes engine='{engine}'")
