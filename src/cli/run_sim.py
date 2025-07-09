from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

import click  # type: ignore
from rich.table import Table  # type: ignore
from rich.console import Console  # type: ignore

from amp_rbc_md import analyse, build_membrane, fasta_to_pdb, gen_mdp, gmx_runner, judge
from amp_rbc_md.martinize_wrap import martinize
from amp_rbc_md.utils import ensure_dir, hash_sequence

CONSOLE = Console()


@click.command("run_sim")
@click.option("--seq", "sequence", type=str, help="Peptidsequenz (FASTA, 1-Letter)")
@click.option(
    "-f",
    "--fasta",
    "fasta_file",
    type=click.Path(exists=True, dir_okay=False),
    help="FASTA-Datei",
)
@click.option("--profile", type=str, default="default", help="Lipid-Profil")
@click.option("--n-replica", "n_replica", type=int, default=1, help="Anzahl Replika")
@click.option(
    "--time", "prod_time", type=str, default="auto", help="Produktionszeit, z.B. 500ns"
)
@click.option("--gpu", "gpu_id", type=str, help="GPU-ID, z.B. '0'")
@click.option(
    "--dry-run", is_flag=True, help="Nur Schritte auflisten, nichts ausführen"
)
def run_sim(
    sequence: Optional[str],
    fasta_file: Optional[str],
    profile: str,
    n_replica: int,
    prod_time: str,
    gpu_id: Optional[str],
    dry_run: bool,
) -> None:  # noqa: D401
    """Führe komplette Pipeline für 1 Sequenz oder FASTA-Batch aus."""

    if not sequence and not fasta_file:
        CONSOLE.print("[red]Fehler: --seq oder --fasta erforderlich.")
        sys.exit(1)

    seqs: list[str]
    if sequence:
        seqs = [sequence]
    else:
        assert fasta_file is not None  # für mypy
        seqs = [
            line.strip()
            for line in Path(fasta_file).read_text().splitlines()
            if line and not line.startswith(">")
        ]
        if len(seqs) > 10:
            CONSOLE.print("[red]Batch-Größe >10 nicht erlaubt.")
            sys.exit(1)

    report_rows = []

    for seq in seqs:
        seq_hash = hash_sequence(seq)
        base_dir = ensure_dir(Path("results") / seq_hash)
        # Strukturvorbereitung
        pdb_path = fasta_to_pdb.fasta_to_pdb(seq, base_dir / "01_pdb")
        cg_gro, _top = martinize(pdb_path, base_dir / "02_cg")
        system_gro = build_membrane.build(profile, cg_gro, base_dir / "03_system")

        # MDP-Dateien
        length_aa = len(seq)
        _mdp_prod = gen_mdp.generate(
            "prod", base_dir / "mdp", length_aa=length_aa, override_nsteps=None
        )
        # Hier würde noch grompp ausgeführt; überspringen wir.
        tpr_stub = base_dir / "04_tpr" / "system.tpr"
        ensure_dir(tpr_stub.parent)
        tpr_stub.write_text("stub tpr")

        for rep in range(n_replica):
            rep_dir = ensure_dir(base_dir / f"rep{rep}")
            runner = gmx_runner.GromacsRunner(
                rep_dir, gpu_id=gpu_id, replica_id=rep, dry_run=dry_run
            )
            trr = runner.run(tpr_stub)
            metrics = analyse.analyse_trajectory(trr, rep_dir)
            verdict = judge.judge_metrics(**metrics)
            report_rows.append(
                (seq_hash, rep, *metrics.values(), verdict.label, verdict.confidence)
            )

    # Zusammenfassung
    table = Table(title="Simulation Summary")
    table.add_column("Hash")
    table.add_column("Rep")
    table.add_column("ΔG")
    table.add_column("CI95")
    table.add_column("Thinning")
    table.add_column("Label")
    table.add_column("Conf")

    for row in report_rows:
        table.add_row(*map(str, row))
    CONSOLE.print(table)


if __name__ == "__main__":
    run_sim()  # type: ignore[misc]
