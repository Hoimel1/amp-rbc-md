from __future__ import annotations

import shutil
from pathlib import Path
from typing import Tuple

from .utils import LOGGER

try:
    import MDAnalysis as mda  # type: ignore
    from moviepy.editor import ImageSequenceClip  # type: ignore
    import matplotlib.pyplot as plt  # type: ignore
    import matplotlib  # type: ignore

    matplotlib.use("Agg")  # Non-interactive backend
except ModuleNotFoundError:  # pragma: no cover
    mda = None  # type: ignore


DEFAULT_FPS = 24
MAX_FRAMES = 240  # Beschränke Video-Länge


def traj_to_mp4(
    traj: str | Path,
    gro: str | Path | None,
    out_mp4: str | Path,
    *,
    fps: int = DEFAULT_FPS,
) -> Path:
    """Rendere ein MP4 der Trajektorie mit echten Daten.

    Diese Funktion erstellt echte Visualisierungen der Trajektorie
    und wirft Fehler bei fehlenden Abhängigkeiten.
    """
    traj = Path(traj)
    out_mp4 = Path(out_mp4)

    # Prüfe Abhängigkeiten
    if mda is None:
        raise RuntimeError(
            "MDAnalysis ist für Rendering erforderlich. "
            "Installieren Sie es mit: pip install MDAnalysis"
        )

    if not traj.exists():
        raise FileNotFoundError(f"Trajektorie {traj} existiert nicht")

    # Finde Gro-Datei
    if gro is None:
        gro_guess = traj.with_suffix(".gro")
        if not gro_guess.exists():
            raise FileNotFoundError(
                f"Gro-Datei {gro_guess} nicht gefunden. "
                "Gro-Datei ist für Rendering erforderlich."
            )
        gro = gro_guess

    gro_path = Path(gro)
    if not gro_path.exists():
        raise FileNotFoundError(f"Gro-Datei {gro_path} existiert nicht")

    LOGGER.info("Rendere Trajektorie: %s", traj)

    # Lade Trajektorie
    try:
        u = mda.Universe(str(gro), str(traj))
    except Exception as e:
        raise RuntimeError(f"Fehler beim Laden der Trajektorie: {e}")

    # Wähle Protein-Atome für Visualisierung
    try:
        protein = u.select_atoms("protein or name CA")
        if len(protein) == 0:
            # Fallback: alle Atome
            protein = u.select_atoms("all")
    except Exception as e:
        raise RuntimeError(f"Fehler bei Atom-Auswahl: {e}")

    # Erstelle temporäres Verzeichnis für Frames
    img_dir = out_mp4.parent / "frames_tmp"
    img_dir.mkdir(parents=True, exist_ok=True)

    try:
        images = []
        frame_count = 0

        for i, ts in enumerate(u.trajectory):
            if i >= MAX_FRAMES:
                LOGGER.info("Maximale Frame-Anzahl erreicht (%d)", MAX_FRAMES)
                break

            try:
                # Erstelle 3D-Plot
                fig = plt.figure(figsize=(8, 6))
                ax = fig.add_subplot(111, projection="3d")

                # Konvertiere Koordinaten von Å zu nm
                coords = protein.positions / 10

                # Zeichne Protein-Atome
                ax.scatter(
                    coords[:, 0],
                    coords[:, 1],
                    coords[:, 2],
                    s=20,
                    c="blue",
                    alpha=0.7,
                    label="Protein",
                )

                # Setze Achsen-Limits basierend auf Daten
                x_range = coords[:, 0].max() - coords[:, 0].min()
                y_range = coords[:, 1].max() - coords[:, 1].min()
                z_range = coords[:, 2].max() - coords[:, 2].min()

                max_range = max(x_range, y_range, z_range)
                center_x = (coords[:, 0].max() + coords[:, 0].min()) / 2
                center_y = (coords[:, 1].max() + coords[:, 1].min()) / 2
                center_z = (coords[:, 2].max() + coords[:, 2].min()) / 2

                ax.set_xlim(center_x - max_range / 2, center_x + max_range / 2)
                ax.set_ylim(center_y - max_range / 2, center_y + max_range / 2)
                ax.set_zlim(center_z - max_range / 2, center_z + max_range / 2)

                # Beschriftungen
                ax.set_xlabel("X (nm)")
                ax.set_ylabel("Y (nm)")
                ax.set_zlabel("Z (nm)")
                ax.set_title(f"Frame {i+1}")

                # Speichere Frame
                frame_path = img_dir / f"frame_{i:04d}.png"
                fig.savefig(frame_path, dpi=150, bbox_inches="tight")
                plt.close(fig)

                images.append(frame_path)
                frame_count += 1

            except Exception as e:
                LOGGER.warning("Fehler beim Rendern von Frame %d: %s", i, e)
                continue

        if not images:
            raise RuntimeError("Keine Frames erfolgreich gerendert")

        LOGGER.info("Gerenderte Frames: %d", frame_count)

        # Erstelle MP4 aus Frames
        try:
            clip = ImageSequenceClip([str(p) for p in images], fps=fps)
            clip.write_videofile(
                str(out_mp4), codec="libx264", audio=False, verbose=False, logger=None
            )
            LOGGER.info("MP4 erfolgreich erstellt: %s", out_mp4)

        except Exception as e:
            raise RuntimeError(f"Fehler beim Erstellen des MP4: {e}")

        return out_mp4

    finally:
        # Aufräumen
        try:
            shutil.rmtree(img_dir, ignore_errors=True)
        except Exception as e:
            LOGGER.warning("Fehler beim Aufräumen der temporären Dateien: %s", e)


__all__ = ["traj_to_mp4"]
