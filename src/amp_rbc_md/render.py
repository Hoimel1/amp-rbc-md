from __future__ import annotations

import shutil
from pathlib import Path
from typing import Tuple

from .utils import LOGGER

try:
    import MDAnalysis as mda  # type: ignore
    from moviepy.editor import ImageSequenceClip  # type: ignore
    import matplotlib.pyplot as plt  # type: ignore
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
    """Rendere ein kurzes MP4 der Trajektorie.

    Bei fehlender Abhängigkeit entsteht eine Dummy-Datei.
    """
    traj = Path(traj)
    out_mp4 = Path(out_mp4)

    if mda is None or not traj.exists():
        LOGGER.warning("Rendering-Backends fehlen, schreibe Dummy-MP4")
        out_mp4.write_bytes(b"dummy mp4")
        return out_mp4

    if gro is None:
        gro_guess = traj.with_suffix(".gro")
        if not gro_guess.exists():
            raise FileNotFoundError("Gro-Datei benötigt für Rendering")
        gro = gro_guess

    u = mda.Universe(str(gro), str(traj))  # type: ignore[arg-type]
    protein = u.select_atoms("protein or name CA")

    img_dir = out_mp4.parent / "frames_tmp"
    img_dir.mkdir(parents=True, exist_ok=True)
    images = []
    for i, ts in enumerate(u.trajectory):
        if i >= MAX_FRAMES:
            break
        fig = plt.figure(figsize=(4, 4))
        ax = fig.add_subplot(111, projection="3d")  # type: ignore[attr-defined]
        coords = protein.positions / 10  # Å → nm
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=5, c="blue")
        ax.set_xlim(-5, 5)
        ax.set_ylim(-5, 5)
        ax.set_zlim(-5, 5)
        ax.axis("off")
        frame_path = img_dir / f"frame_{i:04d}.png"
        fig.savefig(frame_path, dpi=100, bbox_inches="tight")
        plt.close(fig)
        images.append(frame_path)
    if not images:
        out_mp4.write_bytes(b"dummy mp4")
        return out_mp4

    clip = ImageSequenceClip([str(p) for p in images], fps=fps)
    clip.write_videofile(str(out_mp4), codec="libx264", audio=False, verbose=False, logger=None)

    # Aufräumen
    shutil.rmtree(img_dir, ignore_errors=True)
    return out_mp4


__all__ = ["traj_to_mp4"] 