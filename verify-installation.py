#!/usr/bin/env python3
"""
AMP-RBC-MD Installations-Verifikation
Prüft alle kritischen Module und Abhängigkeiten
"""

# SPDX-License-Identifier: MIT

"""AMP-RBC-MD Installations-Verifikation (modernisiert)

Nutzung:
    python verify-installation.py --engine fastfold

Unterstützte Engines:
    fastfold  (Default) – benötigt torch, openfold, fastfold, deepspeed
    colabfold – benötigt jax, colabfold
"""

import sys
import subprocess
from pathlib import Path
import argparse

def check_module(module_name, version_check=None):
    """Prüfe ob ein Modul verfügbar ist."""
    try:
        module = __import__(module_name)
        version = getattr(module, '__version__', 'unbekannt')
        print(f"✅ {module_name}: {version}")
        
        if version_check:
            version_check(module, version)
        return True
    except ImportError as e:
        print(f"❌ {module_name}: Import fehlgeschlagen - {e}")
        return False
    except Exception as e:
        print(f"⚠️  {module_name}: Fehler - {e}")
        return False

def check_jax_cuda():
    """Prüfe JAX CUDA-Unterstützung."""
    try:
        import jax
        import jaxlib
        
        # Prüfe CUDA-Geräte
        try:
            devices = jax.devices()
            gpu_devices = [d for d in devices if d.platform == 'gpu']
            if gpu_devices:
                print(f"✅ JAX CUDA: {len(gpu_devices)} GPU-Geräte gefunden")
                return True
            else:
                print("❌ JAX CUDA: Keine GPU-Devices gefunden")
                return False
        except Exception as e:
            print(f"❌ JAX CUDA: {e}")
            return False
    except ImportError:
        print("❌ JAX CUDA: JAX nicht installiert")
        return False

def check_linear_util():
    """Prüfe JAX linear_util Verfügbarkeit."""
    try:
        from jax._src import linear_util
        print("✅ JAX linear_util: Verfügbar")
        return True
    except ImportError:
        print("❌ JAX linear_util: Nicht verfügbar")
        return False

def check_colabfold():
    """Prüfe ColabFold-Funktionalität."""
    try:
        import colabfold
        print("✅ ColabFold: Import erfolgreich")
        
        # Prüfe batch-Funktion
        try:
            from colabfold import batch
            print("✅ ColabFold batch: Verfügbar")
        except ImportError:
            print("❌ ColabFold batch: Nicht verfügbar")
            
        # Prüfe download_alphafold_params
        try:
            from colabfold.download import download_alphafold_params
            print("✅ ColabFold download_alphafold_params: Verfügbar")
        except ImportError:
            print("❌ ColabFold download_alphafold_params: Nicht verfügbar")
            
        return True
    except Exception as e:
        print(f"❌ ColabFold: {e}")
        return False

def check_gromacs():
    """Prüfe GROMACS-Installation."""
    try:
        import gromacs
        print("✅ gromacs: Import erfolgreich")
        return True
    except ImportError:
        print("❌ gromacs: Import fehlgeschlagen - No module named 'gromacs'")
        return False

def check_amp_rbc_md():
    """Prüfe amp-rbc-md Installation."""
    try:
        import amp_rbc_md
        print("✅ amp-rbc-md: Import erfolgreich")
        
        # Prüfe CLI
        try:
            result = subprocess.run(['amp-rbc-md', '--help'], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                print("✅ amp-rbc-md CLI: Verfügbar")
                return True
            else:
                print("❌ amp-rbc-md CLI: Fehler")
                return False
        except Exception as e:
            print(f"❌ amp-rbc-md CLI: {e}")
            return False
    except ImportError as e:
        print(f"❌ amp-rbc-md: {e}")
        return False

def check_cuda_versions():
    """Prüfe CUDA-Versionskompatibilität."""
    try:
        import torch
        if torch.cuda.is_available():
            cuda_version = torch.version.cuda
            print(f"✅ PyTorch CUDA: {cuda_version}")
            
            # Prüfe GPU-Info
            gpu_count = torch.cuda.device_count()
            if gpu_count > 0:
                gpu_name = torch.cuda.get_device_name(0)
                print(f"   GPU: {gpu_name}")
                print(f"   GPU-Anzahl: {gpu_count}")
            return True
        else:
            print("❌ PyTorch CUDA: Nicht verfügbar")
            return False
    except Exception as e:
        print(f"❌ PyTorch CUDA: {e}")
        return False

def check_cudnn():
    """Prüfe cuDNN-Installation."""
    try:
        import nvidia.cudnn
        print("✅ nvidia-cudnn-cu12: Installiert")
        return True
    except ImportError:
        print("❌ nvidia-cudnn-cu12: Nicht installiert")
        return False

def main() -> None:  # noqa: D401
    """Starte die Verifikations-Checks."""

    parser = argparse.ArgumentParser(description="Verifiziere AMP-RBC-MD-Installation")
    parser.add_argument(
        "--engine",
        choices=["fastfold", "colabfold"],
        default="fastfold",
        help="Vorhersage-Engine, für die die Umgebung geprüft wird (default: fastfold)",
    )

    args = parser.parse_args()

    engine = args.engine.lower()

    print("=== AMP-RBC-MD INSTALLATIONS-VERIFIKATION ===")
    print(f"Engine: {engine}")

    # Python-Version
    print(f"Python Version: {sys.version}")
    print()

    # Allgemeine Basismodule
    print("=== GRUNDMODULE ===")
    pandas_ok = check_module("pandas")
    numpy_ok = check_module("numpy")
    gromacs_ok = check_gromacs()

    # Enginespezifische Prüfungen
    if engine == "fastfold":
        print("\n=== FASTFOLD / OPENFOLD STACK ===")
        torch_ok = check_module("torch")
        cuda_ok = check_cuda_versions()  # nutzt torch
        deepspeed_ok = check_module("deepspeed")
        openfold_ok = check_module("openfold")
        fastfold_ok = check_module("fastfold")

        engine_ok = all([torch_ok, cuda_ok, deepspeed_ok, openfold_ok, fastfold_ok])

    elif engine == "colabfold":
        print("\n=== COLABFOLD / JAX STACK ===")
        jax_ok = check_module("jax")
        jaxlib_ok = check_module("jaxlib")
        jax_cuda_ok = check_jax_cuda()
        linear_util_ok = check_linear_util()
        colabfold_ok = check_colabfold()

        engine_ok = all([jax_ok, jaxlib_ok, jax_cuda_ok, linear_util_ok, colabfold_ok])

    else:  # Fallback – sollte nie erreicht werden, argparse fängt falsch ab
        engine_ok = False

    # Projekt-Import & CLI-Check
    print("\n=== AMP-RBC-MD ===")
    amp_ok = check_amp_rbc_md()

    print("\n=== ZUSAMMENFASSUNG ===")
    all_ok = all([pandas_ok, numpy_ok, gromacs_ok, engine_ok, amp_ok])

    if all_ok:
        print("✅ Alle kritischen Module funktionieren – Installation sieht gut aus!")
    else:
        print("❌ Einige Module fehlen oder haben Probleme.")
        print("💡 Bitte README & TROUBLESHOOTING.md prüfen und fehlende Abhängigkeiten installieren.")

    # Tipps
    print("\n=== WEITERMACHEN ===")
    print("Beispiel-Dry-Run:")
    print("   amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --dry-run")
    print("Echte Simulation:")
    print("   amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --n-replica 1 --profile default -j 1")

if __name__ == "__main__":
    main() 