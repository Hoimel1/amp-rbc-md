#!/usr/bin/env python3
"""
AMP-RBC-MD Installations-Verifikation
Prüft alle kritischen Module und Abhängigkeiten
"""

import sys
import subprocess
from pathlib import Path

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

def main():
    """Hauptfunktion für Verifikation."""
    print("=== AMP-RBC-MD INSTALLATIONS-VERIFIKATION ===")
    
    # Python-Version
    print(f"Python Version: {sys.version}")
    print()
    
    # Kritische Module
    print("=== KRITISCHE MODULE ===")
    jax_ok = check_module("jax")
    jaxlib_ok = check_module("jaxlib")
    torch_ok = check_module("torch")
    colabfold_ok = check_colabfold()
    pandas_ok = check_module("pandas")
    numpy_ok = check_module("numpy")
    gromacs_ok = check_gromacs()
    
    print()
    print("=== CUDA-UNTERSTÜTZUNG ===")
    cuda_ok = check_cuda_versions()
    jax_cuda_ok = check_jax_cuda()
    
    print()
    print("=== JAX LINEAR_UTIL ===")
    linear_util_ok = check_linear_util()
    
    print()
    print("=== COLABFOLD ===")
    colabfold_ok = check_colabfold()
    
    print()
    print("=== CUDNN ===")
    cudnn_ok = check_cudnn()
    
    print()
    print("=== AMP-RBC-MD ===")
    amp_ok = check_amp_rbc_md()
    
    print()
    print("=== ZUSAMMENFASSUNG ===")
    
    all_ok = all([
        jax_ok, jaxlib_ok, torch_ok, pandas_ok, numpy_ok, 
        cuda_ok, cudnn_ok, amp_ok
    ])
    
    if all_ok:
        print("✅ Alle kritischen Module funktionieren")
    else:
        print("❌ Einige Module fehlen oder haben Probleme")
        print("💡 Führen Sie 'bash setup-fixed.sh' aus")
    
    print()
    print("=== NÄCHSTE SCHRITTE ===")
    print("1. Testen Sie eine Simulation:")
    print("   amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --dry-run")
    print("2. Führen Sie eine echte Simulation aus:")
    print("   amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --n-replica 1 --profile default -j 1")

if __name__ == "__main__":
    main() 