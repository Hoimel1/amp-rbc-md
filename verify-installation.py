#!/usr/bin/env python3
"""
Verifikations-Skript für amp-rbc-md Installation
Testet alle kritischen Komponenten und Versionen
"""

import sys
import subprocess
import importlib

def test_import(module_name, version_attr='__version__'):
    """Teste Import und Version eines Moduls"""
    try:
        module = importlib.import_module(module_name)
        version = getattr(module, version_attr, 'unbekannt')
        print(f"✅ {module_name}: {version}")
        return True
    except ImportError as e:
        print(f"❌ {module_name}: Import fehlgeschlagen - {e}")
        return False
    except Exception as e:
        print(f"⚠️  {module_name}: Fehler - {e}")
        return False

def test_cuda_support():
    """Teste CUDA-Unterstützung"""
    print("\n=== CUDA-UNTERSTÜTZUNG ===")
    
    # PyTorch CUDA
    try:
        import torch
        if torch.cuda.is_available():
            print(f"✅ PyTorch CUDA: {torch.version.cuda}")
            print(f"   GPU: {torch.cuda.get_device_name(0)}")
            print(f"   GPU-Anzahl: {torch.cuda.device_count()}")
        else:
            print("❌ PyTorch CUDA: Nicht verfügbar")
    except Exception as e:
        print(f"❌ PyTorch CUDA Test fehlgeschlagen: {e}")
    
    # JAX CUDA
    try:
        import jax
        devices = jax.devices()
        gpu_devices = [d for d in devices if d.platform == 'gpu']
        if gpu_devices:
            print(f"✅ JAX CUDA: {len(gpu_devices)} GPU(s) verfügbar")
            print(f"   Devices: {devices}")
        else:
            print("❌ JAX CUDA: Keine GPU-Devices gefunden")
    except Exception as e:
        print(f"❌ JAX CUDA Test fehlgeschlagen: {e}")

def test_linear_util():
    """Teste JAX linear_util Verfügbarkeit"""
    print("\n=== JAX LINEAR_UTIL ===")
    try:
        import jax
        if hasattr(jax, 'linear_util'):
            print("✅ JAX linear_util: Verfügbar")
            return True
        else:
            print("❌ JAX linear_util: Nicht verfügbar")
            return False
    except Exception as e:
        print(f"❌ JAX linear_util Test fehlgeschlagen: {e}")
        return False

def test_colabfold():
    """Teste ColabFold-Funktionalität"""
    print("\n=== COLABFOLD ===")
    try:
        import colabfold
        print("✅ ColabFold: Import erfolgreich")
        
        # Teste batch-Funktion
        try:
            if hasattr(colabfold, 'batch'):
                print("✅ ColabFold batch: Verfügbar")
            else:
                print("❌ ColabFold batch: Nicht verfügbar")
        except Exception as e:
            print(f"⚠️  ColabFold batch Test: {e}")
        
        # Teste download_alphafold_params
        try:
            if hasattr(colabfold, 'download_alphafold_params'):
                print("✅ ColabFold download_alphafold_params: Verfügbar")
            else:
                print("❌ ColabFold download_alphafold_params: Nicht verfügbar")
        except Exception as e:
            print(f"⚠️  ColabFold download_alphafold_params Test: {e}")
            
        return True
    except Exception as e:
        print(f"❌ ColabFold Test fehlgeschlagen: {e}")
        return False

def test_cudnn():
    """Teste cuDNN Installation"""
    print("\n=== CUDNN ===")
    try:
        result = subprocess.run(['pip', 'list'], capture_output=True, text=True)
        if 'nvidia-cudnn-cu12' in result.stdout:
            print("✅ nvidia-cudnn-cu12: Installiert")
            return True
        else:
            print("❌ nvidia-cudnn-cu12: Nicht gefunden")
            return False
    except Exception as e:
        print(f"❌ cuDNN Test fehlgeschlagen: {e}")
        return False

def test_amp_rbc_md():
    """Teste amp-rbc-md Installation"""
    print("\n=== AMP-RBC-MD ===")
    try:
        import amp_rbc_md
        print("✅ amp-rbc-md: Import erfolgreich")
        
        # Teste CLI
        result = subprocess.run(['amp-rbc-md', '--help'], capture_output=True, text=True)
        if result.returncode == 0:
            print("✅ amp-rbc-md CLI: Verfügbar")
            return True
        else:
            print("❌ amp-rbc-md CLI: Fehler")
            return False
    except Exception as e:
        print(f"❌ amp-rbc-md Test fehlgeschlagen: {e}")
        return False

def main():
    """Hauptfunktion"""
    print("=== AMP-RBC-MD INSTALLATIONS-VERIFIKATION ===")
    print(f"Python Version: {sys.version}")
    
    # Teste kritische Module
    print("\n=== KRITISCHE MODULE ===")
    modules = [
        ('jax', '__version__'),
        ('jaxlib', '__version__'),
        ('torch', '__version__'),
        ('colabfold', None),
        ('pandas', '__version__'),
        ('numpy', '__version__'),
        ('gromacs', None),  # Falls verfügbar
    ]
    
    all_modules_ok = True
    for module_name, version_attr in modules:
        if not test_import(module_name, version_attr):
            all_modules_ok = False
    
    # Spezielle Tests
    test_cuda_support()
    test_linear_util()
    test_colabfold()
    test_cudnn()
    test_amp_rbc_md()
    
    # Zusammenfassung
    print("\n=== ZUSAMMENFASSUNG ===")
    if all_modules_ok:
        print("✅ Alle kritischen Module sind verfügbar")
        print("🎉 Installation erfolgreich!")
    else:
        print("❌ Einige Module fehlen oder haben Probleme")
        print("💡 Führen Sie 'bash setup-clean-environment.sh' aus")
    
    print("\n=== NÄCHSTE SCHRITTE ===")
    print("1. Testen Sie eine Simulation:")
    print("   amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --dry-run")
    print("2. Führen Sie eine echte Simulation aus:")
    print("   amp-rbc-md --seq AAHHIIGGLFSAGKAIHRLIRRRRR --n-replica 1 --profile default -j 1")

if __name__ == "__main__":
    main() 