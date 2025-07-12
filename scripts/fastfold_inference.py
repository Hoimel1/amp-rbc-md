#!/usr/bin/env python3
"""
FastFold Inference ohne MSA
"""
import argparse
import os
import sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='FastFold Inference ohne MSA')
    parser.add_argument('fasta_path', type=str, help='Pfad zur FASTA-Datei')
    parser.add_argument('output_dir', type=str, help='Ausgabeverzeichnis')
    parser.add_argument('--gpus', type=int, default=1, help='Anzahl GPUs')
    args = parser.parse_args()
    
    # Prüfe ob FastFold installiert ist
    try:
        from fastfold.inference import inference
    except ImportError:
        print("❌ FastFold ist nicht installiert. Bitte führen Sie 'conda activate folding-env' aus.")
        sys.exit(1)
    
    # Prüfe Parameter-Pfad
    params_path = os.path.expanduser('~/.fastfold/fastfold_params')
    if not os.path.exists(params_path):
        print(f"❌ FastFold-Parameter nicht gefunden in {params_path}")
        print("Bitte laden Sie die Parameter herunter:")
        print("mkdir -p ~/.fastfold && cd ~/.fastfold")
        print("wget -q https://github.com/hpcaitech/FastFold/releases/download/v0.2.0/fastfold_params.tar.gz")
        print("tar -xf fastfold_params.tar.gz")
        sys.exit(1)
    
    # Erstelle Ausgabeverzeichnis
    os.makedirs(args.output_dir, exist_ok=True)
    
    print(f"🚀 Starte FastFold-Inferenz...")
    print(f"   FASTA: {args.fasta_path}")
    print(f"   Output: {args.output_dir}")
    print(f"   GPUs: {args.gpus}")
    
    try:
        # Führe Inferenz aus
        inference(args.fasta_path, args.output_dir, gpus=args.gpus)
        print("✅ FastFold-Inferenz erfolgreich abgeschlossen!")
    except Exception as e:
        print(f"❌ Fehler bei FastFold-Inferenz: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main() 