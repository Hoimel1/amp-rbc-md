#!/usr/bin/env bash
# Download-Skript für experimentelle PDB-mmCIF-Daten (FastFold/OpenFold)
# ---------------------------------------------------------------
# Verwendet rsync-Mirror des PDBe (EMBL-EBI).
# Legt die Daten unter $OPENFOLD_DATA/pdb_mmcif/mmcif_files ab.
# Aufruf:  bash scripts/download_pdb_mmcif.sh  [/pfad/zu/alphafold_dbs]
# ---------------------------------------------------------------
set -euo pipefail

ROOT=${1:-$HOME/alphafold_dbs}
export OPENFOLD_DATA="$ROOT"

MMCIF_DIR="$OPENFOLD_DATA/pdb_mmcif/mmcif_files"
mkdir -p "$MMCIF_DIR"

echo "[INFO] Lade PDB-mmCIF-Daten in $MMCIF_DIR … (kann Stunden dauern)"
rsync -rlpt -v -z --delete \
  rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/mmCIF/ \
  "$MMCIF_DIR"

# SEQRES-Datei (klein, nützlich für Alignments)
SEQFILE="$OPENFOLD_DATA/pdb_mmcif/pdb_seqres.txt"
mkdir -p "$(dirname "$SEQFILE")"
wget -q -O "$SEQFILE" \
  https://ftp.ebi.ac.uk/pub/databases/pdb/derived_data/pdb_seqres.txt

echo "[DONE] PDB-mmCIF-Mirror vollständig. OPENFOLD_DATA=$OPENFOLD_DATA" 