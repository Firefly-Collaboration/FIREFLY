#!/bin/bash
###################################################################################
#                         FIREFLY — Full Spectral Fitting                         #
###################################################################################
# <> Script: SBATCHidx_Fuji.sh
# <> Author:
#    - Samuel Helps ~ <samuel.helps.sh__at__gmail.com>
# ---------------------------------------------------------------------------------
# Institute of Cosmology and Gravitation - University of Portsmouth, Portsmouth, UK
# ---------------------------------------------------------------------------------

#SBATCH --job-name=idesifirefly
#SBATCH --partition=sciama3-5.q
#SBATCH --mail-type=ALL
#SBATCH --mail-user=up2013158@myport.ac.uk
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
##SBATCH --mem-per-cpu=4G
#SBATCH --output=../slurm_output/firefly_fit_%j.out
#SBATCH --error=../slurm_output/firefly_fit_%j.err
#SBATCH --time=06:00:00

# Usage:
# sbatch SBATCHidx_Fuji.sh DESI_EDR_10000-20000.fits [ /full/path/to/EDR_master_10000-20000.csv ]

set -o errexit
set -o pipefail
set -o nounset

# Check if input file is provided
if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "Usage: sbatch $0 <input_fits_basename> [master_csv_path]"
    exit 1
fi

INPUT_FILE="$1"           # e.g. DESI_EDR_10000-20000.fits
USER_MASTER_CSV="${2:-}"  # optional explicit master csv path

# Directory of the SBATCH script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# Repo root (firefly)
BASE_DIR="$(dirname $(dirname $(dirname $(dirname "$SCRIPT_DIR"))))"

# Extract numeric range (e.g. 10000-20000) from filename
RANGE=$(echo "$INPUT_FILE" | grep -o '[0-9]\+-[0-9]\+' || true)
if [ -z "$RANGE" ]; then
    echo "ERROR: Could not find a range like X-Y in input file name: $INPUT_FILE"
    exit 1
fi
START_IDX=$(echo "$RANGE" | cut -d'-' -f1)
END_IDX=$(echo "$RANGE" | cut -d'-' -f2)

# Default master csv path if not supplied by terminal argument
if [ -z "$USER_MASTER_CSV" ]; then  # for CSV-only runs (example path)
    MASTER_FILE="${USER_MASTER_CSV:-$BASE_DIR/Results/DESI_EDR_VAC/master_files/EDR_master_${RANGE}.csv}"
else                 
    MASTER_FILE="$USER_MASTER_CSV"
fi

if [ ! -f "$MASTER_FILE" ]; then
    echo "Error: Master file not found at: $MASTER_FILE"
    exit 1
fi

echo "Input FITS basename: $INPUT_FILE"
echo "Using master CSV: $MASTER_FILE"
echo "Index range parsed: ${START_IDX}-${END_IDX}"

# Load Anaconda and activate base environment
module load anaconda3/2024.02
source activate base

# Set environment variables
export DESI_INPUT_FILE="$INPUT_FILE"
export PYTHONPATH="$BASE_DIR/python:$PYTHONPATH"

# Navigate to base directory
cd "$BASE_DIR"

# CONFIG: how many parallel jobs at once
# default to SLURM ntasks if set, otherwise 32
MAX_JOBS="${SLURM_NTASKS:-32}"
echo "Will run up to $MAX_JOBS parallel tasks."

# Build missing local indices list once (space separated)
# Uses Python to be robust with header/column names; returns local indices (global-START)
MISSING_LIST=$(python3 - <<PY
import sys, pandas as pd
master = r"${MASTER_FILE}"
start = int(${START_IDX})
end = int(${END_IDX})
try:
    df = pd.read_csv(master)
except Exception as e:
    print("ERROR_READING_CSV:"+str(e))
    sys.exit(2)

# find index-like column: 'spec'|'spectrum'|'index' heuristics, else first column
candidates = [c for c in df.columns if ('spec' in c.lower() or 'spectrum' in c.lower() or 'index' in c.lower())]
idx_col = None
for c in candidates:
    try:
        vals = pd.to_numeric(df[c], errors='coerce').dropna().astype(int).unique().tolist()
        if vals:
            idx_col = c
            break
    except Exception:
        continue
if idx_col is None:
    try:
        vals = pd.to_numeric(df.iloc[:,0], errors='coerce').dropna().astype(int).unique().tolist()
        idx_col = df.columns[0]
    except Exception as e:
        print("ERROR_FINDING_INDEX_COL:"+str(e))
        sys.exit(3)

existing_local = set(pd.to_numeric(df[idx_col], errors='coerce').dropna().astype(int).tolist())
# master CSV contains local indices 0..N-1 
full_local = set(range(0, max(existing_local) + 1))  # not used, we'll use start..end mapping
# compute full set by mapping global->local
missing = []
for g in range(start, end+1):
    l = g - start
    if l not in existing_local:
        missing.append(l)
if not missing:
    print("NONE")
else:
    print(" ".join(str(x) for x in sorted(missing)))
PY
)

# Check for errors in MISSING_LIST retrieval
if [[ "$MISSING_LIST" == ERROR_* ]]; then
    echo "Error determining missing indices: $MISSING_LIST"
    exit 1
fi

# If no missing indices, exit
if [ "$MISSING_LIST" == "NONE" ] || [ -z "$MISSING_LIST" ]; then
    echo "No missing indices found in master CSV. Nothing to retry."
    exit 0
fi

echo "Missing local indices: $MISSING_LIST"

# Semaphore using a FIFO (POSIX): limit concurrency to MAX_JOBS
sem=/tmp/$$.sem
mkfifo "$sem"
# Open file descriptor 3 for reading the fifo
exec 3<> "$sem"
rm "$sem"

# Populate the semaphore with MAX_JOBS tokens
for ((i=0;i<MAX_JOBS;i++)); do
    printf '%s\n' "token" >&3
done

# Launch tasks in background; each consumes a token and returns it on completion
run_task() {
    local local_idx="$1"
    # consume token (read a line)
    read -r -u 3 token || { echo "Failed to acquire sem token"; exit 1; }

    # small randomised stagger to reduce simultaneous writes (0-2s)
    sleep_time=$(awk -v seed="$local_idx" 'BEGIN{srand(seed); print rand()*2}')
    sleep "$sleep_time"

    echo "Starting local index $local_idx (global $((local_idx + START_IDX))) on $(hostname) (stagger ${sleep_time}s)"
    # call the unchanged firefly runner
    srun -n1 --cpus-per-task=1 --mem=4G python3 Launch/DESI/SCIAMA/core/firefly_DESI_EDR.py "$local_idx" "$((local_idx + 1))"
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "Warning: run for local index $local_idx exited with code $rc"
    else
        echo "Completed local index $local_idx (global $((local_idx + START_IDX)))"
    fi

    # Return token to semaphore
    printf '%s\n' "$token" >&3
}

# Iterate missing indices and dispatch
for li in $MISSING_LIST; do
    # start run_task in background
    run_task "$li" &
done

# Wait for all background children to finish
wait

# Close FD 3
exec 3>&-
exec 3<&-

echo "All parallel runs finished for $INPUT_FILE"
