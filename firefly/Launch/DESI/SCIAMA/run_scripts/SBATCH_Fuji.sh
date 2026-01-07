#!/bin/bash
###################################################################################
#                         FIREFLY — Full Spectral Fitting                         #
###################################################################################
# <> Script: SBATCH_Fuji.sh
# <> Author:
#    - Samuel Helps ~ <samuel.helps.sh__at__gmail.com>
# ---------------------------------------------------------------------------------
# Institute of Cosmology and Gravitation - University of Portsmouth, Portsmouth, UK
# ---------------------------------------------------------------------------------

#SBATCH --job-name=desifirefly
#SBATCH --partition=sciama3-5.q
#SBATCH --mail-type=ALL
#SBATCH --mail-user=up2013158@myport.ac.uk
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --output=../slurm_output/firefly_fit_%j.out
#SBATCH --error=../slurm_output/firefly_fit_%j.err

# Usage:
# sbatch SBATCH_Fuji.sh DESI_EDR_10000-20000.fits [ /full/path/to/EDR_master_10000-20000.csv ]

# Check if input file is provided
if [ $# -ne 1 ]; then
    echo "Usage: sbatch SBATCH: <input_file>"
    exit 1
fi

INPUT_FILE=$1

# Directory of the SBATCH script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# Repo root (firefly)
BASE_DIR="$(dirname $(dirname $(dirname $(dirname "$SCRIPT_DIR"))))"

# Extract range from filename (e.g., DESI_EDR_5000-10000.fits -> 5000 10000)
RANGE=$(echo $INPUT_FILE | grep -o '[0-9]*-[0-9]*' | head -n1)
START_IDX=$(echo $RANGE | cut -d'-' -f1)
END_IDX=$(echo $RANGE | cut -d'-' -f2)

# Load Anaconda and activate base environment
module load anaconda3/2024.02
source activate base

# Set environment variables
export PYTHONPATH="$BASE_DIR/python:$PYTHONPATH"
export DESI_INPUT_FILE="$INPUT_FILE"

# Navigate to base directory
cd "$BASE_DIR"

# Always use 0-5000 for internal processing
spectra_start=0
spectra_end=10000

# Export the file range for output directory naming
export DESI_FILE_RANGE=$RANGE

# Number of spectra divided by number of CPU/tasks
chunk_size=313
declare -A running_pids

# Loop over spectra in each chunk
for ((i=spectra_start; i<spectra_end; i+=chunk_size)); do
    chunk_end=$((i + chunk_size))
    if [ $chunk_end -gt $spectra_end ]; then
        chunk_end=$spectra_end
    fi

    # Wait if too many tasks are running
    while [ ${#running_pids[@]} -ge $SLURM_NTASKS ]; do
        for pid in "${!running_pids[@]}"; do
            if ! kill -0 $pid 2>/dev/null; then
                wait $pid
                unset running_pids[$pid]
            fi
        done
        sleep 2
    done

    echo "Processing spectra $i to $chunk_end"
    srun -n1 --mem=10G python3 Launch/DESI/SCIAMA/core/firefly_DESI_EDR.py $i $chunk_end &
    pid=$!
    running_pids[$pid]="$i-$chunk_end"
done

# Wait for remaining tasks
wait

echo "All spectra processed for $INPUT_FILE"