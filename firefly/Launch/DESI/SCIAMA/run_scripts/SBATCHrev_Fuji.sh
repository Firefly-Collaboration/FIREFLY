#!/bin/bash
###################################################################################
#                         FIREFLY — Full Spectral Fitting                         #
###################################################################################
# <> Script: SBATCHrev_Fuji.sh
# <> Author:
#    - Samuel Helps ~ <samuel.helps.sh__at__gmail.com>
# ---------------------------------------------------------------------------------
# Institute of Cosmology and Gravitation - University of Portsmouth, Portsmouth, UK
# ---------------------------------------------------------------------------------

#SBATCH --job-name=Rdesifirefly
#SBATCH --partition=sciama3-5.q
#SBATCH --mail-type=ALL
#SBATCH --mail-user=up2013158@myport.ac.uk
#SBATCH --ntasks=32 
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --output=../slurm_output/firefly_fit_%j.out
#SBATCH --error=../slurm_output/firefly_fit_%j.err

# Usage:
# sbatch SBATCHrev_Fuji.sh DESI_EDR_10000-20000.fits [ /full/path/to/EDR_master_10000-20000.csv ]

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

# Navigate to the base directory
cd "$BASE_DIR"

# Use 0-10000 for efficient DESI galaxies batch processing (in this pipline with the full plots)
spectra_start=0
spectra_end=10000

# Export the file range for output directory naming
export DESI_FILE_RANGE=$RANGE

# Number of spectra divided by number of CPU/tasks
chunk_size=313

# Loop over spectra in each chunk
for ((i=spectra_start; i<spectra_end; i+=chunk_size)); do
    chunk_end=$((i + chunk_size))
    if [ $chunk_end -gt $spectra_end ]; then
        chunk_end=$spectra_end
    fi

    # Reverse the order of indices within the chunk
    for ((j=chunk_end-1; j>=i; j--)); do
        echo "Processing spectrum $j"
        srun -n1 --mem=10G python3 Launch/DESI/SCIAMA/core/firefly_DESI_EDR.py $j $((j + 1)) &
    done
done

# Wait for all tasks to complete
wait

echo "All spectra processed for $INPUT_FILE"