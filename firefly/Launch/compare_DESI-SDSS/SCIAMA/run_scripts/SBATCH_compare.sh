#!/bin/bash
###################################################################################
#                         FIREFLY — Full Spectral Fitting                         #
###################################################################################
# <> Script: SBATCH_compare.sh
# <> Author:
#    - Samuel Helps ~ <samuel.helps.sh__at__gmail.com>
# ---------------------------------------------------------------------------------
# Institute of Cosmology and Gravitation - University of Portsmouth, Portsmouth, UK
# ---------------------------------------------------------------------------------

#SBATCH --job-name=firefly_compare
#SBATCH --partition=sciama3.q
#SBATCH --mail-type=ALL
#SBATCH --mail-user=up2013158@myport.ac.uk
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --output=../slurm_output/firefly_fit_%j.out
#SBATCH --error=../slurm_output/firefly_fit_%j.err

# Usage: sbatch SBATCH_compare.sh <sample_name> [cut]
if [ $# -lt 1 ]; then
    echo "Usage: sbatch SBATCH_compare.sh <sample_name> [cut]"
    exit 1
fi

SAMPLE=$1
CUTMODE=""
if [ $# -eq 2 ] && [ "$2" == "cut" ]; then
    CUTMODE="cut"
fi

# Directory of the SBATCH script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Go up 4 levels to the repo root (firefly)
BASE_DIR="$(dirname $(dirname $(dirname $(dirname "$SCRIPT_DIR"))))"

INPUT_DIR="$BASE_DIR/Data/DESI/DESI-SDSS_samples"
DESI_FILE="$INPUT_DIR/${SAMPLE}_desi_dr1.fits"
SDSS_FILE="$INPUT_DIR/${SAMPLE}_sdss_dr16.fits"

if [ ! -f "$DESI_FILE" ]; then
    echo "DESI input file not found: $DESI_FILE"
    exit 1
fi
if [ ! -f "$SDSS_FILE" ]; then
    echo "SDSS input file not found: $SDSS_FILE"
    exit 1
fi

# Load Anaconda and activate base environment
module load anaconda3/2024.02
source activate base

# Set PYTHONPATH and environment variables
export PYTHONPATH="$BASE_DIR/python:$PYTHONPATH"
export DESI_INPUT_FILE="$DESI_FILE"
export SDSS_INPUT_FILE="$SDSS_FILE"
export SAMPLE_NAME="$SAMPLE"
export CUT_MODE="$CUTMODE"

# Navigate to base directory
cd "$BASE_DIR"

# Define spectral range (0-indexed), Adjust if needed
spectra_start=0
spectra_end=117   # assume 117 spectra (indices 0 to 116) in 'greatwall' sample

chunk_size=4

for ((i=spectra_start; i<spectra_end; i+=chunk_size)); do
    chunk_end=$((i + chunk_size))
    if [ $chunk_end -gt $spectra_end ]; then
        chunk_end=$spectra_end
    fi

    # Reverse order within each chunk (found to have better memory usage on Sciama nodes)
    for ((j=chunk_end-1; j>=i; j--)); do
        echo "Processing spectrum index $j"
        srun -n1 --mem=10G python3 Launch/compare_DESI-SDSS/SCIAMA/core/firefly_DESI-SDSS.py $j $((j + 1)) &
    done
done

wait

echo "All spectra processed for sample $SAMPLE"
