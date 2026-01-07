#!/bin/bash
###################################################################################
#                         FIREFLY — Full Spectral Fitting                         #
###################################################################################
# <> Script: SBATCH_Iron.sh
# <> Author:
#    - Kieran Graham ~ <kieran.graham__at__port.ac.uk>
# <> Contributors:
#	 - Samuel Helps ~ <samuel.helps.sh__at__gmail.com>
# ---------------------------------------------------------------------------------
# Institute of Cosmology and Gravitation - University of Portsmouth, Portsmouth, UK
# ---------------------------------------------------------------------------------

#SBATCH --job-name=IRON_AIO
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --account=desi
#SBATCH --output=../slurm_output/firefly_fit_%j.out
#SBATCH --error=../slurm_output/firefly_fit_%j.err

# --- Job Configuration ---
# Define the start and end indices for this specific job run
# This makes it easy to submit different batches
START_INDEX=2600000
END_INDEX=2700000
SCRIPT_PATH="/global/cfs/cdirs/desi/users/helpss/FIREFLY/firefly/Launch/DESI/NERSC/run_scripts/SBATCH_Iron.sh" 
echo "--- Job Started at $(date) ---"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
echo "Processing rows from $START_INDEX to $END_INDEX"
echo "-----------------------------------"

# --- Environment Setup ---
echo "Setting up DESI environment..."

# STEP 1: Source the main environment file to define $DESI_ROOT
source /global/common/software/desi/desi_environment.sh

# STEP 2: Source the specific versioned environment
source $DESI_ROOT/desi_environment.sh 25.3

# STEP 3: Explicitly activate the conda environment
echo "Activating conda environment: $DESICONDA"
conda activate $DESICONDA

echo "Environment setup complete."
echo "-----------------------------------"

# --- Execution ---
echo "Starting Python script..."

# Run the python script, passing the start and end indices as command-line arguments
# The output is piped through grep to filter out the benign "Frame object" warning
python $SCRIPT_PATH --start $START_INDEX --end $END_INDEX 2>&1 | grep -v "Frame object is constructed without resolution data"

echo "-----------------------------------"
echo "--- Job Finished at $(date) ---"