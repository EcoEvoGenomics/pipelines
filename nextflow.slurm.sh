#!/bin/bash

# CONFIGURABLE: SLURM SETTIGNS
#SBATCH --job-name=submit-pipeline
#SBATCH --output=SLURM-%j-%x.out
#SBATCH --error=SLURM-%j-%x.err
#SBATCH --account=nn10082k
#SBATCH --tasks=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00

# HOUSEKEEPING
set -o errexit
set -o nounset
module --quiet purge

# CONFIGURABLE: NEXTFLOW
module load Miniconda3/22.11.1-1
source ${EBROOTMINICONDA3}/bin/activate
conda activate /cluster/projects/nn10082k/conda_group/Nextflow25.04.6

# CONFIGURABLE: PIPELINE
pipeline=population_structure

# CONFIGURABLE: NEXTFLOW OPTIONS
nextflow -log .nextflow/nextflow.log run pipelines/${pipeline}.nf \
    -params-file nextflow.yaml \
    -profile saga
