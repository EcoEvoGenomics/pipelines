#!/bin/bash

# CONFIGURE: SLURM SETTIGNS
#SBATCH --job-name=submit-pipeline
#SBATCH --output=SLURM-%j-%x.out
#SBATCH --error=SLURM-%j-%x.err
#SBATCH --account=nn10082k
#SBATCH --tasks=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00

# CHOOSE: PIPELINE
pipeline=population_structure

# CONFIGURE: NEXTFLOW, CONDA, SINGULARITY (already enabled on Saga)
module load Miniconda3/23.10.0-1
source ${EBROOTMINICONDA3}/bin/activate
conda activate /cluster/projects/nn10082k/conda_group/Nextflow25.04.6
conda clean --yes --index-cache

# CONFIGURE: NEXTFLOW OPTIONS
nextflow -log .nextflow/nextflow.log run pipelines/${pipeline}.nf \
    -params-file nextflow.yaml \
    -profile saga
