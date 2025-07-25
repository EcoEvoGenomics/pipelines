#!/bin/bash

# ADMIN
#SBATCH --job-name=population_structure
#SBATCH --output=SLURM-%j-%x.out
#SBATCH --error=SLURM-%j-%x.err
#SBATCH --account=nn10082k

# RESOURCE ALLOCATION
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00

# USER INPUT: Define absolute path to pipelines repository
repodir=

# USER CONFIGURATION: Enable Nextflow (if necessary)
module --quiet purge
module load Miniconda3/22.11.1-1
source ${EBROOTMINICONDA3}/bin/activate
conda activate /cluster/projects/nn10082k/conda_group/Nextflow25.04.6

# Begin work
set -o errexit
set -o nounset
cd ${repodir}

nextflow -log .nextflow/nextflow.log \
    run ./pipelines/population_structure.nf \
    -params-file pipelines/population_structure.param.yaml \
    -profile saga

# End work
