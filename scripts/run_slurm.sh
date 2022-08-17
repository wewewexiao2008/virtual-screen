#!/bin/bash

#SBATCH -J virtual_screening_csz_extracting
#SBATCH –p all
#SBATCH -n 200
#SBATCH -o out.%j
#SBATCH -e err.%j

srun ./scripts/job-extract.sh
