#!/bin/bash
#SBATCH  -n 1
#SBATCH --cpus-per-task=106
#SBATCH -t 1-
#SBATCH --mem=80g
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=erisembe@email.unc.edu

module load r/4.4.0

Rscript var_comp.R --args --mode=SLURM
