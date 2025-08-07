#!/bin/bash
#SBATCH --job-name=run_stilt
#SBATCH --ntasks=48
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH -o output.%j
#SBATCH -e error.%j
#SBATCH -p skx
#SBATCH --exclusive
#SBATCH --mail-type=FAIL
#SBATCH --mail-user="sean@belumenus.com"

eval "$(conda shell.bash hook)"
cd /scratch/07351/tg866507/halo 
conda activate stilt
#Rscript r/run_stilt_nyc_hrrr.r --alt=${1} --file=${2} 
Rscript r/run_stilt_nyc_wrfghgd02.r --alt=${1} --file=${2} 
