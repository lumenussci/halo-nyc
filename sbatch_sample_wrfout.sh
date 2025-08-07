#!/bin/bash
#SBATCH --job-name=sample_wrfout
#SBATCH --ntasks=48
#SBATCH -N 1
#SBATCH --time=18:00:00
#SBATCH -o output.%j
#SBATCH -e error.%j
#SBATCH -p skx
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user="sean@belumenus.com"

eval "$(conda shell.bash hook)"
cd /scratch/07351/tg866507/halo-staaqs/ 
conda activate analysis 
python halo_boundary_sampler.py ${1} ${2} 20230726_F1 20230726_F2 20230728_F1 20230728_F2 20230805_F1 20230809_F1 
