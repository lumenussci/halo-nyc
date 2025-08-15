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

#eval "$(conda shell.bash hook)"
cd /home/sean/sata/halo-nyc/stilt/stilt_code
conda activate stilt2
#Rscript r/run_stilt_nyc_hrrr.r --alt=${1} --file=${2} 
for alt in 50 150 250 500 750 1000 1500 2000 3000 5000 7500 10000;
do
    echo ${alt}
    Rscript r/run_stilt_nyc_hrrr.r --alt=${alt} --file=${1} 
done