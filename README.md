# halo-staaqs
 This project is a repository for scripts that are used to analyze the HALO-STAAQS methane data over NYC in July/August of 2023. This includes WRF-STILT modeling as well as inversions using footprints and boundary conditions from STILT.

## Contains:

## sample_wrfout_met.py
description: samples wrfout files to get the pressure variables needed to compare HALO observations to simulations from STILT

inputs: None

outputs: HDF5 files with the sampled variables

## sample_wrfout_ghg.py

purpose: samples wrfout files to pull out boundary conditions for STILT trajectories for the HALO data

input:
- domain: string like d01 or d02
- altitude: receptor altitude to distinguish file names
- flts: flights to read in like 20230726_F1 - folder names

output: 
- boundary point HDF5 file with wrfout samples
