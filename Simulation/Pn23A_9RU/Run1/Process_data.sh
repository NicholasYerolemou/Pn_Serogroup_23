#!/bin/sh
#SBATCH --account=compsci
#SBATCH --partition=ada
#SBATCH --nodes=1 --ntasks=8
#SBATCH --time=02:00:00
#SBATCH --job-name remove_water
#SBATCH --mail-user=YRLNIC001@myuct.ac.za
#SBATCH --mail-type=ALL

module load software/vmd-1.9.3
vmd -psf ../Pn23A_9RU_Min_H2O_Na.psf -e Process_output.tcl #removes the water and creates output pdb and dcd files
	
