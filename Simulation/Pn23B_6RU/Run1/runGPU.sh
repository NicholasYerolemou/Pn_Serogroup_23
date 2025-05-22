#!/bin/sh
#SBATCH --account=gpumk --partition=gpumk
###for hpc node s  --account=compsci --partition=ada etc etc
##SBATCH --nodelist=srvcntgpu008
#SBATCH --nodes=1 --ntasks=32 --gres=gpu:4
#SBATCH --time=200:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name="Pn23B_6RU_Run1"
##SBATCH --dependency=afterok:
#SBATCH --mail-user=YRLNIC001@myuct.ac.za
#SBATCH --mail-type=ALL


cd /scratch/yrlnic001/Pn23B_6RU/Run1
#When running CUDA NAMD always add +idlepoll to the command line. This is needed to poll the GPU for results rather than sleeping while idle.
export LD_LIBRARY_PATH=/opt/exp_soft/NAMD_2.14_Linux-x86_64-multicore-CUDA/:$LD_LIBRARY_PATH
/opt/exp_soft/NAMD_2.14_Linux-x86_64-multicore-CUDA/namd2 +p32 +idlepoll +noAnytimeMigration +setcpuaffinity +isomalloc_sync run_1.conf > run_1.log





