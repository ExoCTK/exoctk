#!/bin/bash
#SBATCH --job-name=HJ_Ray_logPc
#SBATCH --error=err
#SBATCH --time=1-12:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
cd /home/mrline/EXOPLANETS/CHIMERA_TRANSMISSION_CK/OPERATIONAL_HOT_JUPITER_RAY_HZ_GRAY_CLD/
module load openmpi/1.8.5-gnu-stock
module load multinest/3.10
module load anaconda2/4.2.0
mpirun -np 64 python call_pymultinest.py > status

