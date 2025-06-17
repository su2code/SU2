#!/bin/bash


#SBATCH --job-name=build_su2         # Job name
#SBATCH -p genoa
#SBATCH --time=00:10:00              # Maximum runtime (adjust as needed)
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks=1                   # Total tasks
#SBATCH --cpus-per-task=16            # CPU cores per task
#SBATCH --output=su2_build_%j.log    # Output log file

# Load required modules
module load 2022
module load foss/2022a

# Clone SU2 repository if not done already


/home/tberk/SU2/meson.py build -Dwarning_level=2 -Denable-autodiff=true -Denable-directdiff=false -Dwith-mpi=enabled -Denable-cgns=true -Denable-tecio=false -Denable-mlpcpp=true --prefix=/gpfs/home5/tberk/SU2

/home/tberk/SU2/ninja -C build install -j16



