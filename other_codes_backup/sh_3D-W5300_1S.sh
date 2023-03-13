#!/bin/bash
#SBATCH --time=02-00:00       # Time limit dd-hh:mm
#SBATCH --nodes=1             # Number of compute nodes
#SBATCH --cpus-per-task=64    # Number of cores per node
#SBATCH --ntasks-per-node=1   # Do not change
#SBATCH --mem=0

#SBATCH --mail-user=zdeng57@uwo.ca
#SBATCH --mail-type=ALL

source ~/tensorflow/bin/activate
python ./server.py &

module load StdEnv/2020 ansys/2022R1

slurm_hl2hl.py --format ANSYS-FLUENT > machinefile
NCORE=$((SLURM_NTASKS * SLURM_CPUS_PER_TASK))

fluent 3ddp -t $NCORE -cnf=machinefile -mpi=intel -affinity=0 -g -i jou-fluent-start.jou
