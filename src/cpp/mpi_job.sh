#!/bin/sh
#SBATCH -J GRNGAMPI
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --ntasks-per-node=20
#SBATCH --ntasks-per-core=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dennisgwilson@gmail.com
#SBATCH --time=100:00:00

WORK_DIR=/tmpdir/$LOGNAME/dennis/$SLURM_JOBID
mkdir $WORK_DIR
cp main $WORK_DIR
cp ga_log.conf $WORK_DIR
#cp -r mnist $WORK_DIR
#cp -r grns $WORK_DIR
cd $WORK_DIR
mkdir grns
mkdir logs

module purge
module load gcc/5.3.0
module load intelmpi/5.1.2.150
module list
export OMP_NUM_THREADS=1

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
srun ./main
