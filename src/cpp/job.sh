#!/bin/sh
#SBATCH -J GAGRNANN
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-core=1
#SBATCH --mail-user=dennisgwilson@gmail.com
#SBATCH --mail-type=ALL

BINARY=./main

WORK_DIR=$(mktemp -d -p /tmpdir/franklin/dennis/)

cp $BINARY $WORK_DIR
cp ga_log.conf $WORK_DIR
cp -r mnist $WORK_DIR
cp -r grns $WORK_DIR

echo "Working in $WORK_DIR"
cd $WORK_DIR
mkdir logs
./main
