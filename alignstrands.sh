#!/usr/bin/env bash
#
#SBATCH --account austmathjea_slim
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --time 6:00:00            # max time (HH:MM:SS)
#SBATCH --mail-type=END,FAIL
#
# Author: Yuki Koyanagi
# History:
#

FACTOR=1

echo Running on "$(hostname)"
echo Running job: "$SLURM_JOB_NAME"
echo Available nodes: "$SLURM_NODELIST"
echo Number of nodes: "$SLURM_NNODES"
echo Slurm_submit_dir: "$SLURM_SUBMIT_DIR"
echo Start time: "$(date)"
start=$(date +%s)

#echo Enable modules
#module purge
#module add python-intel/3.6.3

echo Mofify python program
sed -i -E "s/FACTOR = [[:digit:]]+/FACTOR = ${FACTOR}/" alignstrands.py

echo Starting Python program
DATADIR=/work/austmathjea/metastr
WORKDIR=/work/austmathjea/metastr/set1
OUTDIR=${WORKDIR}/output${FACTOR}
mkdir ${OUTDIR}
srun python alignstrands.py ${WORKDIR}/test.lst \
     ${DATADIR}/summaries2017 ${DATADIR}/allprot_17 \
     ${WORKDIR}/betapairs.txt -o ${OUTDIR} -m -s

echo Done pairing.

tar cjf ${WORKDIR}/pkls${FACTOR}.tar.bz2 ${OUTDIR}/*
rm -f ${OUTDIR}/*

end=$(date +%s)
echo End time: "$(date)"
dur=$(date -d "0 $end sec - $start sec" +%T)
echo Duration: "$dur"

echo Done.
