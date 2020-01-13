#!/usr/bin/env bash
#
#SBATCH --account austmathjea_slim
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1
#SBATCH --time 4:00:00            # max time (HH:MM:SS)
#SBATCH --mail-type=END,FAIL
#
# Author: Yuki Koyanagi
# History:
#

SUFFIX=v2
GAP=-4
SUBMAT=blosum62plus
LS=Strict

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
sed -i -E "s/GAP_SCORE = -?[[:digit:]]+/GAP_SCORE = ${GAP}/" makepartialmatrix.py
sed -i -E "s/mat = smat\.[[:alnum:]]+/mat = smat\.${SUBMAT}/" makepartialmatrix.py
sed -i -E "s/(Local|Strict)/${LS}/g" makepartialmatrix.py

echo Starting Python program
DATADIR=/work/austmathjea/metastr
WORKDIR=/work/austmathjea/metastr/set3
OUTDIR=${WORKDIR}/output_${SUFFIX}
mkdir ${OUTDIR}
srun python makepartialmatrix.py ${WORKDIR}/test.lst \
     ${DATADIR}/summaries2017 ${WORKDIR}/gammasegs.txt \
     -o ${OUTDIR} -s -u5

echo Done pairing

cd ${OUTDIR}
tar cjf ${WORKDIR}/pkls_${SUFFIX}.tar.bz2 ./*

end=$(date +%s)
echo End time: "$(date)"
dur=$(date -d "0 $end sec - $start sec" +%T)
echo Duration: "$dur"

echo Done.
