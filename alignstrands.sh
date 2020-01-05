#!/usr/bin/env bash
#
#SBATCH --account austmathjea_slim
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --time 4:00:00            # max time (HH:MM:SS)
#SBATCH --mail-type=END,FAIL
#
# Author: Yuki Koyanagi
# History:
#

SUFFIX=v2
FACTOR=4
GAP=-4
BIAS="[1, .2, .2, .1, .1]"
SUBMAT=blosum62plus
LS=Strict
OPER=*
CUTOFF="[0.1,]"

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

sleep 30

if [ "$OPER" == "*" ]; then
    FACTOR=1
fi
echo Mofify python program
sed -i -E "s/FACTOR = [[:digit:]]+/FACTOR = ${FACTOR}/" alignstrands.py
sed -i -E "s/GAP_SCORE = -?[[:digit:]]+/GAP_SCORE = ${GAP}/" alignstrands.py
#sed -i -E "s/BIAS = \[.*\]/BIAS = ${BIAS}/" alignstrands.py
sed -i -E "s/mat = smat\.[[:alnum:]]+/mat = smat\.${SUBMAT}/" alignstrands.py
sed -i -E "s/(Local|Strict)/${LS}/g" alignstrands.py
sed -i -E "s/score (\+|\*)=/score ${OPER}=/" alignstrands.py
sed -i -E "s/^CUTOFF_VALUES = [][[:digit:].,]+/CUTOFF_VALUES = ${CUTOFF}/" alignstrands.py 

echo Starting Python program
DATADIR=/work/austmathjea/metastr
WORKDIR=/work/austmathjea/metastr/set3
OUTDIR=${WORKDIR}/output_${SUFFIX}
mkdir ${OUTDIR}
srun python alignstrands.py ${WORKDIR}/test.lst \
     ${DATADIR}/summaries2017 ${DATADIR}/allprot_17 \
     ${WORKDIR}/gammasegs.txt -o ${OUTDIR} -m -s

echo Done pairing

cd ${OUTDIR}
tar cjf ${WORKDIR}/pkls_${SUFFIX}.tar.bz2 ./*

end=$(date +%s)
echo End time: "$(date)"
dur=$(date -d "0 $end sec - $start sec" +%T)
echo Duration: "$dur"

echo Done.
