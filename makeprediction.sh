#!/usr/bin/env bash
#
#SBATCH --account sdumathjea_slim
#SBATCH --nodes 1                 # number of nodes
#SBATCH --time 3:00:00            # max time (HH:MM:SS)
#SBATCH --mail-type=FAIL,END
#

WORKDIR=/work/austmathjea/metastr/set13
BP_DIR=${WORKDIR}/bp_mat
OUTDIR=${WORKDIR}/prediction1

echo Running on "$(hostname)"
echo Running job: "$SLURM_JOB_NAME"
echo Available nodes: "$SLURM_NODELIST"
echo Slurm_submit_dir: "$SLURM_SUBMIT_DIR"
echo Start time: "$(date)"
start=$(date +%s)

# echo Clearing scratch folder
# rm -f ${SCRATCH}/*

#echo Enable modules
#module purge
#module add python-intel/3.5.2

mkdir ${OUTDIR}

echo Starting Python program
parallel python makeprediction.py ${BP_DIR}/{}.json \
	 -s ${OUTDIR}/{}.json :::: ${WORKDIR}/validation_8.lst

end=$(date +%s)
echo End time: "$(date)"
dur=$(date -d "0 $end sec - $start sec" +%T)
echo Duration: "$dur"

echo Done.
