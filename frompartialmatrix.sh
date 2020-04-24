#!/usr/bin/env bash
#
#SBATCH --account sdumathjea_slim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time 1:00:00            # max time (HH:MM:SS)
#SBATCH --mail-type=END,FAIL
#
# Author: Yuki Koyanagi
# History:
#


echo Running on "$(hostname)"
echo Running job: "$SLURM_JOB_NAME"
echo Available nodes: "$SLURM_NODELIST"
echo Number of nodes: "$SLURM_NNODES"
echo Slurm_submit_dir: "$SLURM_SUBMIT_DIR"
echo Start time: "$(date)"
start=$(date +%s)

echo Enable modules
module purge
module add python-intel/3.6.3

echo Starting Python program
DATADIR=/work/austmathjea/metastr
WORKDIR=/work/austmathjea/metastr/set9
OUTPUTDIR=${WORKDIR}/output_v2
OUTDIR=${WORKDIR}/motifs_v2
mkdir ${OUTDIR}

parallel --lb python frompartialmatrix.py ${OUTPUTDIR}/{}.mat.pkl \
	 -s ${OUTDIR} :::: ${WORKDIR}/validation_v2.lst \
	 > ${WORKDIR}/fpm.v2.log


#cd ${OUTDIR}
#tar cjf ${WORKDIR}/pkls_${SUFFIX}.tar.bz2 ./*

end=$(date +%s)
echo End time: "$(date)"
dur=$(date -d "0 $end sec - $start sec" +%T)
echo Duration: "$dur"

echo Done.
