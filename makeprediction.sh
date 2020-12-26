#!/usr/bin/env bash
#
# Usage: makeprediction.sh {pred_no} {log_no} {lst_file} {top} {#jobs}

WORKDIR=../data/set13
BP_DIR=${WORKDIR}/bp_mat
OUTDIR=${WORKDIR}/prediction${1}
LOGDIR=${WORKDIR}/log
LOG=${LOGDIR}/mp.${1}.${2}.log

mkdir -p ${OUTDIR}
mkdir -p ${LOGDIR}

echo "Start time: $(date)" > ${LOG}
start=$(date +%s)

echo "Preparing environment" >> ${LOG}
sudo yum update -y
sudo yum install -y parallel
pip3 install --user numpy
pip3 install --user permutation
pip3 install --user ../../packages/fatgraph-1.0.2.tar.gz
echo "Done. Time: $(date)" >> ${LOG}

echo "Starting Python program" >> ${LOG}
parallel -j ${5} python3 makeprediction.py ${BP_DIR}/{}.json \
	 -s ${OUTDIR}/{}.json -a0.1 -t${4} :::: ${WORKDIR}/${3}

end=$(date +%s)
echo "End time: $(date)" >> ${LOG}
dur=$(date -d "0 $end sec - $start sec" +%T)
echo "Duration: $dur" >> ${LOG}

exit
