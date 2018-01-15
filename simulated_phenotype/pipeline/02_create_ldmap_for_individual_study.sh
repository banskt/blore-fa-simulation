#!/bin/bash

source CONFIG

JOBSUBDIR="${LDMAPDIR}/jobsub"

if [ ! -d ${JOBSUBDIR} ]; then mkdir -p ${JOBSUBDIR}; fi

cd ${JOBSUBDIR}
for STUDY in ${STUDYNAMES[@]}; do

    RESULTDIR="${LDMAPDIR}/${STUDY}"
    if [ -d ${RESULTDIR} ]; then rm -rf ${RESULTDIR}; fi
    if [ ! -d ${RESULTDIR} ]; then mkdir -p ${RESULTDIR}; fi

    JOBNAME="LD_individual_${STUDY}"
    sed -e "s|_JOBNAME_|${JOBNAME}|g;
            s|_DATADIR_|${DOSAGEDIR}/${STUDY}|g;
            s|_OUTDIR_|${RESULTDIR}|g;
            s|_LOCUSNAMES_|${LOCUSNAMES}|g;
            s|_LDSTORE_|${LDSTORE}|g" ${MASTER_JOBSUBDIR}/ldmap_individual_study.bsub > ${JOBSUBDIR}/${JOBNAME}.bsub

    bsub < ${JOBSUBDIR}/${JOBNAME}.bsub

done
cd ${CURDIR}
