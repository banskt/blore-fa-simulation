#!/bin/bash

START=1
END=50

CURDIR=`pwd`
BASEDIR="/scratch/sbanerj/quasi_laplace_gwas/simulated_phenotype/mcmc_vs_qL"
COPYFROM="/scratch/sbanerj/quasi_laplace_gwas/simulated_phenotype/meta_nold"
DOSAGEDIR="${BASEDIR}/loci_dosages"
STUDYNAMES=('G1' 'G2' 'G3' 'G4' 'G5')
STUDYSAMPLES=('2139' '2420' '2472' '2084' '3967')
LOCUSNAMES="${BASEDIR}/LOCUSNAMES"

SCRIPTDIR="${CURDIR}/../scripts"
MASTER_JOBSUBDIR="${CURDIR}/../jobsubfiles"
QCTOOL="/usr/users/sbanerj/packages/qctool/qctool_v1.4-linux-x86_64/qctool"

JOBSUBDIR="${CURDIR}/mcmc_vs_qL/create_loci"
REF_DOSAGEDIR="/scratch/sbanerj/quasi_laplace_gwas/simulated_phenotype/meta_nold/loci_dosages"

if [ ! -d ${BASEDIR} ]; then mkdir -p ${BASEDIR}; fi
cp ${COPYFROM}/LOCUSNAMES ${LOCUSNAMES}


IFS=$'\r\n' GLOBIGNORE='*' command eval 'LOCIPREFIX=($(cat ${LOCUSNAMES}))'

CLOCICOM="${SCRIPTDIR}/create_loci_common_SNPs.py"

if [ ! -d ${JOBSUBDIR} ]; then mkdir -p ${JOBSUBDIR}; fi

for STUDY in ${STUDYNAMES[@]}; do
    if [ ! -d ${DOSAGEDIR}/${STUDY} ]; then mkdir -p ${DOSAGEDIR}/${STUDY}; fi
    cp ${REF_DOSAGEDIR}/${STUDY}/*.sample ${DOSAGEDIR}/${STUDY}/
done

cd ${JOBSUBDIR}

for LOCUSPREFIX in ${LOCIPREFIX[@]}; do

    JOBNAME="common_SNPs_${LOCUSPREFIX}"
    sed -e "s|_JOBNAME|${JOBNAME}|g;
            s|_SCRIPT_|${CLOCICOM}|g;
            s|_QCTOOL_|${QCTOOL}|g;
            s|_LOCUSP_|${LOCUSPREFIX}|g;
            s|_STUDYN_|\"${STUDYNAMES[*]}\"|g;
            s|_LOCIDO_|${REF_DOSAGEDIR}|g;
            s|_LOCIDN_|${DOSAGEDIR}|g;" ${MASTER_JOBSUBDIR}/create_loci_common_SNPs.bsub > ${JOBNAME}.bsub
    bsub < ${JOBNAME}.bsub
done

cd ${CURDIR}
