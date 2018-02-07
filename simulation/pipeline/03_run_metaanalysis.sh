#!/bin/bash

source CONFIG

SUBDIR="${CURDIR}/submission_scripts"
#CHAIN_JOBS="True"
CHAIN_JOBS="False"

## Do we need to create weighted LD for this simulation?
source ${SUBDIR}/check_weighted_ld # creates RUN_WEIGHTED_LD variable and LOCIPREFIX array

for (( SIM=$START; SIM<=$END; SIM++ )); do 

    INDEX=`echo $SIM | awk '{printf "%03d", $1}'`
    SIMFOLDER="sim${INDEX}"
    RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 10 | head -n 1`

    THIS_JOBSUBDIR="${JOBSUBDIR}/${SIMFOLDER}"
    THIS_SIMDIR="${SIMDIR}/${SIMFOLDER}"

    if [   -d ${THIS_JOBSUBDIR} ]; then echo "${THIS_JOBSUBDIR} exists. Removing..."; rm -rf ${THIS_JOBSUBDIR}; fi
    if [   -d ${THIS_SIMDIR} ];    then echo "${THIS_SIMDIR} exists. Removing...";    rm -rf ${THIS_SIMDIR};    fi

    if [ ! -d ${THIS_JOBSUBDIR} ]; then mkdir -p ${THIS_JOBSUBDIR}; fi
    if [ ! -d ${THIS_SIMDIR} ];    then mkdir -p ${THIS_SIMDIR};    fi
    

    cd ${THIS_JOBSUBDIR}

    #source ${SUBDIR}/create_phenotype_with_features		# PHENO_JOBNAME
    source ${SUBDIR}/create_phenotype_random
    #source ${SUBDIR}/snptest					# SNPTEST_JOBNAME
    source ${SUBDIR}/blore_summary				# BLORE_SUMMARY_JOBNAME
    #source ${SUBDIR}/bimbam_summary				# BIMBAM_SUMMARY_JOBNAME
    #source ${SUBDIR}/meta					# META_JOBNAME
    #source ${SUBDIR}/blore_meta_without_features		# BLORE_META_JOBNAME
    #source ${SUBDIR}/blore_meta_with_features			# BLORE_META_FEAT_JOBNAME
    #source ${SUBDIR}/bimbam_meta				# BIMBAM_META_JOBNAME
    #source ${SUBDIR}/weighted_ld				# WGT_LD_JOBNAME
    #source ${SUBDIR}/paintor					# PAINTOR_JOBNAME
    #source ${SUBDIR}/paintor_fa  				# PAINTORFA_JOBNAME
    #source ${SUBDIR}/caviarbf					# CAVIARBF_JOBNAME
    #source ${SUBDIR}/finemap
    #source ${SUBDIR}/metaskat
    #source ${SUBDIR}/jam

    cd $CURDIR

done
