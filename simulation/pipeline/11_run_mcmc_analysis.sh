#!/bin/bash

source QL_MCMC_CONFIG

SUBDIR="${CURDIR}/submission_scripts"
#CHAIN_JOBS="True"
CHAIN_JOBS="False"

IFS=$'\r\n' GLOBIGNORE='*' command eval 'LOCIPREFIX=($(cat ${LOCUSNAMES}))'

for (( SIM=$START; SIM<=$END; SIM++ )); do 

    INDEX=`echo $SIM | awk '{printf "%03d", $1}'`
    SIMFOLDER="sim${INDEX}"
    RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 10 | head -n 1`

    THIS_JOBSUBDIR="${JOBSUBDIR}/${SIMFOLDER}"
    THIS_SIMDIR="${SIMDIR}/${SIMFOLDER}"

    #if [   -d ${THIS_JOBSUBDIR} ]; then rm -rf ${THIS_JOBSUBDIR}; fi
    #if [   -d ${THIS_SIMDIR} ];    then rm -rf ${THIS_SIMDIR};    fi

    if [ ! -d ${THIS_JOBSUBDIR} ]; then mkdir -p ${THIS_JOBSUBDIR}; fi
    if [ ! -d ${THIS_SIMDIR} ];    then mkdir -p ${THIS_SIMDIR};    fi
    

    cd ${THIS_JOBSUBDIR}

    #source ${SUBDIR}/create_phenotype 				# PHENO_JOBNAME
    #source ${SUBDIR}/snptest					# SNPTEST_JOBNAME
    ### source ${SUBDIR}/gen_inflation # now included in meta
    #source ${SUBDIR}/blore_summary				# BLORE_SUMMARY_JOBNAME
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
    #source ${SUBDIR}/metacca
    source ${SUBDIR}/gemma

    cd $CURDIR

done
