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

    ##if [   -d ${THIS_JOBSUBDIR} ]; then echo "${THIS_JOBSUBDIR} exists. Removing..."; rm -rf ${THIS_JOBSUBDIR}; fi
    ##if [   -d ${THIS_SIMDIR} ];    then echo "${THIS_SIMDIR} exists. Removing...";    rm -rf ${THIS_SIMDIR};    fi

    if [ ! -d ${THIS_JOBSUBDIR} ]; then mkdir -p ${THIS_JOBSUBDIR}; fi
    if [ ! -d ${THIS_SIMDIR} ];    then mkdir -p ${THIS_SIMDIR};    fi
    

    cd ${THIS_JOBSUBDIR}
    #mkdir snptest bimbam blore caviarbf finemap paintor metaskat
    rm -rf metaskat
    mkdir metaskat

    #source ${SUBDIR}/create_phenotype_with_features		# PHENO_JOBNAME
    #cd snptest
    #source ${SUBDIR}/snptest					# SNPTEST_JOBNAME
    #cd ../blore
    #source ${SUBDIR}/blore_summary				# BLORE_SUMMARY_JOBNAME
    #cd ../bimbam
    #source ${SUBDIR}/bimbam_summary				# BIMBAM_SUMMARY_JOBNAME
    #cd ../snptest
    #source ${SUBDIR}/meta					# META_JOBNAME
    #cd ../blore
    #source ${SUBDIR}/blore_meta_without_features		# BLORE_META_JOBNAME
    #source ${SUBDIR}/blore_meta_with_features			# BLORE_META_FEAT_JOBNAME
    #cd ../bimbam
    #source ${SUBDIR}/bimbam_meta				# BIMBAM_META_JOBNAME
    #cd ..
    #source ${SUBDIR}/weighted_ld				# WGT_LD_JOBNAME
    #cd paintor
    #source ${SUBDIR}/paintor					# PAINTOR_JOBNAME
    #source ${SUBDIR}/paintor_fa  				# PAINTORFA_JOBNAME
    #cd ../caviarbf
    #source ${SUBDIR}/caviarbf					# CAVIARBF_JOBNAME
    #cd ../finemap
    #source ${SUBDIR}/finemap
    #cd ../metaskat
    cd metaskat
    source ${SUBDIR}/metaskat

    cd $CURDIR

done
