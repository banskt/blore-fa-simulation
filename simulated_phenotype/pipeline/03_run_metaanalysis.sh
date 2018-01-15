#!/bin/bash

source CONFIG

SUBDIR="submission_scripts"

## Do we need to create weighted LD for this simulation?
source ${SUBDIR}/check_weighted_ld # creates RUN_WEIGHTED_LD variable and LOCIPREFIX array

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

    source ${SUBDIR}/create_phenotype 				# PHENO_JOBNAME
    source ${SUBDIR}/snptest					# SNPTEST_JOBNAME
    ## source ${SUBDIR}/gen_inflation # now included in meta
    source ${SUBDIR}/blore_summary				# BLORE_SUMMARY_JOBNAME
    source ${SUBDIR}/bimbam_summary				# BIMBAM_SUMMARY_JOBNAME
    source ${SUBDIR}/meta					# META_JOBNAME
    source ${SUBDIR}/blore_meta_without_features		# BLORE_META_JOBNAME
    source ${SUBDIR}/blore_meta_with_features			# BLORE_META_FEAT_JOBNAME
    source ${SUBDIR}/bimbam_meta				# BIMBAM_META_JOBNAME
    source ${SUBDIR}/weighted_ld				# WGT_LD_JOBNAME
    #source ${SUBDIR}/paintor
    #source ${SUBDIR}/paintorfa
    #source ${SUBDIR}/caviarbf
    #source ${SUBDIR}/finemap
    #source ${SUBDIR}/metacca

    cd $CURDIR

done

    ## # Run CAVIARBF
    ## CAVIARBF_C1_JOBNAME="caviarbf_c1_${SIM}_${RANDSTRING}"
    ## sed "s|_JOBNAME|${CAVIARBF_C1_JOBNAME}|g;
    ##     9s|_SIMDIR_|${SIMFOLDER}|;
    ##    10s|_NCAUSAL|1|"         ${CAVIARBF_MASTER_JOBSUB} > ${CAVIARBF_C1_JOBNAME}.bsub
    ## bsub -w "done(${META_JOBNAME})" < ${CAVIARBF_C1_JOBNAME}.bsub

    ## CAVIARBF_C2_JOBNAME="caviarbf_c2_${SIM}_${RANDSTRING}"
    ## sed "s|_JOBNAME|${CAVIARBF_C2_JOBNAME}|g;
    ##     9s|_SIMDIR_|${SIMFOLDER}|;
    ##    10s|_NCAUSAL|2|"         ${CAVIARBF_MASTER_JOBSUB} > ${CAVIARBF_C2_JOBNAME}.bsub
    ## bsub -w "done(${META_JOBNAME})" < ${CAVIARBF_C2_JOBNAME}.bsub

    ## # Run PAINTOR =====================================================================
    ## PAINTOR_JOBNAME="paintor_combinedLD_${SIM}_${RANDSTRING}"
    ## sed "s|_JOBNAME|${PAINTOR_JOBNAME}|g;
    ##     9s|_SIMDIR_|${SIMFOLDER}|"   ${PAINTOR_MASTER_JOBSUB} > ${PAINTOR_JOBNAME}.bsub
    ## bsub -w "done(${META_JOBNAME})" < ${PAINTOR_JOBNAME}.bsub

    ## # Run metaCCA =====================================================================
    ## METACCA_JOBNAME="metaCCA_${SIM}_${RANDSTRING}"
    ## sed "s|_JOBNAME|${METACCA_JOBNAME}|g;
    ##     9s|_SIMDIR_|${SIMFOLDER}|"   ${METACCA_MASTER_JOBSUB} > ${METACCA_JOBNAME}.bsub
    ## bsub -w "done(${META_JOBNAME})" < ${METACCA_JOBNAME}.bsub

    ## # Run FINEMAP =====================================================================
    ## FINEMAP_JOBNAME="finemap_c${NCAUSAL}_combined_${SIM}_${RANDSTRING}"
    ## sed "s|_JOBNAME|${FINEMAP_JOBNAME}|g;
    ##     9s|_SIMDIR_|${SIMFOLDER}|;
    ##    10s|_NCAUSAL|${NCAUSAL}|"     ${FINEMAP_MASTER_JOBSUB} > ${FINEMAP_JOBNAME}.bsub
    ## bsub -w "done(${META_JOBNAME})" < ${FINEMAP_JOBNAME}.bsub

    ## # Run PAINTOR over all annotations independently ==================================
    ## PAINTORFA_SEARCH_JOBNAME="paintor_combinedLD_FA_${SIM}_${RANDSTRING}"
    ## sed "s|_JOBNAME|${PAINTORFA_SEARCH_JOBNAME}|g;
    ##     9s|_SIMDIR_|${SIMFOLDER}|" \
    ##             ${PAINTORFA_SEARCHALL_MASTER_JOBSUB} > ${PAINTORFA_SEARCH_JOBNAME}.bsub
    ## bsub -w "done(${PAINTOR_JOBNAME})" < ${PAINTORFA_SEARCH_JOBNAME}.bsub

    ## # Run PAINTOR with selected annotations ===========================================
    ## PAINTORFA_SELECT_JOBNAME="paintor_combinedLD_FA_selected_${SIM}_${RANDSTRING}"
    ## sed "s|_JOBNAME|${PAINTORFA_SELECT_JOBNAME}|g;
    ##     9s|_SIMDIR_|${SIMFOLDER}|" \
    ##              ${PAINTORFA_COMBINED_MASTER_JOBSUB} > ${PAINTORFA_SELECT_JOBNAME}.bsub
    ## bsub -w "done(${PAINTORFA_SEARCH_JOBNAME})"  < ${PAINTORFA_SELECT_JOBNAME}.bsub
