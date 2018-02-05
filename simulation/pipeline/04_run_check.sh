#!/bin/bash

S1="discarding /usr/users/sbanerj/miniconda3/envs/py35/bin from PATH"
S2="prepending /usr/users/sbanerj/miniconda3/envs/py35/bin to PATH"
S3="This problem is unconstrained."
S4=" Bad direction in the line search;"
S5="   refresh the lbfgs memory and restart the iteration."
S6=" Warning:  more than 10 function and gradient"
S7="   evaluations in the last line search.  Termination"
S8="   may possibly be caused by a bad search direction."
S9=" Line search cannot locate an adequate point after 20 function"
S10=" and gradient evaluations.  Previous x, f and g restored."
S11="Possible causes: 1 error in function or gradient evaluation;"
S12="                 2 rounding error dominate computation."
S13="^$"
S14="discarding /usr/users/sbanerj/miniconda3/bin from PATH"


source CONFIG

for (( SIM=$START; SIM<=$END; SIM++ )); do 

    INDEX=`echo $SIM | awk '{printf "%03d", $1}'`
    SIMFOLDER="sim${INDEX}"

    THIS_JOBSUBDIR="${JOBSUBDIR}/${SIMFOLDER}"
    THIS_SIMDIR="${SIMDIR}/${SIMFOLDER}"

    cd ${THIS_JOBSUBDIR}
    echo "Error files in ${SIMFOLDER}..."

    cat *.err | grep -v -e "${S1}" -e "${S2}"

    #for dir in snptest bimbam caviarbf finemap paintor; do
    for dir in blore; do
        cd ${dir};
        for file in `ls blore_meta_*_4_mu0.err`; do
        #for file in `ls *.err`; do
            cat ${file} | grep -v -e "${S1}" -e "${S2}" -e "${S3}" -e "${S4}" -e "${S5}" -e "${S6}" -e "${S7}" -e "${S8}" -e "${S9}" -e "${S10}" -e "${S11}" -e "${S12}" -e "${S13}" -e "${S14}"
        done
        cd ..
    done

    cd $CURDIR

done
