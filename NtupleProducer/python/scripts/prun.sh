#!/bin/bash
CODE=${1/.py/}; shift

if [[ "$CODE" == "" || "$CODE" == "--help" || "$CODE" == "-h" || "$CODE" == "-?" ]]; then
    echo "scripts/prun.sh config.py [-j <N>] --<version>  <Sample> <OutputName> [ --noclean ] [ --inline-customize \"<options>\" ]"
    exit 0;
fi

N=8
if [[ "$1" == "-j" ]]; then N=$2; shift; shift; fi;

if [[ "$1" == "--106X_v0" ]]; then
    shift;
    MAIN=/eos/cms/store/cmst3/group/l1tr/gpetrucc/11_1_X/NewInputs106X/200420.done/$1
    PREFIX="inputs106X_"
elif [[ "$1" == "--110X_v0" ]]; then
    shift;
    MAIN=/eos/cms/store/cmst3/group/l1tr/gpetrucc/11_1_X/NewInputs110X/280420.done/$1
    PREFIX="inputs110X_"
    if echo $CODE | grep -q 110X; then
        echo "Assume $CODE is ready for 110X";
    else
        echo "Convert ${CODE}.py to ${CODE/_106X/}_110X.py automatically updating era and geometry";
        sed -e 's+Phase2C8+Phase2C9+g' -e 's+GeometryExtended2026D41+GeometryExtended2026D49+' ${CODE}.py > ${CODE/_106X/}_110X.py;
        CODE=${CODE/_106X/}_110X
    fi

else 
    echo "You mush specify the version of the input samples to run on "
    echo "   --106X_v0 : 106X L1T MC inputs, with default HGCal TPs "
    echo "   --110X_v0 : 110X HLT MC inputs, with default HGCal TPs "

fi;
 
if [[ "$L1TPF_LOCAL_INPUT_DIR" != "" ]] && test -d $L1TPF_LOCAL_INPUT_DIR; then
    L1TPF_LOCAL_MAIN=$L1TPF_LOCAL_INPUT_DIR/$(basename $(dirname $MAIN))/$(basename $MAIN)
    if test -d $L1TPF_LOCAL_MAIN; then
        echo "Will use local directory $L1TPF_LOCAL_MAIN";
        MAIN=$L1TPF_LOCAL_MAIN;
    fi;
fi;

INPUT=$1; shift
if [[ "$1" != "" ]]; then
    OUTPUT="$1";
    shift;
else
    OUTPUT="${INPUT}";
fi;

clean=true
if [[ "$1" == "--noclean" ]]; then
    clean=false;
    shift
fi

PSCRIPT=$CMSSW_BASE/src/FastPUPPI/NtupleProducer/python/scripts/
$PSCRIPT/cmsSplit.pl --files "$MAIN/${PREFIX}*root" --label ${OUTPUT} ${CODE}.py --bash --n $N $* && bash ${CODE}_${OUTPUT}_local.sh 

if $clean; then
    REPORT=$(grep -o "tee \S\+_${OUTPUT}.report.txt" ${CODE}_${OUTPUT}_local.sh  | awk '{print $2}');
    if [[ "$REPORT" != "" ]] && test -f ${REPORT}; then
        if grep -q ', 0 failed' $REPORT; then
           bash ${CODE}_${OUTPUT}_cleanup.sh
        else
           echo "Failed jobs... will not delete anything."
           exit 1;
        fi;
    else
        bash ${CODE}_${OUTPUT}_cleanup.sh
    fi;
fi
