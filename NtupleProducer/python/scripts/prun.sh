CODE=${1/.py/}; shift
MAIN=/eos/cms/store/cmst3/group/l1tr/gpetrucc/106X/NewInputs104X/240719.done/$1
PREFIX="inputs104X_"
N=8

if [[ "$1" == "-j" ]]; then N=$2; shift; shift; fi;

if [[ "$1" == "--v1" ]]; then # this is the default but we keep anyway
    shift;
    MAIN=/eos/cms/store/cmst3/group/l1tr/gpetrucc/106X/NewInputs104X/240719.done/$1
    PREFIX="inputs104X_"
elif [[ "$1" == "--v1_fat" ]]; then # this is the default but we keep anyway
    shift;
    MAIN=/eos/cms/store/cmst3/group/l1tr/gpetrucc/106X/NewInputs104X/240719_fat.done/$1
    PREFIX="inputs104X_"
elif [[ "$1" == "--v2_fat" ]]; then # this is the default but we keep anyway
    shift;
    MAIN=/eos/cms/store/cmst3/group/l1tr/gpetrucc/106X/NewInputs104X/240719_newhgc_try2_fat.done/$1
    PREFIX="inputs104X_"
elif [[ "$1" == "--v3_fat" ]]; then # this is the default but we keep anyway
    shift;
    MAIN=/eos/cms/store/cmst3/group/l1tr/gpetrucc/106X/NewInputs104X/300819_hgcv3_fat.done/$1
    PREFIX="inputs104X_"
elif [[ "$1" == "--v0" ]]; then # this is the default but we keep anyway
    shift;
    MAIN=/eos/cms/store/cmst3/group/l1tr/gpetrucc/106X/NewInputs104X/240719_oldhgc.done/$1
    PREFIX="inputs104X_"
elif [[ "$1" == "--106X_v0" ]]; then # this is the default but we keep anyway
    shift;
    MAIN=/eos/cms/store/cmst3/group/l1tr/gpetrucc/106X/NewInputs106X/250819.done/$1
    PREFIX="inputs106X_"
    if echo $CODE | grep -q 106X; then
        echo "Assume $CODE is ready for 106X";
    else
        echo "Convert ${CODE}.py to ${CODE/_104X/}_106X.py automatically updating era and geometry";
        sed -e 's+Phase2C4+Phase2C8+g' -e 's+GeometryExtended2023D35+GeometryExtended2023D41+' ${CODE}.py > ${CODE/_104X/}_106X.py;
        CODE=${CODE/_104X/}_106X
    fi
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
