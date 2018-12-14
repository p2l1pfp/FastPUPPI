CODE=${1/.py/}; shift
MAIN=/eos/cms/store/cmst3/group/l1tr/gpetrucc/101X/NewInputs/231018/$1
PREFIX="inputs_"

if [[ "$1" == "--3D" ]]; then
    shift;
    MAIN=/eos/cms/store/cmst3/group/l1tr/gpetrucc/101X/NewInputs/011118/$1
    PREFIX="reinputs_"
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
$PSCRIPT/cmsSplit.pl --files "$MAIN/${PREFIX}*${INPUT}*root" --label ${OUTPUT} ${CODE}.py --bash --n 8 $* && bash ${CODE}_${OUTPUT}_local.sh 

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
