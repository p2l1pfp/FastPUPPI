CODE=${1/.py/}; shift
#MAIN=/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/$1/
#MAIN=/eos/cms/store/cmst3/user/gpetrucc/l1phase2/93X/Inputs/$1/
MAIN=/eos/cms/store/cmst3/user/gpetrucc/l1phase2/101X/NewInputs/040418/$1
PREFIX="inputs_"

INPUT=$1; shift
if [[ "$1" != "" ]]; then
    OUTPUT="$1";
    shift;
else
    OUTPUT="${INPUT}";
fi;

if [[ "$1" == "--92X" ]]; then
    MAIN=/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/$INPUT/
    shift;
fi;
clean=true
if [[ "$1" == "--noclean" ]]; then
    clean=false;
    shift
fi

#~gpetrucc/pl/cmsSplit.pl --files "$MAIN/${INPUT}/inputs_17D*root" --label ${OUTPUT} ${CODE}.py --bash --n 8 $* && bash ${CODE}_${OUTPUT}_local.sh && bash ${CODE}_${OUTPUT}_cleanup.sh
~gpetrucc/pl/cmsSplit.pl --files "$MAIN/${PREFIX}*${INPUT}*root" --label ${OUTPUT} ${CODE}.py --bash --n 8 $* && bash ${CODE}_${OUTPUT}_local.sh 

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
