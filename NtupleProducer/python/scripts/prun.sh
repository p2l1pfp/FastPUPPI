CODE=${1/.py/}; shift
#MAIN=/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/$1/
MAIN=/eos/cms/store/cmst3/user/gpetrucc/l1phase2/93X/Inputs/$1/
PREFIX="inputs_17D"

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
fi

#~gpetrucc/pl/cmsSplit.pl --files "$MAIN/${INPUT}/inputs_17D*root" --label ${OUTPUT} ${CODE}.py --bash --n 8 $* && bash ${CODE}_${OUTPUT}_local.sh && bash ${CODE}_${OUTPUT}_cleanup.sh
~gpetrucc/pl/cmsSplit.pl --files "$MAIN/${PREFIX}*${INPUT}*root" --label ${OUTPUT} ${CODE}.py --bash --n 8 $* && bash ${CODE}_${OUTPUT}_local.sh && $clean && bash ${CODE}_${OUTPUT}_cleanup.sh
