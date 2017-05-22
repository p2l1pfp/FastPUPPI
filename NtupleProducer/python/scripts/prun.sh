CODE=${1/.py/}; shift
MAIN=/eos/cms/store/cmst3/user/gpetrucc/l1phase2/Spring17D/200517
PREFIX="inputs_17D"

INPUT=$1; shift
if [[ "$1" != "" ]]; then
    OUTPUT="$1";
    shift;
else
    OUTPUT="${INPUT}";
fi;

if [[ "$1" == "--01" ]]; then
    MAIN=/eos/cms/store/cmst3/user/gpetrucc/l1phase2/Spring17D/010517
    shift;
elif [[ "$1" == "--phil" ]]; then
    MAIN=/eos/cms/store/cmst3/user/pharris/L1PF/inputs/$INPUT/
    PREFIX=""
    shift;
fi;

#~gpetrucc/pl/cmsSplit.pl --files "$MAIN/${INPUT}/inputs_17D*root" --label ${OUTPUT} ${CODE}.py --bash --n 8 $* && bash ${CODE}_${OUTPUT}_local.sh && bash ${CODE}_${OUTPUT}_cleanup.sh
~gpetrucc/pl/cmsSplit.pl --files "$MAIN/${PREFIX}*${INPUT}*root" --label ${OUTPUT} ${CODE}.py --bash --n 8 $* && bash ${CODE}_${OUTPUT}_local.sh && bash ${CODE}_${OUTPUT}_cleanup.sh
