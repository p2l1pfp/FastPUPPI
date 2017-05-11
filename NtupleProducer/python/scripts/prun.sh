CODE=${1/.py/}; shift
INPUT=$1; shift
if [[ "$1" != "" ]]; then
    OUTPUT="$1";
    shift;
else
    OUTPUT="${INPUT}";
fi;
MAIN=/eos/cms/store/cmst3/user/gpetrucc/l1phase2/Spring17D/010517
#~gpetrucc/pl/cmsSplit.pl --files "$MAIN/${INPUT}/inputs_17D*root" --label ${OUTPUT} ${CODE}.py --bash --n 8 $* && bash ${CODE}_${OUTPUT}_local.sh && bash ${CODE}_${OUTPUT}_cleanup.sh
~gpetrucc/pl/cmsSplit.pl --files "$MAIN/inputs_17D*${INPUT}*root" --label ${OUTPUT} ${CODE}.py --bash --n 8 $* && bash ${CODE}_${OUTPUT}_local.sh && bash ${CODE}_${OUTPUT}_cleanup.sh
