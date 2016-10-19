#!/bin/bash

eosdir=$1
file=$2
dir=$3
label=$4 
cfg=$5

cp $dir/$cfg .
sed "s@XXXX@$eosdir/$file@g" $cfg > runNtupleProducer_cfg.py
cd $dir
eval `scramv1 runtime -sh`
cd -
cmsRun runNtupleProducer_cfg.py
#mv ntuple.root  $dir/F${label}$file
mv MetFile.root $dir/M${label}$file