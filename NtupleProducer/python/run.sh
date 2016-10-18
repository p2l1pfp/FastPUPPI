#!/bin/bash

eosdir=$1
file=$2
dir=$3
label=$4 

cp $dir/runNtupleProducer_cfg_tmp.py .
sed "s@XXXX@$eosdir/$file@g" runNtupleProducer_cfg_tmp.py > runNtupleProducer_cfg.py
cd $dir
eval `scramv1 runtime -sh`
cd -
cmsRun runNtupleProducer_cfg.py
#mv ntuple.root  $dir/F${label}$file
mv MetFile.root $dir/M${label}$file
#mv puppi.root $dir/PUP${label}$file
rm puppi.root