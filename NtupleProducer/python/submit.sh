#!/bin/bash

#for x in `eos ls eos/cms/store/mc/TTI2023Upg14D/SinglePionPlusFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/`; do
#    bsub -q 8nh run.sh /store/mc/TTI2023Upg14D/SinglePionPlusFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/ $x $PWD
#done

#for x in `eos ls eos/cms/store/mc/TTI2023Upg14D/SinglePionMinusFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/`; do
#    bsub -q 8nh run.sh /store/mc/TTI2023Upg14D/SinglePionMinusFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/ $x $PWD
#done

#for  x in `eos ls  eos/cms/store/mc/TTI2023Upg14D/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/`; do
#     bsub -q 8nh run.sh /store/mc/TTI2023Upg14D/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/ $x $PWD
#done

#for  x in `eos ls  eos/cms/store/mc/TTI2023Upg14D/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/`; do
#     bsub -q 8nh run.sh /store/mc/TTI2023Upg14D/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/ $x $PWD
#done
#cfg=$1
#dir=/store/group/cmst3/user/pharris/L1PF/Pion/
#for  x in `eos ls  $dir `; do
#    bsub -q 8nh -o out.%J run.sh $dir $x $PWD Pion $cfg
#done
#dir=/store/cmst3/user/gpetrucc/l1phase2/081116/DR_PU0/2023D3_guns/
#for  x in `eos ls  $dir `; do
#    bsub -q 8nh -o out.%J run.sh $dir $x $PWD Photon $cfg
#done
#dir=/store/relval/CMSSW_8_1_0_pre15/RelValSingleMuPt100Extended/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D4PU140-v1
#cfg=runNtupleProducer_cfg_tmp_pt0.py
#for  x in `eos ls  $dir `; do
#    bsub -q 1nh -o out.%J run.sh $dir $x $PWD Mu $cfg
#done
#exit
cfg=runNtupleProducer_cfg_tmp_pt0.py
dir=/store/relval/CMSSW_8_1_0_pre15/RelValZMM_13/GEN-SIM-RECO/81X_upgrade2023_realistic_v3_2023D3Timing13TeV-v1/10000/
for  x in `eos ls  $dir `; do
    bsub -q 1nh -o out.%J run.sh $dir $x $PWD G0PU $cfg
done
cfg=runNtupleProducer_cfg_tmp_pt4.py
dir=/store/relval/CMSSW_8_1_0_pre15/RelValZMM_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/
for  x in `eos ls  $dir `; do
    bsub -q 1nh -o out.%J run.sh $dir $x $PWD G200PU $cfg
done
dir=/store/relval/CMSSW_8_1_0_pre15/RelValZMM_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/
for  x in `eos ls  $dir `; do
    bsub -q 1nh -o out.%J run.sh $dir $x $PWD G140PU $cfg
done
cfg=runNtupleProducer_cfg_tmp_pt0.py
dir=/store/relval/CMSSW_8_1_0_pre15/RelValTTbar_13/GEN-SIM-RECO/81X_upgrade2023_realistic_v3_2023D3Timing13TeV-v1/10000/
for  x in `eos ls  $dir`; do
     bsub -q 8nh -o out.%J run.sh $dir $x $PWD Top0PU $cfg
done
cfg=runNtupleProducer_cfg_tmp_pt4.py
dir=/store/relval/CMSSW_8_1_0_pre15/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/
for  x in `eos ls  $dir`; do
     bsub -q 8nh -o out.%J run.sh $dir $x $PWD Top140PU $cfg
done
dir=/store/relval/CMSSW_8_1_0_pre15/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/
for  x in `eos ls  $dir`; do
     bsub -q 8nh -o out.%J run.sh $dir $x $PWD Top200PU $cfg
done
