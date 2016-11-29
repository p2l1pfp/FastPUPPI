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
cfg=$1
dir=/store/group/cmst3/user/pharris/L1PF/Pion/
for  x in `eos ls  $dir `; do
    bsub -q 8nh -o out.%J run.sh $dir $x $PWD Pion $cfg
done

dir=/store/cmst3/user/gpetrucc/l1phase2/081116/DR_PU0/2023D3_guns/
for  x in `eos ls  $dir `; do
    bsub -q 8nh -o out.%J run.sh $dir $x $PWD Photon $cfg
done

exit
dir=/store/relval/CMSSW_8_1_0_pre15/RelValZMM_13/GEN-SIM-RECO/81X_upgrade2023_realistic_v3_2023D3Timing13TeV-v1/10000/
for  x in `eos ls  $dir `; do
    bsub -q 8nh -o out.%J run.sh $dir $x $PWD G0PU $cfg
done
dir=/store/relval/CMSSW_8_1_0_pre15/RelValZMM_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/
for  x in `eos ls  $dir `; do
    bsub -q 8nh -o out.%J run.sh $dir $x $PWD G200PU $cfg
done
dir=/store/relval/CMSSW_8_1_0_pre15/RelValZMM_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/
#dir=/store/group/cmst3/user/pharris/L1PF/ZMM/
#dir=/store/group/cmst3/user/pharris/L1PF/Pion/
#dir=/store/cmst3/user/gpetrucc/l1phase2/081116/DR_PU0/2023D3_guns/
#for  x in `eos ls  $dir | grep Photon | grep FEVT `; do
for  x in `eos ls  $dir `; do
    bsub -q 8nh -o out.%J run.sh $dir $x $PWD G140PU $cfg
done
exit
cfg=$1
#dir=/store/mc/TTI2023Upg14D/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/
dir=
dir=/store/relval/CMSSW_8_1_0_pre15/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/
dir=/store/relval/CMSSW_8_1_0_pre15/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/
for  x in `eos ls  $dir`; do
     bsub -q 8nh -o out.%J run.sh $dir $x $PWD Top $cfg
done

dir=/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/
for  x in `eos ls  $dir`; do
    bsub -q 8nh -o out.%J run.sh $dir $x $PWD Zmm $cfg
done
exit
#dir=/store/mc/TTI2023Upg14D/SinglePionMinusFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/
#for  x in `eos ls  $dir`; do
#     bsub -q 8nh -o out.%J run.sh $dir $x $PWD PiMinus
#done

#dir=/store/mc/TTI2023Upg14D/SinglePionPlusFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/
#for  x in `eos ls  $dir`; do
#     bsub -q 8nh run.sh $dir $x $PWD PiPlus
#done