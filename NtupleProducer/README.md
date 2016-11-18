cmsrel CMSSW_8_1_0_pre16   #  or SLHC12_patch1, with the merge-topic below, this is equivalent
cd CMSSW_8_1_0_pre16/src
cmsenv

#git cms-addpkg SLHCUpgradeSimulations/L1TrackTrigger
#git cms-merge-topic EmanuelPerez:TTI_62X_TrackTriggerObjects
#scramv1 b -j 8
