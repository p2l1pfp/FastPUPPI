Basic Instructions

```
cmsrel CMSSW_9_1_0_pre2
cd CMSSW_9_1_0_pre2/src
cmsenv
git cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline
git cms-merge-topic -u cms-l1t-offline:phase2-l1t-integration-CMSSW_9_1_0_pre2
git-cms-addpkg L1Trigger/L1THGCal
git clone https://github.com/cms-data/L1Trigger-L1THGCal.git L1Trigger/L1THGCal/data
git cms-merge-topic gpetruc:L1EGC_rebased-CMSSW_9_1_0_pre2
git cms-addpkg DataFormats/Phase2L1CaloTrig
git cms-addpkg L1Trigger/L1CaloTrigger
git clone git@github.com:nhanvtran/FastPUPPI.git -b 91X_newClusterer
scram b -j8
```

The first step is to produce the inputs:
`FastPUPPI/NtupleProducer/python/runInputs.py`

The trigger MC can be found on DAS `dataset=/*/*PhaseIISpring17D*/*`
