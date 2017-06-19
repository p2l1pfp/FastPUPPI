Basic Instructions

```
cmsrel CMSSW_9_2_0
cd CMSSW_9_2_0/src
cmsenv
git cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline
git cms-merge-topic -u cms-l1t-offline:l1t-phase2-v1.12
git cms-merge-topic skinnari:Tracklet_92X
git clone git@github.com:p2l1pfp/FastPUPPI.git -b 91X_newClusterer
scram b -j8
```

The first step is to produce the inputs:
`FastPUPPI/NtupleProducer/python/runInputs.py`

The trigger MC can be found on DAS `dataset=/*/*PhaseIISpring17D*/*`

Other resources:
`https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TPhase2Instructions#CMSSW_9_2_0_and_l1t_phase2_v1_10`
`https://twiki.cern.ch/twiki/bin/view/CMS/L1Tracklet90X#Recipe_for_CMSSW_9_2_0`