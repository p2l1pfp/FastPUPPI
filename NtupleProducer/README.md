Basic Instructions

```
cmsrel CMSSW_9_2_0
cd CMSSW_9_2_0/src
cmsenv
git cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline
git cms-merge-topic -u cms-l1t-offline:l1t-phase2-v1.14.1
git clone git@github.com:p2l1pfp/FastPUPPI.git -b 92X
scram b -j8
```

The first step is to produce the inputs:
```
cd FastPUPPI/NtupleProducer/python/runInputs.py
cmsRun runInputs.py
```

Example how to submit jobs to the lxbatch:
```
cd FastPUPPI/NtupleProducer/python
./scripts/cmsSplit.pl --dataset /SingleE_FlatPt-8to100/PhaseIISpring17D-NoPU_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW --jobs 1 --events-per-job 10 --label test runInputs.py --lsf 1nh --eosoutdir /eos/cms/store/cmst3/user/jngadiub/L1PFInputs/
./runInputs_test_bsub.sh
```
For more options for job splitting:
```
./scripts/cmsSplit.pl --help
```
The job outputs are copied to the eos directory --eosoutdir

When jobs are finished you can clean up your local dir
```
./runInputs_test_cleanup.sh
```

The second step runs the algorithm and creates ntuples which can be used to do analysis:
```
cmsRun runNtupleProducer_cfg.py
cmsRun runJetNTupler.py
cmsRun runRespNTupler.py
```
`Should we try to collate all of these?`

The plotting scripts are in:
```FastPUPPI/NtupleProducer/python/scripts```


The trigger MC can be found on DAS `dataset=/*/*PhaseIISpring17D*/*`

Other resources: <br>
[Correlator code](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TPhase2Instructions#CMSSW_9_2_0_and_l1t_phase2_v1_10) <br>
[Track trigger code](https://twiki.cern.ch/twiki/bin/view/CMS/L1Tracklet90X#Recipe_for_CMSSW_9_2_0) <br>
