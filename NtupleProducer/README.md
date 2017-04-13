Basic Instructions

```
cmsrel CMSSW_9_1_0_pre2
cd CMSSW_9_1_0_pre2/src
cmsenv
scram b -j8
```

The current basic config is:
`FastPUPPI/NtupleProducer/python/runNtupleProducer_cfg.py`

Relval samples are located in:
`/store/relval/CMSSW_9_1_0_pre1/RelValZMM_14/`
n.b. the \_14 samples are for upgrade, and they are not in every pre-release so we are using pre1 instead of pre2
