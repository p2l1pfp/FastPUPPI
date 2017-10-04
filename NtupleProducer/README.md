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

1) Ntuple for single particle and jets response plots and calibrations:

```
cmsRun python/runRespNTupler.py
```

NB: For single particle add "goGun()" at the end of the script, remove it for jets.

To run the ntuplizer over many files do for instance:

```
source python/scripts/prun.sh python/runRespNTupler.py SinglePion_PU0 SinglePion_PU0
```

2) Ntuple for jet HT and MET studies

```
cmsRun python/runJetMetNTupler.py
```

To run the ntuplizer over many files do for instance:

```
source python/scripts/prun.sh python/runJetNTupler.py TTbar_PU140 TTbar_PU140
```

The third step is to produce the plots from the ntuple. The plotting scripts are in:
```FastPUPPI/NtupleProducer/python/scripts```

1) For single particle or jet response:

```
python python/scripts/respPlots.py respTuple_SinglePion_PU0.root plots_dir -w l1pf -p pion
python python/scripts/respPlots.py respTuple_TTbar_PU140.root plots_dir -w l1pf -p jet
```

2) For jet HT plots:

```
python python/scripts/jetHtRateTurnOnPlots.py jetmetTuple_TTbar_PU140.root jetmetTuple_SingleNeutrino_PU140.root plots_dir eff -w l1pf
python python/scripts/jetHtRateTurnOnPlots.py jetmetTuple_TTbar_PU140.root jetmetTuple_SingleNeutrino_PU140.root plots_dir isorate -w l1pf
python python/scripts/jetHtRateTurnOnPlots.py jetmetTuple_TTbar_PU140.root jetmetTuple_SingleNeutrino_PU140.root plots_dir rate -w l1pf
```

3) For MET plots:

```
python python/scripts/met/compareMET.py
```

How to derive the JECs for each algo:

1) run the script ```respCorrSimple.py``` for each algo

```
python python/scripts/respCorrSimple.py respTuple_TTbar_PU140.root plots_dir -p jet -w L1Calo_pt -e L1Calo_pt
python python/scripts/respCorrSimple.py respTuple_TTbar_PU140.root plots_dir -p jet -w L1TK_pt -e L1TK_pt
python python/scripts/respCorrSimple.py respTuple_TTbar_PU140.root plots_dir -p jet -w L1TKV_pt -e L1TKV_pt
python python/scripts/respCorrSimple.py respTuple_TTbar_PU140.root plots_dir -p jet -w L1PF_pt -e L1PF_pt
python python/scripts/respCorrSimple.py respTuple_TTbar_PU140.root plots_dir -p jet -w L1Puppi_pt -e L1Puppi_pt
```

2) and copy the results for the corresponding algo in the runJetMetNTupler.py

The trigger MC can be found on DAS `dataset=/*/*PhaseIISpring17D*/*`

Other resources: <br>
[Correlator code](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TPhase2Instructions#CMSSW_9_2_0_and_l1t_phase2_v1_10) <br>
[Track trigger code](https://twiki.cern.ch/twiki/bin/view/CMS/L1Tracklet90X#Recipe_for_CMSSW_9_2_0) <br>
