Basic Instructions

```
cmsrel CMSSW_10_5_0_pre1
cd CMSSW_10_5_0_pre1/src
cmsenv
git cms-init
git cms-checkout-topic -u p2l1pfp:L1PF_10_5_X

# calibrations
( cd L1Trigger/Phase2L1ParticleFlow/data && \
  rm *.root && \
  git clone git@github.com:p2l1pfp/l1pf-calibrations.git . -b 105X && \
  git checkout 105X_v1 )

git clone git@github.com:p2l1pfp/FastPUPPI.git -b 105X

scram b -j8
```

The first step is to produce the inputs:
```
cd FastPUPPI/NtupleProducer/python/
cmsRun runInputs104X.py # or runInputs93X.py depending on the version of the input files
```
By default the input files contain detailed HGC information, allowing to re-run clustering and to look at individual depths for the towers.
This however takes up a substantial disk space especially at PU 200. 
You can change this by dropping `*_hgcalConcentratorProducer*_*_IN` and `*_hgcalTowerMapProducer_*_IN` in the output commands, or adding a call to `goSlim()` at the bottom of the cfg to do it.

The second step runs the algorithm and creates ntuples which can be used to do analysis:

1) Ntuple for single particle and jets response plots and calibrations:

```
cmsRun python/runRespNTupler.py
```

NB: 
   * For single particle add `goGun()` at the end of the script, remove it for jets.
   * For 93X samples, add a `goOld()` at the end of the script to select 93X calibrations

To run the ntuplizer over many files do for instance:
```
source python/scripts/prun.sh python/runRespNTupler.py --104X TTbar_PU0 TTbar_PU0
source python/scripts/prun.sh python/runRespNTupler.py --104X SinglePion_PU0 SinglePion_PU0  --inline-customize 'goGun()'
```
Look into the prun.sh script to check the paths to the input files and the corresponding options.

2) Ntuple for jet HT and MET studies (BEING REDEFINED, DOCUMENTATION WILL BE UPDATED SOON)

```
cmsRun python/runJetMetNTupler.py
```

To run the ntuplizer over many files do for instance:

```
source python/scripts/prun.sh python/runJetNTuplerNew.py TTbar_PU200 TTbar_PU200
```

The third step is to produce the plots from the ntuple. The plotting scripts are in:
```FastPUPPI/NtupleProducer/python/scripts```

1) For single particle or jet response:

```
python python/scripts/respPlots.py respTupleNew_SinglePion_PU0.root plots_dir -w l1pf -p pion
python python/scripts/respPlots.py respTupleNew_TTbar_PU200.root plots_dir -w l1pf -p jet
```

2) For jet HT plots:

```
python python/scripts/jetHtRateTurnOnPlots.py jetmetTuple_TTbar_PU200.root jetmetTuple_SingleNeutrino_PU200.root plots_dir eff -w l1pf
python python/scripts/jetHtRateTurnOnPlots.py jetmetTuple_TTbar_PU200.root jetmetTuple_SingleNeutrino_PU200.root plots_dir isorate -w l1pf
python python/scripts/jetHtRateTurnOnPlots.py jetmetTuple_TTbar_PU200.root jetmetTuple_SingleNeutrino_PU200.root plots_dir rate -w l1pf
```

3) For MET plots:

```
python python/scripts/met/compareMET.py
```

How to derive the JECs for each algo:

1) run the script ```respCorrSimple.py``` for each algo

```
python python/scripts/respCorrSimple.py respTupleNew_TTbar_PU200.root plots_dir -p jet -w L1Calo_pt -e L1Calo_pt
python python/scripts/respCorrSimple.py respTupleNew_TTbar_PU200.root plots_dir -p jet -w L1TK_pt -e L1TK_pt
python python/scripts/respCorrSimple.py respTupleNew_TTbar_PU200.root plots_dir -p jet -w L1TKV_pt -e L1TKV_pt
python python/scripts/respCorrSimple.py respTupleNew_TTbar_PU200.root plots_dir -p jet -w L1PF_pt -e L1PF_pt
python python/scripts/respCorrSimple.py respTupleNew_TTbar_PU200.root plots_dir -p jet -w L1Puppi_pt -e L1Puppi_pt
```

2) and copy the results for the corresponding algo in the runJetMetNTupler.py

The trigger MC can be found on DAS `dataset=/*/*PhaseIISpring17D*/*`

Other resources: <br>
[Correlator code](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TPhase2Instructions) <br>
[Track trigger code](https://twiki.cern.ch/twiki/bin/view/CMS/L1Tracklet90X) <br>
