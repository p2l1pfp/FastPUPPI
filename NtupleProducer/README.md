Basic Instructions

```
cmsrel CMSSW_10_5_0_pre1
cd CMSSW_10_5_0_pre1/src
cmsenv
git cms-init
git cms-checkout-topic -u p2l1pfp:L1PF_10_5_X

# calibrations
git clone git@github.com:p2l1pfp/l1pf-calibrations.git L1Trigger/Phase2L1ParticleFlow/data -b 105X

# scripts
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
You can change this by dropping `hgcalConcentratorProducer*` and `hgcalTowerMapProducer` in the output commands, or adding a call to `goSlim()` at the bottom of the cfg to do it.

The second step runs the algorithm and creates ntuples which can be used to do analysis.
All python configuration files for CMSSW are under `NtupleProducer/python`, while standalone python scripts or fwlite macros are under `NtupleProducer/python/scripts`.

1) Ntuple for single particle and jets response plots and calibrations:

```
cmsRun runRespNTupler.py
```

NB: 
   * For single particle add `goGun()` at the end of the script, remove it for jets.
   * For 93X samples, add a `goOld()` at the end of the script to select 93X calibrations

To run the ntuplizer over many files, from within `NtupleProducer/python` do for instance:
```
./scripts/prun.sh runRespNTupler.py --104X TTbar_PU0 TTbar_PU0
./scripts/prun.sh runRespNTupler.py --104X SinglePion_PU0 SinglePion_PU0  --inline-customize 'goGun()'
```
Look into the prun.sh script to check the paths to the input files and the corresponding options.

2) Ntuple for jet HT and MET studies

```
cmsRun runPerformanceNTuple.py
```
This produces both a response ntuple like the one for the runRespNTupler.py, but by default without the inputs for response calibration, and a NanoAOD-like file with the jets and MET.

To run the ntuplizer over many files do for instance:

```
./scripts/prun.sh runPerformanceNTuple.py TTbar_PU200 TTbar_PU200
```

The third step is to produce the plots from the ntuple. The plotting scripts are in:
```FastPUPPI/NtupleProducer/python/scripts```

1) For single particle or jet response:

```
python scripts/respPlots.py respTupleNew_SinglePion_PU0.root plots_dir -w l1pfnew -p pion
python scripts/respPlots.py respTupleNew_TTbar_PU200.root plots_dir -w l1pfnew -p jet
```

2) For jet, MET and HT plots, the first step is to produce JECs
```
python scripts/makeJecs.py perfNano_TTbar_PU200.root -A -o jecs.root
```
then you can use `jetHtSuite.py` to make the plots

```
python scripts/jetHtSuite.py perfNano_TTbar_PU200.root perfNano_SingleNeutrino_PU200.root plots_dir -w oldcomp,newcomp -v ht
python scripts/jetHtSuite.py perfNano_TTbar_PU200.root perfNano_SingleNeutrino_PU200.root plots_dir -w oldcomp,newcomp -v met
python scripts/jetHtSuite.py perfNano_TTbar_PU200.root perfNano_SingleNeutrino_PU200.root plots_dir -w oldcomp,newcomp -v jet4
```
