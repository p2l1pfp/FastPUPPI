Basic Instructions

```
cmsrel CMSSW_10_6_1_patch2
cd CMSSW_10_6_1_patch2/src
cmsenv
git cms-init
git cms-checkout-topic -u p2l1pfp:L1PF_10_6_1p2_X

# scripts
git clone git@github.com:p2l1pfp/FastPUPPI.git -b 106X

scram b -j8
```

The first step is to produce the inputs:
```
cd FastPUPPI/NtupleProducer/python/
cmsRun runInputs104X.py # or runInputs106X.py or runInputs93X.py depending on the version of the input files
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

To run the ntuplizer over many files, from within `NtupleProducer/python` do for instance:
```
./scripts/prun.sh runRespNTupler.py --v0 TTbar_PU0 TTbar_PU0.v0
./scripts/prun.sh runRespNTupler.py --v0 ParticleGun_PU0 ParticleGun_PU0.v0  --inline-customize 'goGun()'
```
Look into the prun.sh script to check the paths to the input files and the corresponding options.

2) Ntuple for jet HT and MET studies

```
cmsRun runPerformanceNTuple.py
```
This produces both a response ntuple like the one for the runRespNTupler.py, but by default without the inputs for response calibration, and a NanoAOD-like file with the jets and MET.

To run the ntuplizer over many files do for instance:

```
./scripts/prun.sh runPerformanceNTuple.py --v0  TTbar_PU200 TTbar_PU200.v0
```

The third step is to produce the plots from the ntuple. The plotting scripts are in:
```FastPUPPI/NtupleProducer/python/scripts```

1) For single particle or jet response:

```
python scripts/respPlots.py respTupleNew_SinglePion_PU0.v0.root plots_dir -w l1pfw -p pion
python scripts/respPlots.py respTupleNew_TTbar_PU200.v0.root plots_dir -w l1pf -p jet
```

2) For jet, MET and HT plots, the first step is to produce JECs
```
python scripts/makeJecs.py perfNano_TTbar_PU200.v0.root -A -o jecs.root
```
then you can use `jetHtSuite.py` to make the plots

```
python scripts/jetHtSuite.py perfNano_TTbar_PU200.root perfNano_SingleNeutrino_PU200.root plots_dir -w l1pfpu -v ht
python scripts/jetHtSuite.py perfNano_TTbar_PU200.root perfNano_SingleNeutrino_PU200.root plots_dir -w l1pfpu -v met
python scripts/jetHtSuite.py perfNano_TTbar_PU200.root perfNano_SingleNeutrino_PU200.root plots_dir -w l1pfpu -v jet4
```
