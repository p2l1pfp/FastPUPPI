Basic Instructions

```
cmsrel CMSSW_11_1_0_pre6
cd CMSSW_11_1_0_pre6/src
cmsenv
git cms-init
git cms-checkout-topic -u p2l1pfp:L1PF_11_1_X

# scripts
git clone git@github.com:p2l1pfp/FastPUPPI.git -b 11_1_X

scram b -j8
```

The first step is to produce the inputs:
```
cd FastPUPPI/NtupleProducer/python/
cmsRun runInputs106X.py # 
```
Currently only two kinds of samples are supported:
 * `10_6_X` : L1T TDR, Phase2C8, Geometry D41 using HGCal v10
 * `11_0_X` : HLT TDR, Phase2C9, Geometry D49, HGCal v11.
The older `10_4_X` (MTD TDR) and `9_3_X`samples are no longer supported.

By default the input files contain detailed HGC information, allowing to re-run clustering and to look at individual depths for the towers.
This however takes up a substantial disk space especially at PU 200. 
You can change this by dropping `hgcalConcentratorProducer*` and `hgcalTowerMapProducer` in the output commands, or adding a call to `goSlim()` at the bottom of the cfg to do it.

The second step runs the algorithm and creates ntuples which can be used to do analysis.
All python configuration files for CMSSW are under `NtupleProducer/python`, while standalone python scripts or fwlite macros are under `NtupleProducer/python/scripts`. 
In order to run the python configuration file on many events locally, a driver script `scripts/prun.sh` can be used to run locally the python configuration files, which takes care of selecting the input files, splitting the task to run on multiple CPUs and merge the result.

1) Ntuple for single particle calibrations (and possibly jet studies):

```
cmsRun runRespNTupler.py
```

NB: 
   * For single particle add `goGun()` at the end of the script, remove it for jets.

To run the ntuplizer over many files, from within `NtupleProducer/python` do for instance:
```
./scripts/prun.sh runRespNTupler.py --110X_v0 MultiPion_PT0to200_PU0 MultiPion_PT0to200_PU0.110X_v0  --inline-customize 'goGun()'
./scripts/prun.sh runRespNTupler.py --110X_v0 TTbar_PU0 TTbar_PU0.110X_v0
```
Look into the prun.sh script to check the paths to the input files and the corresponding options.

2) Ntuple for jet HT and MET studies

```
cmsRun runPerformanceNTuple.py
```
This produces both a response ntuple like the one for the runRespNTupler.py (but by default without the detector-level inputs for response calibration) and a NanoAOD-like file with the jets and MET.

To run the ntuplizer over many files do for instance:

```
./scripts/prun.sh runPerformanceNTuple.py --110X_v0  TTbar_PU200 TTbar_PU200.110X_v0 --inline-customize 'addCHS();addTKs()';
```

The third step is to produce the plots from the ntuple. The plotting scripts are in:
```FastPUPPI/NtupleProducer/python/scripts```

1) For single particle or jet response:

```
python scripts/respPlots.py respTupleNew_SinglePion_PU0.110X_v0.root plots_dir -w l1pfw -p pion
python scripts/respPlots.py respTupleNew_TTbar_PU200.110X_v0.root plots_dir -w l1pf -p jet
```

2) For jet, MET and HT plots, the first step is to produce JECs
```
python scripts/makeJecs.py perfNano_TTbar_PU200.root -A -o jecs.root
```
then you can use `jetHtSuite.py` to make the plots

```
python scripts/jetHtSuite.py perfNano_TTbar_PU200.root perfNano_SingleNeutrino_PU200.root plots_dir -w l1pfpu -v ht
python scripts/jetHtSuite.py perfNano_TTbar_PU200.root perfNano_SingleNeutrino_PU200.root plots_dir -w l1pfpu -v met
python scripts/jetHtSuite.py perfNano_TTbar_PU200.root perfNano_SingleNeutrino_PU200.root plots_dir -w l1pfpu -v jet4
```

3) For object multiplicitly studies:

Plot the number of elements in all subdetectors: 
```
python scripts/objMultiplicityPlot.py perfTuple_TTbar_PU200.root  plotdir -s TTbar_PU200  
```
Notes:
 * you should add `goRegional()` to the options in `runPerformanceNTuple.py` if you want to get meaningful numbers per region and not just overall
 * you can select only some subdetector with `-d`, or some kind of object with `-c`.
