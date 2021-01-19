#!/bin/bash

V="110X_v2.0"
PLOTDIR="plots/11_1_X/from110X/${V}/corr"
SAMPLES="--110X_v2";

PU=PU0
if [[ "$1" == "--pu200" ]]; then 
    PU=PU200; shift;
fi;

TUPLE="respTupleNew";
if [[ "$1" == "--perf" ]]; then
    TUPLE="perfTuple"; shift;
fi




if [[ "$1" == "--backup" ]]; then # ,DoublePhoton_FlatPt-1To100
    NTUPLES=$(ls ${TUPLE}_{DoubleElectron_FlatPt-1To100,DoublePhoton_FlatPt-1To100,MultiPion_PT0to200,K0Long_PT0to200_gun}_${PU}.${V}.0.root); 
    shift;
else
    NTUPLES=$(ls ${TUPLE}_{DoubleElectron_FlatPt-1To100,DoublePhoton_FlatPt-1To100,MultiPion_PT0to200,K0Long_PT0to200_gun}_${PU}.${V}.root); 
fi;

W=$1; shift;
case $W in
    em-barrel)
         python scripts/respCorrSimple.py $NTUPLES $PLOTDIR/$W -p emmix -w L1RawBarrelEcal_pt02 -e L1RawBarrelEcal_pt02 --fitrange 15 80 --root emcorr_barrel_110X.root  --barrel-eta --ptmax 120 && \
         cp -v emcorr_barrel_110X.root $CMSSW_BASE/src/L1Trigger/Phase2L1ParticleFlow/data 
         ;;
    em-hgc) # currenlt only used for old-style PF for the moment
         python scripts/respCorrSimple.py $NTUPLES $PLOTDIR/$W -p emmix  -w L1RawHGCalEM_pt -e L1RawHGCalEM_pt --fitrange 15 80 --root emcorr_hgc_110X.root --hgcal-eta --ptmax 120 && \
         cp -v emcorr_hgc_110X.root  $CMSSW_BASE/src/L1Trigger/Phase2L1ParticleFlow/data
         ;;
    had-barrel)
         python scripts/respCorrSimple.py $NTUPLES $PLOTDIR/$W -p pion -w L1RawBarrelCalo_pt -e L1RawBarrelCalo_pt --fitrange 15 80  --root hadcorr_barrel_110X.root \
                --emf-slices L1RawBarrelCaloEM_pt/L1RawBarrelCalo_pt 0.125,0.50,0.875,1.125 --ptmax 150  --barrel-eta && \
         cp -v hadcorr_barrel_110X.root  $CMSSW_BASE/src/L1Trigger/Phase2L1ParticleFlow/data
         ;;
    hgc)
         python scripts/respCorrSimple.py $NTUPLES $PLOTDIR/$W -p mix  -w L1RawHGCal_pt -e L1RawHGCal_pt --fitrange 15 100 --root hadcorr_HGCal3D_TC_110X.root \
                    --emf-slices L1RawHGCalEM_pt02/L1RawHGCal_pt02 0.125,0.250,0.50,0.875,1.125  --ptmax 100 --hgcal-eta && \
         cp -v hadcorr_HGCal3D_TC_110X.root  $CMSSW_BASE/src/L1Trigger/Phase2L1ParticleFlow/data
         ;;
    hf)
         python scripts/respCorrSimple.py  $NTUPLES $PLOTDIR/$W  -p mixmix -w L1RawHFCalo_pt -e L1RawHFCalo_pt  --root hfcorr_110X.root  --ptmin 5 --hf-eta-fine --ptmax 100 && \
         cp -v hfcorr_110X.root $CMSSW_BASE/src/L1Trigger/Phase2L1ParticleFlow/data
         ;;
    res-em-barrel)
         python scripts/respCorrSimple.py $NTUPLES $PLOTDIR/$W -p emmix -w L1BarrelEcal_pt02 -e L1BarrelEcal_pt02 --fitrange 15 80  --barrel-eta -r && \
         echo "Put this into pfClustersFromL1EGClusters_cfi.py" 
         ;;
    res-had-barrel)
         python scripts/respCorrSimple.py $NTUPLES $PLOTDIR/$W -p pion -w L1BarrelCalo_pt -e L1BarrelCalo_pt --fitrange 15 80  --barrel-eta -r && \
         echo "Put this into l1ParticleFlow_cff.py under pfClustersFromCombinedCaloHCal" 
         ;;
    res-em-hgc) # only used for old-style PF
         python scripts/respCorrSimple.py $NTUPLES $PLOTDIR/$W -p emmix -w L1HGCalEM_pt -e L1HGCalEM_pt --fitrange 15 80  --hgcal-eta -r && \
         echo "Put this into pfClustersFromHGC3DClustersEM_cfi.py" 
         ;;
    res-hgc)
         python scripts/respCorrSimple.py $NTUPLES $PLOTDIR/$W -p pion -w L1HGCal_pt -e L1HGCal_pt --fitrange 15 150  --hgcal-eta -r && \
         echo "Put this into pfClustersFromHGC3DClusters_cfi.py" 
         ;;
    res-hf)
         python scripts/respCorrSimple.py $NTUPLES $PLOTDIR/$W -p pimix -w L1HFCalo_pt -e L1HFCalo_pt --fitrange 15 60  --hf-eta -r && \
         echo "Put this into l1ParticleFlow_cff.py under pfClustersFromCombinedCaloHF" 
         ;;
    backup)
         for f in $NTUPLES; do cp -v $f ${f/$V/$V.0}; done
         ;;
    rerun)
         python runRespNTupler.py || exit 1;
         for f in $NTUPLES; do test -f $f && mv -v $f ${f/$V/$V.0}; done
         for X in {MultiPion_PT0to200,DoublePhoton_FlatPt-1To100,DoublePhoton_FlatPt-1To100,DoubleElectron_FlatPt-1To100,K0Long_PT0to200_gun}_${PU}; do
             ./scripts/prun.sh runRespNTupler.py $SAMPLES $X ${X}.${V}  --inline-customize 'goGun();noPU();' ;  # --maxfiles 999;
             grep -q '^JOBS:.*, 0 passed' respTupleNew_${X}.${V}.report.txt && echo " === ERROR. Will stop here === " && break || true;
         done
         ;;
    plot-ecal)
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -p electron,photon --gauss -w debug-ecal --eta 0 1.3 --no-eta --ymaxRes 0.35 --ptmax 150 --ptdef pt02,ptbest
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -p electron,photon         -w debug-ecal --eta 0 1.7 --no-pt  --ptdef pt02,ptbest
         ;;
    plot-hcal)
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -p photon,pion --gauss -w debug-hcal --eta 0 1.3 --no-eta --ymax 2.1 --ymaxRes 1.2 --ptmax 150 --ptdef pt02,ptbest
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -p pion                -w debug-hcal --eta 0 1.6 --no-pt  --ymax 2.1 --ptdef pt02,ptbest
         ;;
    plot-hgc)
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -p photon,electron --gauss -w debug-hgc --eta 1.7 2.5 --no-eta --ymax 1.8 --ymaxRes 0.35 --ptmax 100 --ptdef pt02,ptbest
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -p photon,electron --gauss -w debug-hgc --eta 2.5 3.0 --no-eta --ymax 1.8 --ymaxRes 0.35 --ptmax 100 --ptdef pt02,ptbest
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -p pion,klong,mixmix  --gauss -w debug-hgc --eta 1.7 2.5 --no-eta --ymax 1.8 --ymaxRes 1.2  --ptmax 100 --ptdef pt02,ptbest
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -p pion,klong,mixmix  --gauss -w debug-hgc --eta 2.5 3.0 --no-eta --ymax 1.8 --ymaxRes 1.2  --ptmax 100 --ptdef pt02,ptbest
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -p photon,electron,pion,klong,mixmix  -w debug-hgc --eta 1.3 3.3 --no-pt  --ymax 1.8 --ptdef pt02,ptbest
         ;;
    plot-hf)
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -p photon,pion,mixmix --gauss -w debug-hf --eta 3.0 3.5  --ymax 1.8 --ymaxRes 1.5 --ptmax 100 --ptdef pt02,ptbest --no-eta
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -p photon,pion,mixmix --gauss -w debug-hf --eta 3.0 3.2  --ymax 1.8 --ymaxRes 1.5 --ptmax 100 --ptdef pt02,ptbest --no-eta
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -p photon,pion,mixmix --gauss -w debug-hf --eta 3.2 3.5  --ymax 1.8 --ymaxRes 1.5 --ptmax 100 --ptdef pt02,ptbest --no-eta
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -p photon,pion,mixmix --gauss -w debug-hf --eta 3.5 4.5  --ymax 1.8 --ymaxRes 1.5 --ptmax 100 --ptdef pt02,ptbest --no-eta
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -p photon,pion,mixmix --gauss -w debug-hf --eta 4.5 5.0  --ymax 1.8 --ymaxRes 1.5 --ptmax 100 --ptdef pt02,ptbest --no-eta
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -p photon,pion,mixmix --gauss -w debug-hf --eta 3.0 5.0  --ymax 1.8  --ptdef pt02,ptbest --no-pt
         ;;
    plots-pf)
         if [[ "$PU" == "PU0" ]]; then
             python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pf -p electron,photon  -g  --ymaxRes 0.35 --ptmax 100 -E 3.0
             python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pf -p pion,mixmix,klong      -g  --ymaxRes 1.2  --ptmax 100 -E 3.0
             python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pf -p pion,photon,mixmix -g  --eta 3.0 5.0  --ymax 3 --ymaxRes 1.5 --label hf  --no-fit 
             python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pf -p pion,photon,mixmix -g  --eta 3.5 4.5  --ymax 3 --ymaxRes 1.5 --label hf  --no-fit --no-eta
         else
             python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pf -p electron,photon,pizero -g  --ymax 2.5 --ymaxRes 0.6 --ptmax 80 -E 3.0
             python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pf -p pion,klong,mixmix       -g  --ymax 2.5 --ymaxRes 1.5 --ptmax 80 -E 3.0
             python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pf -p pion,pizero,mixmix      -g  --ymax 3.0 --ymaxRes 1.5 --ptmax 80 --eta 3.0 5.0 --label hf  --no-fit 
         fi;
         ;;
    plots-pfcomp)
         if [[ "$PU" == "PU0" ]]; then
             python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pfcomp -p electron,photon,pizero -g  --ymaxRes 0.35 --ptmax 150 -E 3.0
             python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pfcomp -p pion,klong,mixmix       -g  --ymaxRes 1.2  --ptmax 150 -E 3.0
             python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pfcomp -p electron,photon,pizero -g  --ymaxRes 0.35 --ptmax 50 -E 3.0 --label lowpt
             python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pfcomp -p pion,klong,mixmix       -g  --ymaxRes 1.2  --ptmax 50 -E 3.0 --label lowpt
             python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pfcomp -p pion,pizero,mixmix      -g  --eta 3.0 5.0  --ymax 3 --ymaxRes 1.5 --label hf  --no-fit 
         else
             python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pfcomp -p electron,photon,pizero -g  --ymax 2.5 --ymaxRes 0.6 --ptmax 80 -E 3.0
             python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pfcomp -p pion,klong,mixmix       -g  --ymax 2.5 --ymaxRes 1.5 --ptmax 80 -E 3.0
             python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pfcomp -p pion,pizero,mixmix      -g  --ymax 3.0 --ymaxRes 1.5 --ptmax 80 --eta 3.0 5.0 --label hf  --no-fit 
         fi;
         ;;
    plots-pfold)
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pfold -p electron,photon,pizero -g  --ymaxRes 0.35 --ptmax 150 -E 3.0
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pfold -p pion,klong,mixmix       -g  --ymaxRes 1.2  --ptmax 150 -E 3.0
         python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pfold -p pion,pizero,mixmix      -g  --eta 3.0 5.0  --ymax 3 --ymaxRes 1.5 --label hf  --no-fit 
         ;;

esac;
