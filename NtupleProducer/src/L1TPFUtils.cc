#include <FastPUPPI/NtupleProducer/interface/L1TPFUtils.h>

#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "FastSimulation/Particle/interface/RawParticle.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

void l1tpf::propagate(int iOption,std::vector<double> &iVars,const math::XYZTLorentzVector& iMom,const math::XYZTLorentzVector& iVtx,double iCharge,double iBField) {
    BaseParticlePropagator particle = BaseParticlePropagator(RawParticle(iMom,iVtx),0.,0.,iBField);
    particle.setCharge(iCharge);
    double ecalShowerDepth=0;
    if(iOption == 0 || iOption == 1) particle.propagateToEcalEntrance(false);
    if(iOption == 2 || iOption == 3) particle.propagateToHcalEntrance(false);
    if(iOption == 4) particle.propagateToHcalExit    (false);
    if(iOption == 1) ecalShowerDepth = reco::PFCluster::getDepthCorrection(particle.momentum().E(),false,false);
    math::XYZVector point = math::XYZVector(particle.vertex())+math::XYZTLorentzVector(particle.momentum()).Vect().Unit()*ecalShowerDepth;
    iVars.push_back(particle.momentum().px());
    iVars.push_back(particle.momentum().py());
    iVars.push_back(particle.momentum().pz());
    iVars.push_back(particle.momentum().energy());
    iVars.push_back(point.eta());
    iVars.push_back(point.phi());
    iVars.push_back(point.rho());
}

