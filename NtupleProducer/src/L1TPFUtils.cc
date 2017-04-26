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
std::pair<float,float> l1tpf::towerEtaBounds(int ieta) {
  //const float towerEtas[33] = {0,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696,0.783,0.870,0.957,1.044,1.131,1.218,1.305,1.392,1.479,1.566,1.653,1.740,1.830,1.930,2.043,2.172,2.322,2.5,2.650,3.000,3.5,4.0,4.5,5.0}; 
  const float towerEtas[41] = {0,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696,0.783,0.870,0.957,1.044,1.131,1.218,1.305,1.392,1.479,1.566,1.653,1.740,1.830,1.930,2.043,2.172,2.322,2.5,2.650,2.853,3.139,3.314,3.489,3.664,3.839,4.013,4.191,4.363,4.538,4.716,4.889,5.191};
  //float towerEtas[82];
  //for(int i0 = 0; i0 < 82; i0++) towerEtas[i0] = tmpTowerEtas[i0/2]*0.5+tmpTowerEtas[(i0+1)/2]*0.5;
  return std::make_pair( towerEtas[abs(ieta)-1],towerEtas[abs(ieta)] );
}
float l1tpf::towerEtaSize(int ieta) {
  std::pair<float,float> bounds = towerEtaBounds(ieta);
  float size = (bounds.second-bounds.first);
  return size;
}
float l1tpf::towerPhiSize(int ieta) {
  const int kNphi = 72;
  return 2.*M_PI/kNphi;
}
int l1tpf::towerNEta() {
  return 82;
}
int l1tpf::towerNPhi(int ieta) {
  const int kNphi = 72;
  return kNphi;
}
int l1tpf::translateIEta(float eta) { 
  if(eta == 0) return 1;
  const float towerEtas[41] = {0,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696,0.783,0.870,0.957,1.044,1.131,1.218,1.305,1.392,1.479,1.566,1.653,1.740,1.830,1.930,2.043,2.172,2.322,2.5,2.650,2.853,3.139,3.314,3.489,3.664,3.839,4.013,4.191,4.363,4.538,4.716,4.889,5.191};
  int ieta = 0;
  for(int i0 = 1; i0 < 41; i0++) if(fabs(eta) > towerEtas[i0-1] && fabs(eta) < towerEtas[i0]) ieta = i0;
  float sign = eta>0 ? 1. : -1.;
  return sign*ieta;
}
int l1tpf::translateIPhi(float phi,float eta) {
  int iphi = int(phi/towerPhiSize(eta))+0.5;
  if (iphi < 1)                 iphi = iphi + towerNPhi(eta);
  if (iphi > towerNPhi(eta))    iphi = iphi - towerNPhi(eta);
  if (iphi > towerNPhi(eta) || iphi < 0) return 0;
  return iphi;
}
int l1tpf::translateAEta(int ieta,bool iInvert) {
  int lEta = ieta+l1tpf::towerNEta()/2;
  if(iInvert) lEta = ieta-l1tpf::towerNEta()/2;
  return lEta;
}
int l1tpf::translateAPhi(int iphi,bool iInvert) {
  int lPhi = iphi-1;
  if(iInvert) lPhi = iphi+1;
  return lPhi;
}
float l1tpf::towerEta(int ieta) {
  std::pair<float,float> bounds = towerEtaBounds(ieta);
  float eta = (bounds.second+bounds.first)/2.;
  float sign = ieta>0 ? 1. : -1.;
  return sign*eta; 
}
float l1tpf::towerPhi(int ieta, int iphi){
  float phi = (float(iphi)-0.5)*towerPhiSize(ieta);
  if (phi > M_PI) phi = phi - (2*M_PI);
  return phi;
}

