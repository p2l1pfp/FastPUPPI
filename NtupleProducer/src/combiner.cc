#include "../interface/combiner.hh"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iostream>
#include <sstream>
#include <TFile.h>


combiner::combiner(const std::string iPionFile,const std::string iElectronFile,const std::string iTrackFile) {
  fDRMatch = 0.1;
  fNEta  = 61;
  loadFile(fPionRes    ,iPionFile);
  loadFile(fElectronRes,iElectronFile);
  loadFile(fTrackRes   ,iTrackFile);
}
void combiner::loadFile(TGraph** &iF1, std::string iFile) { 
  TFile *lFile = new TFile(iFile.c_str());
  iF1 = new TGraph*[fNEta];
  for(int i0 = 0; i0 < fNEta; i0++) { 
    std::stringstream pSS; pSS << "eta_res_" << i0;
    iF1[i0] = (TGraph*) lFile->Get(pSS.str().c_str());
  }
  lFile->Close();
}
//Add Calo deposit
void combiner::addCalo(double iCalo,double iEcal,double iCaloEta,double iCaloPhi,double iEcalEta,double iEcalPhi) { 
  if(iCalo < 1 && iEcal < 1) return;
  double lEt  = iCalo;
  double lEta = iCaloEta;
  double lPhi = iCaloPhi;
  if(iCalo < 1.1*iEcal) lEt  = iEcal;
  if(iCalo < 1.1*iEcal) lEta = iEcalEta;
  if(iCalo < 1.1*iEcal) lPhi = iEcalPhi;
  double lSigma = lEt;
  iCalo < 1.1*iEcal? lSigma = getEleRes(lEt,lEta,lPhi) : lSigma = getPionRes(lEt,lEta,lPhi);
  int lId = 2; if(iCalo < 1.1*iEcal) lId = 3;
  Particle lParticle(lEt,lEta,lPhi,0.,lId,lSigma,0.,lEta,lPhi);
  insert(lParticle,fParticles);
}
//Add Tracks
void combiner::addTrack(double iPt,double iEta,double iPhi,double idZ,double iCaloEta,double iCaloPhi) { 
  if(iPt < 1) return;
  double lSigma = getTrkRes(iPt,iEta,iPhi);
  Particle lParticle(iPt,iEta,iPhi,0.137,0.,lSigma,idZ,iCaloEta,iCaloPhi);
  insert(lParticle,fTkParticles);
}
//Iterate down in Track Pt
void combiner::link() { 
  for(unsigned int i0   = 0; i0 < fTkParticles.size(); i0++) { 
    for(unsigned int i1 = 0; i1 < fParticles.size();   i1++) { 
      if(deltaR(fTkParticles[i0],fParticles[i1]) > fDRMatch) continue;
      if(fParticles[i1].id == 0 || fParticles[i1].id == 1) continue;
      merge(fTkParticles[i0],fParticles[i1],fParticles);
    }
  }
}
//Merge Particles
void combiner::merge(Particle &iTkParticle,Particle &iParticle1,std::vector<Particle> &iCollection) { 
  double lTotSigma = sqrt(iTkParticle.sigma*iTkParticle.sigma + iParticle1.sigma*iParticle1.sigma);
  //Case 1 calo to large
  if(iParticle1.Et-iTkParticle.Et > 2*lTotSigma) { 
    TLorentzVector pVec0; pVec0.SetPtEtaPhiM(iParticle1.Et ,iParticle1.Eta ,iParticle1.Phi ,0.);
    TLorentzVector pVec1; pVec1.SetPtEtaPhiM(iTkParticle.Et,iTkParticle.Eta,iTkParticle.Phi,0.);
    pVec0 -= pVec1;
    iParticle1.Et = pVec0.Pt(); iParticle1.Eta = pVec0.Eta(); iParticle1.Phi = pVec0.Phi(); iParticle1.M = 0;
    iCollection.push_back(iTkParticle);
    return;
  }
  //Case 2 combined 
  if(fabs(iParticle1.Et-iTkParticle.Et) < 2*lTotSigma) { 
    double pSigma   = iParticle1.sigma;
    double pTkSigma = iTkParticle.sigma;
    double pSigTot  = 1./pSigma/pSigma + 1./pTkSigma/pTkSigma;
    double pAvgEt   = (iParticle1.Et/pSigma/pSigma + iTkParticle.Et/pTkSigma/pTkSigma)/pSigTot;
    iParticle1.Et   = pAvgEt;
    iParticle1.Eta  = iTkParticle.Eta;
    iParticle1.Phi  = iTkParticle.Phi;
    if(iParticle1.id == 3) iParticle1.M    = 0.137;
    iParticle1.id -= 2;
    return;
  }
  //Case 3 MIP
  if(iParticle1.Et-iTkParticle.Et < -2*lTotSigma) { 
    iParticle1.Et -= 2; //2 is a bullshit guess at the MIP energy lost
    iCollection.push_back(iTkParticle);
  }
  return;
}
//Order by pt
void combiner::insert(Particle &iParticle,std::vector<Particle> &iParticles) { 
  if(iParticles.size() == 0) iParticles.push_back(iParticle);
  for(std::vector<Particle>::iterator pParticle = iParticles.begin(); pParticle != iParticles.end(); pParticle++) { 
    if(pParticle->Et > iParticle.Et) continue; 
    iParticles.insert(pParticle,iParticle);
    break;
  }
}
//Delta R
double  combiner::deltaR(Particle &iParticle1,Particle &iParticle2) {
  double pDPhi=fabs(iParticle1.caloPhi-iParticle2.caloPhi); 
  if(pDPhi > 2.*TMath::Pi()-pDPhi) pDPhi = 2.*TMath::Pi()-pDPhi;
  double pDEta=iParticle1.caloEta-iParticle2.caloEta;
  return sqrt(pDPhi*pDPhi+pDEta*pDEta);
}
