#include "../interface/combiner.hh"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iostream>
#include <sstream>
#include <TFile.h>
#include <TH1F.h>


combiner::combiner(const std::string iPionFile,const std::string iElectronFile,const std::string iTrackFile) {
  fDRMatch = 0.2;
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
  if(iCalo < 1.5*iEcal) lEt  = iEcal;
  if(iCalo < 1.5*iEcal) lEta = iEcalEta;
  if(iCalo < 1.5*iEcal) lPhi = iEcalPhi;
  double lSigma = lEt;
  iCalo < 1.5*iEcal ? lSigma = getEleRes(lEt,lEta,lPhi) : lSigma = getPionRes(lEt,lEta,lPhi);
  int lId = 2; if(iCalo < 1.5*iEcal) lId = 3;
  Particle lParticle(lEt,lEta,lPhi,0.,lId,lSigma,0.,0.,lEta,lPhi);
  insert(lParticle,fParticles);
}

//Add Tracks
void combiner::addTrack(double iPt,double iEta,double iPhi,double idZ,double iCaloEta,double iCaloPhi, double iCharge) { 
  if(iPt < 1) return;
  double lSigma = getTrkRes(iPt,iEta,iPhi);
  Particle lParticle(iPt,iEta,iPhi,0.137,0,lSigma,idZ,0,iCaloEta,iCaloPhi,iCharge);
  insert(lParticle,fTkParticles);
}

//Add muons
void combiner::addMuon(double iPt, double iEta, double iPhi, double charge, double quality){
  // for now don't worry about the HCAL MIP matching for muons, just the track matching, so that means that the Calorimeter iEta,iPhi don't matter
  float curPhi = iPhi;
  if (iPhi > TMath::Pi()) curPhi = iPhi - 2*TMath::Pi();
  double lSigma = getTrkRes(iPt,iEta,curPhi); // this needs to be updated with the muon resolutions!
  Particle lParticle(iPt,iEta,curPhi,0.105,4,lSigma,0.,0,iEta,iPhi,charge,quality); // id is 4 for muons
  // insert(lParticle,fMuParticles); // this does something weird to muons
  fMuParticles.push_back(lParticle);
}

//Iterate down in Track Pt
void combiner::link() { 
  // first do the muon/tracking matching
  for(unsigned int i0   = 0; i0 < fMuParticles.size(); i0++) { 
    float minPtDiff = 9999; // find the track that best matches the muon in pT
    for(unsigned int i1 = 0; i1 < fTkParticles.size(); i1++) { 
      float curDR = deltaRraw(fTkParticles[i1],fMuParticles[i0]);
      float curDPT = (fMuParticles[i0].Et/fTkParticles[i1].Et);
      if (curDPT < 1){ // invert it
        float tmptmp = curDPT;
        curDPT = 1/tmptmp;
      }
      if (curDR < 0.2 && (curDPT) < 4 && fabs(curDPT) < minPtDiff){
        //std::cout << "Matched! " << fMuParticles[i0].Et << ", " << fTkParticles[i1].Et << "," << curDR << std::endl;
        fTkParticles[i1].id = 4;
        // fTkParticles[i1].charge = fMuParticles[i0].charge; // use tracks charge, more better like
        minPtDiff = fabs(curDPT);
      }
    } 
  }

  // then do track + calo linking
  for(unsigned int i0   = 0; i0 < fTkParticles.size(); i0++) { 
    if(fTkParticles[i0].id == 4) continue; // skip muons for now, add them at the end
    bool pFill = false;
    for(unsigned int i1 = 0; i1 < fParticles.size();   i1++) { 
      // what happens if there is no matching cluster?  does it throw out the track? it should to reduce fake tracks
      if(deltaR(fTkParticles[i0],fParticles[i1]) > fDRMatch) continue;
      if(fParticles[i1].id == 0 || fParticles[i1].id == 1)   continue;
      merge(fTkParticles[i0],fParticles[i1],fParticles);
      pFill = true;
    }
    //Remove high pT fakes
    if(fTkParticles[i0].Et > 10.) pFill = true;
    if(!pFill) fParticles.push_back(fTkParticles[i0]);
  }

  // now do muons...
  for(unsigned int i0   = 0; i0 < fTkParticles.size(); i0++) { 
    if (fTkParticles[i0].id == 4) fParticles.push_back(fTkParticles[i0]);
  }

}

//Merge Particles
void combiner::merge(Particle &iTkParticle,Particle &iParticle1,std::vector<Particle> &iCollection) { 
  double lTotSigma = sqrt(iTkParticle.sigma*iTkParticle.sigma + iParticle1.sigma*iParticle1.sigma);
  //Case 1 calo to large
  if(iParticle1.Et-iTkParticle.Et > 3.*lTotSigma) { 
    TLorentzVector pVec0; pVec0.SetPtEtaPhiM(iParticle1.Et ,iParticle1.Eta ,iParticle1.Phi ,0.);
    TLorentzVector pVec1; pVec1.SetPtEtaPhiM(iTkParticle.Et,iTkParticle.Eta,iTkParticle.Phi,0.);
    pVec0 -= pVec1;
    iParticle1.Et = pVec0.Pt(); iParticle1.Eta = pVec0.Eta(); iParticle1.Phi = pVec0.Phi(); iParticle1.M = 0;
    iCollection.push_back(iTkParticle);
    return;
  }
  //Case 2 combined 
  if(fabs(iParticle1.Et-iTkParticle.Et) < 3.*lTotSigma) { 
    double pSigma   = iParticle1.sigma;
    double pTkSigma = iTkParticle.sigma;
    double pSigTot  = 1./pSigma/pSigma + 1./pTkSigma/pTkSigma;
    double pAvgEt   = (iParticle1.Et/pSigma/pSigma + iTkParticle.Et/pTkSigma/pTkSigma)/pSigTot;
    iParticle1.Et   = pAvgEt;
    iParticle1.Eta  = iTkParticle.Eta;
    iParticle1.Phi  = iTkParticle.Phi;
    if(iParticle1.id == 3) iParticle1.M    = 0.137;
    iParticle1.id -= 2;
    iParticle1.charge = iTkParticle.charge;
    return;
  }
  //Case 3 MIP
  if(iParticle1.Et-iTkParticle.Et < -3*lTotSigma) { 
    //Now a cut to remove fake tracks
    iCollection.push_back(iTkParticle);
    iParticle1.Et -= 2; //2 is a bullshit guess at the MIP energy lost
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

//Delta R
double  combiner::deltaRraw(Particle &iParticle1,Particle &iParticle2) {
  // std::cout << iParticle1.Phi << "," << iParticle2.Phi << ", " << iParticle1.Eta << "," << iParticle2.Eta << std::endl;
  double pDPhi=fabs(iParticle1.Phi-iParticle2.Phi); 
  if(pDPhi > 2.*TMath::Pi()-pDPhi) pDPhi = 2.*TMath::Pi()-pDPhi;
  double pDEta=iParticle1.Eta-iParticle2.Eta;
  return sqrt(pDPhi*pDPhi+pDEta*pDEta);
}

void combiner::doVertexing(){
  fDZ  = -999;
  // std::vector<Particle> fTkParticlesWVertexing;

  TH1F *h_dz = new TH1F("h_dz","h_dz",200,-20,20); // 0.5mm binning
  for (int i = 0; i < h_dz->GetXaxis()->GetNbins(); ++i) h_dz->SetBinContent(i+1,0.); //initialize all to 0.

  // std::cout << "------ fTkParticles.size(): " << fTkParticles.size() << std::endl;
  for(unsigned int i0   = 0; i0 < fTkParticles.size(); i0++) { 
    // std::cout << "fTkParticles[i0].dZ " << i0 << ": " << fTkParticles[i0].dZ << "," << fTkParticles[i0].Et << std::endl;
    float curdz  = fTkParticles[i0].dZ;
    float curbin = h_dz->GetXaxis()->FindBin(curdz);
    h_dz->SetBinContent( curbin, h_dz->GetBinContent(curbin) + fTkParticles[i0].Et );
  }

  int imaxbin = h_dz->GetMaximumBin();
  float pvdz = h_dz->GetXaxis()->GetBinCenter(imaxbin);
  float binwidth = h_dz->GetXaxis()->GetBinWidth(imaxbin);
  float pvdz_lo = pvdz - 1.5*binwidth;
  float pvdz_hi = pvdz + 1.5*binwidth;
  fDZ = pvdz;

  for(unsigned int i0   = 0; i0 < fTkParticles.size(); i0++) { 
    float curdz  = fTkParticles[i0].dZ;
    if( (curdz < pvdz_hi && curdz > pvdz_lo)) fTkParticles[i0].pvid = 1;
    if(!(curdz < pvdz_hi && curdz > pvdz_lo)) fTkParticles[i0].pvid = 2;
    if (curdz < pvdz_hi && curdz > pvdz_lo)  fTkParticlesWVertexing.push_back( fTkParticles[i0] );
    if(fTkParticles[i0].id == 4 && fTkParticles[i0].Et > 20) std::cout << " ---> " << curdz << " --> " << pvdz_hi << " -- " << pvdz_lo << std::endl;
  }

  for(unsigned int i0   = 0; i0 < fParticles.size(); i0++) { 
    float curdz  = fParticles[i0].dZ;
    int   curid  = fParticles[i0].id;
    if( (curdz < pvdz_hi && curdz > pvdz_lo) || curid == 2 || curid == 3)  fCHSParticles.push_back( fParticles[i0] );
    if( (curdz < pvdz_hi && curdz > pvdz_lo) && curid != 2 && curid != 3)  fParticles[i0].pvid = 1;
    if(!(curdz < pvdz_hi && curdz > pvdz_lo) && curid != 2 && curid != 3)  fParticles[i0].pvid = 2;
    if(                                         curid == 2 || curid == 3)  fParticles[i0].pvid = 0;
  }
}
