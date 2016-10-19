#include "../interface/combiner.hh"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iostream>
#include <sstream>
#include <TFile.h>
#include <TH1F.h>
#include <algorithm>
#include "TMath.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"

combiner::combiner(const std::string iPionFile,const std::string iElectronFile,const std::string iTrackFile,std::string iFile = "puppi.root") {
  fDRMatch = 0.15;
  fNEta  = 61;
  loadFile(fPionRes    ,iPionFile);
  loadFile(fElectronRes,iElectronFile);
  loadFile(fTrackRes   ,iTrackFile);

  fFile = new TFile(iFile.c_str(),"RECREATE");
  fTree = new TTree("puppi","puppi");
  fTree->Branch("pt",&b_pt);
  fTree->Branch("alphaF",&b_alphaF);
  fTree->Branch("alphaC",&b_alphaC);
  fTree->Branch("eta",&b_Eta);
  fTree->Branch("phi",&b_Phi);
  fTree->Branch("et",&b_Et);
  fTree->Branch("puppiWeight",&b_PuppiWeight);
  fTree->Branch("alphaFMed",&b_alphaFMed);
  fTree->Branch("alphaCMed",&b_alphaCMed);
  fTree->Branch("alphaFRms",&b_alphaFRms);
  fTree->Branch("alphaCRms",&b_alphaCRms);
  fTree->Branch("nPars",&b_nParticles,"nPars/I");
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
  Particle lParticle(lEt,lEta,lPhi,0.,lId,lSigma,0.,lEta,lPhi);
  insert(lParticle,fParticles);
}

//Add Tracks
void combiner::addTrack(double iPt,double iEta,double iPhi,double idZ,double iCaloEta,double iCaloPhi, double iCharge) { 
  if(iPt < 1) return;
  double lSigma = getTrkRes(iPt,iEta,iPhi);
  Particle lParticle(iPt,iEta,iPhi,0.137,0,lSigma,idZ,iCaloEta,iCaloPhi,iCharge);
  insert(lParticle,fTkParticles);
}

//Add muons
void combiner::addMuon(double iPt, double iEta, double iPhi, double charge, double quality){
  // for now don't worry about the HCAL MIP matching for muons, just the track matching, so that means that the Calorimeter iEta,iPhi don't matter
  float curPhi = iPhi;
  if (iPhi > TMath::Pi()) curPhi = iPhi - 2*TMath::Pi();
  double lSigma = getTrkRes(iPt,iEta,curPhi); // this needs to be updated with the muon resolutions!
  Particle lParticle(iPt,iEta,curPhi,0.105,4,lSigma,0.,0,0,charge,quality); // id is 4 for muons
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
    if(fTkParticles[i0].Et > 300.) pFill = true;
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
    TLorentzVector pVec1; pVec1.SetPtEtaPhiM(iTkParticle.Et,iTkParticle.caloEta,iTkParticle.caloPhi,0.);
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
  if(iParticles.size() == 0) {iParticles.push_back(iParticle); return;}
  bool lFill = false;
  for(std::vector<Particle>::iterator pParticle = iParticles.begin(); pParticle != iParticles.end(); pParticle++) { 
    if(pParticle->Et > iParticle.Et) continue; 
    iParticles.insert(pParticle,iParticle);
    lFill = true;
    break;
  }
  if(!lFill) iParticles.push_back(iParticle);
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
  
  // std::vector<Particle> fTkParticlesWVertexing;
  
  TH1F *h_dz = new TH1F("h_dz","h_dz",40,-20,20); // 1cm binning
  for (int i = 0; i < h_dz->GetXaxis()->GetNbins(); ++i) h_dz->SetBinContent(i+1,0.); //initialize all to 0.

  // std::cout << "------ fTkParticles.size(): " << fTkParticles.size() << std::endl;
  int ntracks = 0;
  for(unsigned int i0   = 0; i0 < fParticles.size(); i0++) { 
    if (fParticles[i0].charge == 0) continue;
    ntracks++;
    // std::cout << i0 << ": fParticles[i0].dZ = " << fParticles[i0].dZ << "," << fParticles[i0].Et << ", " << fParticles[i0].id << std::endl;
    float curdz  = fParticles[i0].dZ;
    float curbin = h_dz->GetXaxis()->FindBin(curdz);
    h_dz->SetBinContent( curbin, h_dz->GetBinContent(curbin) + fParticles[i0].Et );
  }

  int imaxbin = h_dz->GetMaximumBin();
  float pvdz = h_dz->GetXaxis()->GetBinCenter(imaxbin);
  float binwidth = h_dz->GetXaxis()->GetBinWidth(imaxbin);
  float pvdz_lo = pvdz - 1.5*binwidth;
  float pvdz_hi = pvdz + 1.5*binwidth;

  for(unsigned int i0   = 0; i0 < fParticles.size(); i0++) { 
    float curdz  = fParticles[i0].dZ;
    if (curdz < pvdz_hi && curdz > pvdz_lo){
      fTkParticlesWVertexing.push_back( fParticles[i0] );
      fParticles[i0].isPV = 1;
    }
  }

  // std::cout << "ntracks = " << ntracks << ", " << fTkParticlesWVertexing.size() << ", " << fTkParticles.size() << std::endl;
  // return fTkParticlesWVertexing;
}


////////////////////////////////////////////////////////////////////////////
// PUPPI TIME!
////////////////////////////////////////////////////////////////////////////

void combiner::fetchPuppi(){

  // std::cout << "fetching puppi..." << std::endl;

  // separate particles out into neutrals and charged PU and charged PV
  std::vector<Particle> puppi_neutrals;
  std::vector<Particle> puppi_chargedPU;
  std::vector<Particle> puppi_chargedPV;
  for(unsigned int i0   = 0; i0 < fParticles.size(); i0++) { 
    if (fParticles[i0].charge == 0) puppi_neutrals.push_back(fParticles[i0]);
    if (fParticles[i0].charge != 0 && fParticles[i0].isPV != 0 && fParticles[i0].id != 4) puppi_chargedPV.push_back(fParticles[i0]);
    if (fParticles[i0].charge != 0 && fParticles[i0].isPV == 0) puppi_chargedPU.push_back(fParticles[i0]);
  }
  std::cout << fParticles.size() << ", " << puppi_neutrals.size() << ", " << puppi_chargedPU.size() << ", " << puppi_chargedPV.size() << std::endl;

  // compute alphas (not efficient, but we don't have that many particles)
  computeAlphas(fParticles,0);
  computeAlphas(puppi_chargedPV,1);

  // compute median and RMS
  computeMedRMS();

  // compute weights 
  computeWeights();

}

void combiner::computeAlphas(std::vector<Particle> neighborParticles, int isCentral){

  for(unsigned int i0   = 0; i0 < fParticles.size(); i0++) { 
    float curalpha = 0;
    for (unsigned int j0 = 0; j0 < neighborParticles.size(); j0++) {

      float curdr = deltaRraw(fParticles[i0],neighborParticles[j0]);
      if (curdr < 0.5 && curdr != 0.) curalpha += neighborParticles[j0].Et * neighborParticles[j0].Et / curdr / curdr;

    }

    float logcuralpha = -99;
    if (curalpha != 0.) logcuralpha = log(curalpha);
    if (isCentral == 0) fParticles[i0].alphaF = logcuralpha;
    if (isCentral == 1) fParticles[i0].alphaC = logcuralpha;
  }

}

void combiner::computeMedRMS(){

  // std::cout << "computing alpha med/rms..." << std::endl;
  std::vector<float> v_alphaFs;
  std::vector<float> v_alphaCs;
  for(unsigned int i0   = 0; i0 < fParticles.size(); i0++) { 
    // if (fParticles[i0].charge != 0 && fParticles[i0].isPV == 0){
    // if (fParticles[i0].alphaF > -98) v_alphaFs.push_back(fParticles[i0].alphaF);
    // if (fParticles[i0].alphaC > -98) v_alphaCs.push_back(fParticles[i0].alphaC);

    //Use non-pv tracks for alphaF in central
    if (fabs(fParticles[i0].Eta) < 2.5 && (fParticles[i0].id == 2 || fParticles[i0].id == 3 || fParticles[i0].isPV == 1)) continue;
    if (fParticles[i0].alphaF > -98 && fabs(fParticles[i0].Eta) > 2.5) v_alphaFs.push_back(fParticles[i0].alphaF); // no eta extrap
    if (fParticles[i0].alphaC > -98 && fabs(fParticles[i0].Eta) < 2.5) v_alphaCs.push_back(fParticles[i0].alphaC); // no eta extrap
    // }
  }
  std::sort(v_alphaCs.begin(),v_alphaCs.end());
  std::sort(v_alphaFs.begin(),v_alphaFs.end());

  // std::cout << v_alphaCs.size() << ", " << v_alphaFs.size() << std::endl;

  if (v_alphaCs.size() > 0){
    int medIndex = v_alphaCs.size()/2+1;
    alphaCMed = v_alphaCs[medIndex];
    double sum = 0.0;
    for(unsigned int i = 0; i < v_alphaCs.size(); i++)
      sum += v_alphaCs[i] * v_alphaCs[i];
    alphaCRms = sqrt(sum / float(v_alphaCs.size()));
  }
  else{
    alphaCMed = 8.;
    alphaCRms = 8.;
  }

  if (v_alphaFs.size() > 0){
    int medIndex = v_alphaFs.size()/2+1;
    alphaFMed = v_alphaFs[medIndex];
    double sum = 0.0;
    for(unsigned int i = 0; i < v_alphaFs.size(); i++)
      sum += v_alphaFs[i] * v_alphaFs[i];
    alphaFRms = sqrt(sum / float(v_alphaFs.size()));
  }
  else{
    alphaFMed = 6.;
    alphaFRms = 6.;
  }

  // fudging this for now
  alphaCRms /= 2;
  alphaFRms /= 2;

}

void combiner::computeWeights(){

  // std::cout << "computing weights..." << std::endl;

  for(unsigned int i0   = 0; i0 < fParticles.size(); i0++) { 

    float curmed = alphaFMed;
    float currms = alphaFRms;
    float curalpha = fParticles[i0].alphaF;
    if (fabs(fParticles[i0].Eta) < 2.5){ // central!
      curmed = alphaCMed;
      currms = alphaCRms;
      curalpha = fParticles[i0].alphaC;
    }
    //PH: Getting good performance
    //if (fParticles[i0].isPV == 1 && fParticles[i0].charge != 0){ fParticles[i0].puppiWeight = 1; continue;}
    //if (fParticles[i0].isPV == 0 && fParticles[i0].charge != 0){ fParticles[i0].puppiWeight = 0; continue;}
    if (curalpha < -98){ fParticles[i0].puppiWeight = 0; continue; }

    float lVal = (curalpha - curmed)*fabs((curalpha - curmed))/currms/currms;
    float lPval = ROOT::Math::chisquared_cdf(lVal,1);
    // std::cout << "curmed = " << curmed << ", currms = " << currms << ", curalpha = " << curalpha << "..." << lVal << ", " << lPval << std::endl;
    
    fParticles[i0].puppiWeight = lPval;
    // if (curalpha - curmed) fParticles[i0].puppiWeight = 1;
  }

  // std::cout << "write out PUPPI collection..." << std::endl;
  int puppictr = 0;
  for(unsigned int i0   = 0; i0 < fParticles.size(); i0++) { 

    float ptcutC = 4.0;//Tight cuts for high PU
    float ptcutF = 4.0;
    if (fParticles[i0].puppiWeight <= 0.01) continue;
    if (fParticles[i0].puppiWeight == 1) { fParticlesPuppi.push_back(fParticles[i0]); puppictr++; }
    if (fParticles[i0].puppiWeight < 1 && fParticles[i0].puppiWeight > 0.01){
      if (fabs(fParticles[i0].Eta < 2.5)){
        if (fParticles[i0].Et*fParticles[i0].puppiWeight > ptcutC){
          fParticlesPuppi.push_back(fParticles[i0]);
          fParticlesPuppi[puppictr].Et = fParticles[i0].Et*fParticles[i0].puppiWeight;          
          puppictr++;
        }
      }
      else{
        if (fParticles[i0].Et*fParticles[i0].puppiWeight > ptcutF){
          fParticlesPuppi.push_back(fParticles[i0]);
          fParticlesPuppi[puppictr].Et = fParticles[i0].Et*fParticles[i0].puppiWeight;
          puppictr++;
        }
      }
    }
  }

  // std::cout << "puppi particles = " << fParticlesPuppi.size() << std::endl;

}

void combiner::fill(){

  fFile->cd();

  // std::cout << "filling..." << std::endl;
  b_pt.clear();
  b_alphaC.clear();
  b_alphaF.clear();
  b_Phi.clear();
  b_Eta.clear();
  b_Et.clear();
  b_PuppiWeight.clear();

  b_nParticles = fParticles.size();
  b_alphaFMed = alphaFMed;
  b_alphaFRms = alphaFRms;
  b_alphaCMed = alphaCMed;
  b_alphaCRms = alphaCRms;

  for(unsigned int i0   = 0; i0 < fParticles.size(); i0++) { 

    b_pt.push_back(fParticles[i0].Et);
    b_alphaC.push_back(fParticles[i0].alphaC);
    b_alphaF.push_back(fParticles[i0].alphaF);
    b_Eta.push_back(fParticles[i0].Eta);
    b_Phi.push_back(fParticles[i0].Phi);
    b_Et.push_back(fParticles[i0].Et);
    b_PuppiWeight.push_back(fParticles[i0].puppiWeight);

    // std::cout << "alphaF = " << fParticles[i0].alphaF << ", puppiWeight = " << fParticles[i0].puppiWeight << std::endl;
  }

  fTree->Fill();
}
