#include "../interface/combiner.hh"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iostream>
#include <sstream>
#include <TFile.h>
#include <TH1F.h>
#include <algorithm>
#include <cassert>
#include "TMath.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

combiner::combiner(const std::string & iPionFile,const std::string & iElectronFile,const std::string & iTrackFile,const std::string & iFile,double iEtaCharged,double iPuppiPt,double iPuppiDr,double iVtxRes,int debug) {
  fDebug = debug; // 0 = nothing; 1 = something; 2 = lots; ...
  fEta     = iEtaCharged;
  fPuppiPt = iPuppiPt;
  fPuppiDr = iPuppiDr;
  fDRMatch = 0.15; fPtMatchLow = 2.0; fPtMatchHigh = 2.0; fUseTrackCaloSigma = false; fRescaleUnmatchedTrack = false; fMaxInvisiblePt = 20.0;
  fVtxRes  = iVtxRes;
  fNEta  = l1tpf::towerNEta();
  loadFile(fPionRes    ,iPionFile);
  loadFile(fElectronRes,iElectronFile);
  loadFile(fTrackRes   ,iTrackFile);

  if (iFile != "") {
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
  } else {
      fFile = nullptr;
      fTree = nullptr;
  }
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
l1tpf::Particle combiner::makeCalo(double iCalo,double iEcal,double iCaloEta,double iCaloPhi,double iEcalEta,double iEcalPhi) const { 
  double lEt  = iCalo;
  double lEta = iCaloEta;
  double lPhi = iCaloPhi;
  if(iCalo < 1.1*iEcal) lEt  = iEcal;
  if(iCalo < 1.1*iEcal) lEta = iEcalEta;
  if(iCalo < 1.1*iEcal) lPhi = iEcalPhi;
  double lSigma = lEt;
  iCalo < 1.1*iEcal ? lSigma = getEleRes(lEt,lEta,lPhi) : lSigma = getPionRes(lEt,lEta,lPhi);
  int lId = Particle::NH; if(iCalo < 1.1*iEcal) lId = Particle::GAMMA;
  return Particle(lEt,lEta,lPhi,0.,lId,lSigma,0.,lEta,lPhi);
}
void combiner::addCalo(const l1tpf::Particle & particle) { 
  if(particle.pt() < 0.01) return;
  assert(particle.charge() == 0);
  assert(particle.pdgId() == 2 || particle.pdgId() == 3);
  insert(particle,fParticles);
}

//Add Tracks
void combiner::addTrack(const l1tpf::Particle &particle) { 
  if(particle.pt() < 0.1) return; //Just to avoid calling things that will crash later
  assert(particle.charge() != 0);
  assert(particle.pdgId() == 0);
  insert(particle,fTkParticles);
}

//Add muons
void combiner::addMuon(const l1tpf::Particle &particle) {
  assert(particle.charge() != 0);
  assert(particle.pdgId() == 4);
  fMuParticles.push_back(particle);
}

//Iterate down in Track Pt
void combiner::link(bool iMetRate) { 
  // first do the muon/tracking matching
  for(unsigned int i0   = 0; i0 < fMuParticles.size(); i0++) { 
    float minPtDiff = 9999; // find the track that best matches the muon in pT
    int imatch = -1;
    for(unsigned int i1 = 0; i1 < fTkParticles.size(); i1++) { 
      float curDR = deltaRraw(fTkParticles[i1],fMuParticles[i0]);
      float curDPT = (fMuParticles[i0].pt()/fTkParticles[i1].pt());
      if (curDPT < 1){ // invert it
        float tmptmp = curDPT;
        curDPT = 1/tmptmp;
      }
      if (curDR < 0.2 && (curDPT) < 4 && fabs(curDPT) < minPtDiff){
        if (imatch > -1) { // we've found a new match, the previous one should be unmatched
            fTkParticles[i1].setPdgId(Particle::CH);
            fTkParticles[i1].setMass(0.135);
        }
        fTkParticles[i1].setPdgId(Particle::MU);
        fTkParticles[i1].setMass(0.105);
        minPtDiff = fabs(curDPT);
        imatch = i1;
      }
    } 
  }
  // then do track + calo linking
  if (fDebug) printf("Trying to link. I have %d tracks, %d calo\n", int(fTkParticles.size()), int(fParticles.size()));
  for(unsigned int i0   = 0; i0 < fTkParticles.size(); i0++) { 
    if(fTkParticles[i0].pdgId() == Particle::MU) continue; // skip muons for now, add them at the end
    if (fDebug>1) printf("\t track %d (pt %7.2f +- %5.2f, est. calo unc +- %5.2f)\n", i0, fTkParticles[i0].pt(), fTkParticles[i0].sigma(), fTkParticles[i0].caloSigma());
    bool pFilled = false; // if true, the track has already been added in the PF candidates; if false, it's still to be added
    int pIMatch = -1; double pPtMatch = -1; int pIMatchPtTooLow = -1;
    for(unsigned int i1 = 0; i1 < fParticles.size();   i1++) { 
      if (fDebug>1) printf("\t\t calo %d (pt %7.2f +- %5.2f): ", i1, fParticles[i1].pt(), fParticles[i1].sigma());
      // what happens if there is no matching cluster?  does it throw out the track? it should to reduce fake tracks (see below=> still debating)
      float caloSigma = fUseTrackCaloSigma ? fTkParticles[i0].caloSigma() : fParticles[i1].sigma();
      if(combiner::deltaR(fTkParticles[i0],fParticles[i1]) > fDRMatch) { 
            if (fDebug>1) printf("outside dR (%.3f).\n",combiner::deltaR(fTkParticles[i0],fParticles[i1])); 
            continue; }
      if(fParticles[i1].pdgId() == Particle::CH || fParticles[i1].pdgId() == Particle::EL)  { 
            if (fDebug>1) printf("already linked.\n");  
            continue; }
      if(fParticles[i1].pt()+fPtMatchLow*caloSigma < fTkParticles[i0].pt()) { 
            if (pIMatchPtTooLow == -1) pIMatchPtTooLow = i1; // calo are pt sorted, so this is always the highest pt unlinked
            if (fDebug>1) printf("calo pt is too low.\n"); 
            continue; }
      if(pIMatch != -1 && fabs(fParticles[i1].pt()-fTkParticles[i0].pt()) > pPtMatch) { 
            if (fDebug>1) printf("new match, but worse than %d.\n", pIMatch); 
            continue; }
      pIMatch = i1; pPtMatch = fabs(fParticles[i1].pt()-fTkParticles[i0].pt());
      if (fDebug>1) printf("new best match.\n");
    }
    if(pIMatch != -1) { 
      if (fDebug) printf("Linking track of pt %7.2f eta %+5.2f phi %+5.2f (calo eta %+5.2f phi %+5.2f) to calo of pt %7.2f eta %+5.2f phi %+5.2f (dR = %.3f, deta = %.3f, dphi = %.3f)\n",
                    fTkParticles[i0].pt(), fTkParticles[i0].eta(), fTkParticles[i0].phi(), fTkParticles[i0].caloEta(), fTkParticles[i0].caloPhi(),
                    fParticles[pIMatch].pt(), fParticles[pIMatch].eta(), fParticles[pIMatch].phi(),
                    combiner::deltaR(fTkParticles[i0],fParticles[pIMatch]), (fTkParticles[i0].caloEta()-fParticles[pIMatch].eta()), ::deltaPhi(fTkParticles[i0].caloPhi(),fTkParticles[i0].phi()));
      merge(fTkParticles[i0],fParticles[pIMatch],fParticles);
      pFilled = true;
    }
    //Remove high pT fakes
    if(fTkParticles[i0].pt() > fMaxInvisiblePt) {
        if (fRescaleUnmatchedTrack && pIMatchPtTooLow != -1) {
            Particle & iParticle1 = fParticles[pIMatchPtTooLow];
            iParticle1.setPtEtaPhiM(iParticle1.pt(), fTkParticles[i0].eta(), fTkParticles[i0].phi(), 0.137);
            iParticle1.setPdgId(iParticle1.pdgId() == Particle::GAMMA ? Particle::EL : Particle::CH); 
            iParticle1.setDz(fTkParticles[i0].dz());
        }
        pFilled = true; // mark it as used so it's not turned into a PF candidate
    }
    if(!pFilled) fParticles.push_back(fTkParticles[i0]);
  }
  // now do muons... when not using muon gun+PU as a neutrino gun
  if(!iMetRate) { 
    for(unsigned int i0   = 0; i0 < fTkParticles.size(); i0++) { 
      if (fTkParticles[i0].pdgId() == Particle::MU) fParticles.push_back(fTkParticles[i0]);
    }
  }
}

//Merge Particles
void combiner::merge(Particle &iTkParticle,Particle &iParticle1,std::vector<Particle> &iCollection) { 

  float caloSigma = fUseTrackCaloSigma ? iTkParticle.caloSigma() : iParticle1.sigma();
  float lTotSigma = hypot(iTkParticle.sigma(), caloSigma);
  //Case 1 calo to large
  if(iParticle1.pt()-iTkParticle.pt() > fPtMatchHigh*lTotSigma) {
    //TLorentzVector pVec0 = iParticle1.tp4();
    //TLorentzVector pVec1 = iTkParticle.tp4();
    //pVec0 -= pVec1;
    //iParticle1.setPtEtaPhiM(pVec0.Pt(),pVec0.Eta(),pVec0.Phi(),0);
    //- scalar subtraction works slightly better than vector one, in addition to be simpler
    iParticle1.setPt(iParticle1.pt()-iTkParticle.pt());
    iCollection.push_back(iTkParticle);
    if (fDebug) printf("   case 1: add track, reduce calo pt to %7.2f\n", iParticle1.pt());
  } else {
    double pSigma   = caloSigma; //iParticle1.sigma();
    double pTkSigma = iTkParticle.sigma();
    double pSigTot  = 1./pSigma/pSigma + 1./pTkSigma/pTkSigma;
    double pAvgEt   = (iParticle1.pt()/pSigma/pSigma + iTkParticle.pt()/pTkSigma/pTkSigma)/pSigTot;
    iParticle1.setPtEtaPhiM(pAvgEt, iTkParticle.eta(), iTkParticle.phi(), 0.137);
    iParticle1.setPdgId(iParticle1.pdgId() == Particle::GAMMA ? Particle::EL : Particle::CH); 
    iParticle1.setDz(iTkParticle.dz());
    iParticle1.setCharge(iTkParticle.charge());
    if (fDebug) printf("   case 2: merge, avg pt %7.2f\n", iParticle1.pt());
  }
}

//Order by pt // FIXME: GP: why do you need to do this? => PH: the ordering affects the choice of what gets linked first if there are multiple links present. This changes what happens downstream.
// Perhaps there is a better way :)
void combiner::insert(const Particle &iParticle,std::vector<Particle> &iParticles) { 
  if(iParticles.size() == 0) {iParticles.push_back(iParticle); return;}
  bool lFill = false;
  for(std::vector<Particle>::iterator pParticle = iParticles.begin(); pParticle != iParticles.end(); pParticle++) { 
    if(pParticle->pt() > iParticle.pt()) continue; 
    iParticles.insert(pParticle,iParticle);
    lFill = true;
    break;
  }
  if(!lFill) iParticles.push_back(iParticle);
}

//Delta R
double  combiner::deltaR(Particle &iParticle1,Particle &iParticle2) {
  return ::deltaR(iParticle1.caloEta(),iParticle1.caloPhi(),iParticle2.caloEta(),iParticle2.caloPhi());
}

//Delta R
double  combiner::deltaRraw(Particle &iParticle1,Particle &iParticle2) {
  return ::deltaR(iParticle1.eta(),iParticle1.phi(),iParticle2.eta(),iParticle2.phi());
}

void combiner::doVertexing(){
  
  // std::vector<Particle> fTkParticlesWVertexing;
  int lNBins = int(40./fVtxRes);
  TH1F *h_dz = new TH1F("h_dz","h_dz",lNBins,-20,20); // 1cm binning
  for (int i = 0; i < h_dz->GetXaxis()->GetNbins(); ++i) h_dz->SetBinContent(i+1,0.); //initialize all to 0.
  int ntracks = 0;
  for(unsigned int i0   = 0; i0 < fParticles.size(); i0++) { 
    if (fParticles[i0].charge() == 0) continue;
    ntracks++;
    //if(fParticles[i0].quality() < 7) continue;
    float curdz  = fParticles[i0].dz();
    float curbin = h_dz->GetXaxis()->FindBin(curdz);
    h_dz->SetBinContent( curbin, h_dz->GetBinContent(curbin) + std::min(fParticles[i0].pt(),50.) );
   }
  int imaxbin = h_dz->GetMaximumBin();
  float pvdz = h_dz->GetXaxis()->GetBinCenter(imaxbin);
  float binwidth = h_dz->GetXaxis()->GetBinWidth(imaxbin);
  float pvdz_lo = pvdz - 1.5*binwidth;
  float pvdz_hi = pvdz + 1.5*binwidth;
  fDZ = pvdz;
  for(unsigned int i0   = 0; i0 < fParticles.size(); i0++) { 
    float curdz  = fParticles[i0].dz();
    if(fParticles[i0].charge() != 0) fParticles[i0].setIsPV(0);
    if (curdz < pvdz_hi && curdz > pvdz_lo){
      if(fParticles[i0].charge() != 0) fTkParticlesWVertexing.push_back( fParticles[i0] );
      if(fParticles[i0].charge() != 0) fParticles[i0].setIsPV(1);
   }
  }
  delete h_dz;
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
    if (fParticles[i0].charge() == 0) puppi_neutrals.push_back(fParticles[i0]);
    if (fParticles[i0].charge() != 0 && fParticles[i0].isPV() != 0 && fParticles[i0].pdgId() < 2) puppi_chargedPV.push_back(fParticles[i0]);
    if (fParticles[i0].charge() != 0 && fParticles[i0].isPV() == 0) puppi_chargedPU.push_back(fParticles[i0]);
  }
  if (fDebug) std::cout << fParticles.size() << ", " << puppi_neutrals.size() << ", " << puppi_chargedPU.size() << ", " << puppi_chargedPV.size() << std::endl;

  // compute alphas (not efficient, but we don't have that many particles)
  computeAlphas(fParticles,0);
  computeAlphas(puppi_chargedPV,1);
  //computeAlphas(fParticles,1);

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
      if (curdr < fPuppiDr && curdr != 0) curalpha += neighborParticles[j0].pt() * neighborParticles[j0].pt() / curdr / curdr;
    }
    float logcuralpha = -99;
    if (curalpha != 0.) logcuralpha = log(curalpha);
    if (isCentral == 0) fParticles[i0].setAlphaF(logcuralpha);
    if (isCentral == 1) fParticles[i0].setAlphaC(logcuralpha);
  }

}

void combiner::computeMedRMS(){

  // std::cout << "computing alpha med/rms..." << std::endl;
  std::vector<float> v_alphaFs;
  std::vector<float> v_alphaCs;
  for(unsigned int i0   = 0; i0 < fParticles.size(); i0++) { 
    // if (fParticles[i0].charge() != 0 && fParticles[i0].isPV == 0){
    // if (fParticles[i0].alphaF() > -98) v_alphaFs.push_back(fParticles[i0].alphaF());
    // if (fParticles[i0].alphaC() > -98) v_alphaCs.push_back(fParticles[i0].alphaC());

    //Use non-pv tracks for alphaF in central
    if (fabs(fParticles[i0].eta()) < fEta && (fParticles[i0].pdgId() == 2 || fParticles[i0].pdgId() == 3 || fParticles[i0].isPV() == 1)) continue;
    if (fParticles[i0].alphaF() > -98 && fabs(fParticles[i0].eta()) > fEta) v_alphaFs.push_back(fParticles[i0].alphaF()); // no eta extrap
    if (fParticles[i0].alphaC() > -98 && fabs(fParticles[i0].eta()) < fEta) v_alphaCs.push_back(fParticles[i0].alphaC()); // no eta extrap
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
      sum += (v_alphaCs[i]) * (v_alphaCs[i]);
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
      sum += (v_alphaFs[i]) * (v_alphaFs[i]);
    alphaFRms = sqrt(sum / float(v_alphaFs.size()));
  }
  else{
    alphaFMed = 6.;
    alphaFRms = 6.;
  }

  // fudging this for now
  alphaCRms /= 2;
  alphaFRms /= 2;

  // fudging this for now
  //alphaCRms *= 10.*sqrt(2.);
  //alphaFRms *= 10.*sqrt(2.);

}

void combiner::computeWeights(){

  // std::cout << "computing weights..." << std::endl;

  for(unsigned int i0   = 0; i0 < fParticles.size(); i0++) { 

    float curmed = alphaFMed;
    float currms = alphaFRms;
    float curalpha = fParticles[i0].alphaF();
    if (fabs(fParticles[i0].eta()) < fEta){ // central!
      curmed = alphaCMed;
      currms = alphaCRms;
      curalpha = fParticles[i0].alphaC();
    }
    //PH: Getting good performance
    if (fParticles[i0].charge() != 0){ 
        fParticles[i0].setPuppiWeight(fParticles[i0].isPV() ? 1.0 : 0.0); 
        continue;
    }
    if (curalpha < -98){ fParticles[i0].setPuppiWeight(0); continue; }

    float lVal = (curalpha - curmed)*fabs((curalpha - curmed))/currms/currms;
    float lPval = ROOT::Math::chisquared_cdf(lVal,1);
    fParticles[i0].setPuppiWeight(lPval);
  }

  // std::cout << "write out PUPPI collection..." << std::endl;
  float ptcutC = fPuppiPt;//Tight cuts for high PU
  float ptcutF = fPuppiPt*2.0;
  for (const l1tpf::Particle & particle : fParticles) {
      if (particle.pdgId() == 4) {
          fParticlesPuppi.push_back(particle);
      } else if (particle.charge() != 0) {
          assert(particle.puppiWeight() == particle.isPV());
          if (particle.isPV()) fParticlesPuppi.push_back(particle);
      } else {
          if (particle.puppiWeight() > 0.01) {
              float rescpt = particle.puppiWeight() * particle.pt();
              if (rescpt > (std::abs(particle.eta()) < fEta ? ptcutC : ptcutF)) {
                  fParticlesPuppi.push_back(particle); 
                  fParticlesPuppi.back().setPt(rescpt); 
              }
          }
      }
  }

}

void combiner::fill(){
  if (!fFile) return;

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

    b_pt.push_back(fParticles[i0].pt());
    b_alphaC.push_back(fParticles[i0].alphaC());
    b_alphaF.push_back(fParticles[i0].alphaF());
    b_Eta.push_back(fParticles[i0].eta());
    b_Phi.push_back(fParticles[i0].phi());
    b_Et.push_back(fParticles[i0].pt());
    b_PuppiWeight.push_back(fParticles[i0].puppiWeight());

    // std::cout << "alphaF = " << fParticles[i0].alphaF() << ", puppiWeight = " << fParticles[i0].puppiWeight() << std::endl;
  }

  fTree->Fill();
}
