#include "../interface/isoanalyzer.hh"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/JetDefinition.hh>

isoanalyzer::isoanalyzer(std::string iFile, int debug) {

  fFile = new TFile(iFile.c_str(),"RECREATE");
  fTree = new TTree("iso","iso");

  fTree->Branch("gmu1_pt"    ,&_gmu1_pt   ,"gmu1_pt/F");
  fTree->Branch("gmu1_eta"   ,&_gmu1_eta  ,"gmu1_eta/F");
  fTree->Branch("gmu1_phi"   ,&_gmu1_phi  ,"gmu1_phi/F");
  fTree->Branch("gmu2_pt"    ,&_gmu2_pt   ,"gmu2_pt/F");
  fTree->Branch("gmu2_eta"   ,&_gmu2_eta  ,"gmu2_eta/F");
  fTree->Branch("gmu2_phi"   ,&_gmu2_phi  ,"gmu2_phi/F");

  fTree->Branch("mu1_pt"    ,&_mu1_pt   ,"mu1_pt/F");
  fTree->Branch("mu1_eta"   ,&_mu1_eta  ,"mu1_eta/F");
  fTree->Branch("mu1_phi"   ,&_mu1_phi  ,"mu1_phi/F");
  fTree->Branch("mu1_iso_pf"      ,&_mu1_iso_pf  ,"mu1_iso_pf/F");
  fTree->Branch("mu1_iso_tk"      ,&_mu1_iso_tk  ,"mu1_iso_tk/F");
  fTree->Branch("mu1_iso_tkvtx"   ,&_mu1_iso_tkvtx  ,"mu1_iso_tkvtx/F");
  fTree->Branch("mu1_iso_pup"     ,&_mu1_iso_pup  ,"mu1_iso_pup/F");

  fTree->Branch("mu2_pt"    ,&_mu2_pt   ,"mu2_pt/F");
  fTree->Branch("mu2_eta"   ,&_mu2_eta  ,"mu2_eta/F");
  fTree->Branch("mu2_phi"   ,&_mu2_phi  ,"mu2_phi/F");
  fTree->Branch("mu2_iso_pf"      ,&_mu2_iso_pf     ,"mu2_iso_pf/F");
  fTree->Branch("mu2_iso_tk"      ,&_mu2_iso_tk     ,"mu2_iso_tk/F");
  fTree->Branch("mu2_iso_tkvtx"   ,&_mu2_iso_tkvtx  ,"mu2_iso_tkvtx/F");
  fTree->Branch("mu2_iso_pup"     ,&_mu2_iso_pup    ,"mu2_iso_pup/F");

  // debug
  fDebug = debug;
}

void isoanalyzer::setGenMuons(const reco::GenParticleCollection &iGenParticles,int iIndex) { 
  
  //save the top 2 muons

  for (reco::GenParticleCollection::const_iterator itGenP = iGenParticles.begin(); itGenP!=iGenParticles.end(); ++itGenP) {
    if(itGenP->status() == 23 && abs(itGenP->pdgId()) == 13){ // pythia codes??? use status 23 for now
      // std::cout << "muon found!" << itGenP->status() << " " <<  itGenP->pt() << std::endl;
      float curpt = itGenP->pt();
      if (curpt > _gmu1_pt){
        _gmu1_pt = curpt;
        _gmu1_eta = itGenP->eta();
        _gmu1_phi = itGenP->phi();
      }
      else if (curpt > _gmu2_pt){
        _gmu2_pt = curpt;
        _gmu2_eta = itGenP->eta();
        _gmu2_phi = itGenP->phi();
      }
      else continue;
    }
  }

  // std::cout << "gen muon info = " << _gmu1_pt << " " << _gmu1_eta << " " << _gmu1_phi << ", " << _gmu2_pt << " " << _gmu2_eta << " " << _gmu2_phi << std::endl;

}

void isoanalyzer::matchMuons(std::vector<combiner::Particle> &iParticle){
  
  TLorentzVector gmu1;
  TLorentzVector gmu2;
  gmu1.SetPtEtaPhiM(_gmu1_pt,_gmu1_eta,_gmu1_phi,0.109);
  gmu2.SetPtEtaPhiM(_gmu2_pt,_gmu2_eta,_gmu2_phi,0.109);

  _mu1.SetPtEtaPhiM(0.,0.,0.,0.);
  _mu2.SetPtEtaPhiM(0.,0.,0.,0.);

  // std::cout << "computing isolation..." << gmu1.Pt() << "," << gmu2.Pt() << std::endl; 

  // first find the muon matched to GEN
  for(unsigned   int i0 = 0; i0 < iParticle.size(); i0++) {
    if ((gmu1.Pt() > 0)&&(iParticle[i0].pdgId() == combiner::MU)){
      if (gmu1.DeltaR(iParticle[i0].tp4()) < 0.4){
        // std::cout << "found a matching muon1!" << iParticle[i0].tp4().Pt() << std::endl;
        _mu1 = iParticle[i0].tp4();
      }
    }
    if ((gmu2.Pt() > 0)&&(iParticle[i0].pdgId() == combiner::MU)){
      if (gmu2.DeltaR(iParticle[i0].tp4()) < 0.4){
        // std::cout << "found a matching muon2!" << iParticle[i0].tp4().Pt() << std::endl;
        _mu2 = iParticle[i0].tp4();        
      }
    }    
  }

  // setting info
  _mu1_pt  = _mu1.Pt();
  _mu1_eta = _mu1.Eta();
  _mu1_phi = _mu1.Phi();
  _mu2_pt  = _mu2.Pt();
  _mu2_eta = _mu2.Eta();
  _mu2_phi = _mu2.Phi();

}

void isoanalyzer::computeIso(std::vector<combiner::Particle> &iParticle, float coneSize, std::string type){

  TLorentzVector mu1 = _mu1;
  TLorentzVector mu2 = _mu2;

  float iso1 = 0;
  float iso2 = 0;

  // compute mu1 isolation
  if (mu1.Pt() > 0){
    for(unsigned   int i0 = 0; i0 < iParticle.size(); i0++) {
      float curdr = mu1.DeltaR(iParticle[i0].tp4());
      if ((curdr < coneSize) && (curdr > 0.001)){ // protection from the muon itself
        iso1 += iParticle[i0].tp4().Pt();
      }
    }
  }
  if (mu2.Pt() > 0){
    for(unsigned   int i0 = 0; i0 < iParticle.size(); i0++) {
      float curdr = mu2.DeltaR(iParticle[i0].tp4());
      if ((curdr < coneSize) && (curdr > 0.001)){ // protection from the muon itself
        iso2 += iParticle[i0].tp4().Pt();
      }
    }
  }

  if (mu1.Pt() > 0) iso1 /= mu1.Pt();
  if (mu2.Pt() > 0) iso2 /= mu2.Pt();
  // std::cout << type << ", isolations = " << iso1 << ", " << iso2 << std::endl;

  if (type == "pf"){
    _mu1_iso_pf = iso1; _mu2_iso_pf = iso2;
  }
  else if (type == "tk"){
    _mu1_iso_tk = iso1; _mu2_iso_tk = iso2;
  }
  else if (type == "tkvtx"){
    _mu1_iso_tkvtx = iso1; _mu2_iso_tkvtx = iso2;
  }
  else if (type == "pup"){
    _mu1_iso_pup = iso1; _mu2_iso_pup = iso2;
  }

}

// double isoanalyzer::genmatch(int iId, fastjet::PseudoJet &matchjet,std::vector < fastjet::PseudoJet > &genjets){
//   double lDR = 1000;
//   double lPt = -1;
//   for(unsigned int i0 = 0; i0 < genjets.size(); i0++) { 
//     double pDR = matchjet.delta_R(genjets[i0]);
//     if(pDR > 0.25) continue;
//     if(lDR < pDR)  continue;
//     lDR = pDR;
//     lPt = genjets[i0].pt();
//   }
//   if(iId == 0) return lDR;
//   return lPt;
// }

void isoanalyzer::clear(){

    _gmu1_pt  = -999.;
    _gmu1_eta = -999.;
    _gmu1_phi = -999.;
    _gmu2_pt  = -999.;
    _gmu2_eta = -999.;
    _gmu2_phi = -999.;

    _mu1_pt  = -999.;
    _mu1_eta = -999.;
    _mu1_phi = -999.;
    _mu1_iso_pf = -999.;
    _mu1_iso_tk = -999.;
    _mu1_iso_tkvtx = -999.;
    _mu1_iso_pup = -999.;    
    
    _mu2_pt  = -999.;
    _mu2_eta = -999.;
    _mu2_phi = -999.;
    _mu2_iso_pf = -999.;
    _mu2_iso_tk = -999.;
    _mu2_iso_tkvtx = -999.;
    _mu2_iso_pup = -999.;    


}





