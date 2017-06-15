#ifndef FASTPUPPI_NTUPLERPRODUCER_ISOANALYZER_HH
#define FASTPUPPI_NTUPLERPRODUCER_ISOANALYZER_HH
#include "L1TPFParticle.h"
#include "TFile.h"
#include "TTree.h"
#include <algorithm>  
#include <iostream>
#include <fastjet/PseudoJet.hh>
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class isoanalyzer { 
public:
  isoanalyzer(std::string iFile, int debug=0);
  void setGenMuons(const reco::GenParticleCollection &iGenParticles,int iIndex);
  void matchMuons(std::vector<l1tpf::Particle> &iParticle);
  void computeIso(std::vector<l1tpf::Particle> &iParticle, float coneSize, std::string type);
  // double genmatch(int iId, fastjet::PseudoJet &matchjet,std::vector < fastjet::PseudoJet > &genjets);
  void clear();
  inline void fill()  {fFile->cd(); fTree->Fill();}
  inline void write() {fFile->cd(); fTree->Write(); fFile->Write();}   
private:
  TFile  *fFile;
  TTree  *fTree;
  int     fBase;
  int     fSize;
  int     fNVars;

  float     _gmu1_pt;
  float     _gmu1_eta;
  float     _gmu1_phi;
  float     _gmu2_pt;
  float     _gmu2_eta;
  float     _gmu2_phi;

  TLorentzVector _mu1;
  TLorentzVector _mu2;

  float     _mu1_pt;
  float     _mu1_eta;
  float     _mu1_phi;
  float     _mu1_iso_pf;
  float     _mu1_iso_tk;
  float     _mu1_iso_tkvtx;
  float     _mu1_iso_pup;

  float     _mu2_pt;
  float     _mu2_eta;
  float     _mu2_phi;
  float     _mu2_iso_pf;
  float     _mu2_iso_tk;
  float     _mu2_iso_tkvtx;
  float     _mu2_iso_pup;


  int     fDebug;
};
#endif
