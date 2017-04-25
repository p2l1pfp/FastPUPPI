#ifndef FASTPUPPI_NTUPLERPRODUCER_JETANALYZER_HH
#define FASTPUPPI_NTUPLERPRODUCER_JETANALYZER_HH
#include "combiner.hh"
#include "TFile.h"
#include "TTree.h"
#include <algorithm>  
#include <iostream>
#include <fastjet/PseudoJet.hh>
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class jetanalyzer { 
public:
  jetanalyzer(std::string iFile, int debug=0);
  void setZ(std::vector<combiner::Particle> &iParticle,double iDZ);
  void setJets(std::vector<combiner::Particle> &iParticles,int iIndex);
  void setGenJets(const reco::GenParticleCollection &iGenParticles,int iIndex);
  std::vector<fastjet::PseudoJet> cluster(std::vector < fastjet::PseudoJet > &particles, double iRadius,double iPt);
  double dz(fastjet::PseudoJet &iJet,std::vector<combiner::Particle> &iParticles);
  double genmatch(int iId, fastjet::PseudoJet &matchjet,std::vector < fastjet::PseudoJet > &genjets);
  inline void clear() {for(int i0 = 0; i0 < fNVars; i0++) fVar[i0] = 0;}
  inline void fill()  {fFile->cd(); fTree->Fill();}
  inline void write() {fFile->cd(); fTree->Write(); fFile->Write();}   
private:
  TFile  *fFile;
  TTree  *fTree;
  int     fBase;
  int     fSize;
  int     fNVars;
  double  *fVar;
  std::vector<fastjet::PseudoJet> fGenJets;
  int     fDebug;
};
#endif
