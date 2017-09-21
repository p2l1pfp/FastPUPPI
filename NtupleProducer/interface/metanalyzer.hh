#ifndef FASTPUPPI_NTUPLERPRODUCER_METANALYZER_HH
#define FASTPUPPI_NTUPLERPRODUCER_METANALYZER_HH
#include "L1TPFParticle.h"
#include "TFile.h"
#include "TTree.h"
#include <algorithm>  
#include <iostream>
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class metanalyzer { 
public:
  metanalyzer(std::string iFile);
  void setZ(std::vector<l1tpf::Particle> &iParticle,double iDZ);
  void setMETRecoil(int iId,std::vector<l1tpf::Particle> &iParticle,bool iAdd);
  void setGenMET(const reco::GenParticleCollection &iGenParticles);
  inline void clear() {for(int i0 = 0; i0 < fNVars; i0++) fVar[i0] = 0;}
  inline void fill()  {fFile->cd(); fTree->Fill();}
  inline void write() {fFile->cd(); fTree->Write(); fFile->Write();}   
private:
  TFile  *fFile;
  TTree  *fTree;
  int     fNVars;
  double  fVar[40];
};
#endif
