#ifndef FASTPUPPI_NTUPLERPRODUCER_METANALYZER_HH
#define FASTPUPPI_NTUPLERPRODUCER_METANALYZER_HH
#include "combiner.hh"
#include "TFile.h"
#include "TTree.h"
#include <algorithm>  
#include <iostream>

class metanalyzer { 
public:
  metanalyzer(std::string iFile);
  void setZ(std::vector<combiner::Particle> &iParticle);
  void setMETRecoil(int iId,std::vector<combiner::Particle> &iParticle,bool iAdd);
  inline void clear() {for(int i0 = 0; i0 < fNVars; i0++) fVar[i0] = 0;}
  inline void fill()  {fFile->cd(); fTree->Fill();}
  inline void write() {fFile->cd(); fTree->Write(); fFile->Write();}   
private:
  TFile  *fFile;
  TTree  *fTree;
  int     fNVars;
  double  fVar[20];
};
#endif
