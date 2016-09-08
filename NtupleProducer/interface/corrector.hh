#ifndef FASTPUPPI_NTUPLERPRODUCER_CORRECTOR_HH
#define FASTPUPPI_NTUPLERPRODUCER_CORRECTOR_HH
#include "TGraph.h"

class corrector { 
public:
  corrector(const std::string iFile,int iNFrac=11);
  double correct(double iHcal,double iEcal,int iEta);
  double ecalFrac();
  double hcalFrac();
private:
  TGraph ***fGraph;
  int fNEta;
  int fNFrac;
  int fIEta;
  double fFrac;
  double fEcal;
  double fHcal;
  double fPtCorr;
};
#endif
