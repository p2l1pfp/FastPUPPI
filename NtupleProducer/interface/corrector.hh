#ifndef FASTPUPPI_NTUPLERPRODUCER_CORRECTOR_HH
#define FASTPUPPI_NTUPLERPRODUCER_CORRECTOR_HH
#include "TGraph.h"

class corrector { 
public:
  corrector(const std::string iFile,int iNFrac=11);
  double correct(double iHcal,double iEcal,int iEta,int iPhi);
  double correct(double iCorr,double iHcal,double iEcal,int iEta,int iPhi);
  double ecalFrac();
  double hcalFrac();
private:
  TGraph ****fGraph;
  int fNEta;
  int fNPhi;
  int fNFrac;
  int fIEta;
  int fIPhi;
  double fFrac;
  double fEcal;
  double fHcal;
  double fPtCorr;
};
#endif
