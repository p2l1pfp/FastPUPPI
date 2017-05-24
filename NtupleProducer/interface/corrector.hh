#ifndef FASTPUPPI_NTUPLERPRODUCER_CORRECTOR_HH
#define FASTPUPPI_NTUPLERPRODUCER_CORRECTOR_HH
#include "TGraph.h"

class corrector { 
public:
  corrector(const std::string iFile,int iNFrac=11,int debug=0);
  double correct(double iHcal,double iEcal,int iEta,int iPhi);
private:
  TGraph ****fGraph;
  const int fNEta;
  const int fNPhi;
  const int fNFrac;
};
#endif
