#include "../interface/corrector.hh"
#include <iostream>
#include <sstream>
#include <TFile.h>

corrector::corrector(const std::string iFile,int iNFrac) {
  TFile *lFile = new TFile(iFile.c_str());
  fNEta  = 61;
  fNFrac = iNFrac;
  fGraph = new TGraph**[fNEta];
  for(int i0 = 0; i0 < fNEta; i0++) { 
    fGraph[i0] = new TGraph*[fNFrac];
    for(int i1 = 0; i1 < fNFrac; i1++) { 
      std::stringstream pSS; pSS << "eta_" << i0 << "_frac_" << i1;
      fGraph[i0][i1] = (TGraph*) lFile->Get(pSS.str().c_str());
    }
  }
}
double corrector::correct(double iHcal,double iEcal,int iEta) { 
  if(fNFrac == 1  && iEcal       < 1. ) return 0;
  if(fNFrac == 11 && iHcal+iEcal < 10.) return 0; 
  fEcal = 0; if(iEcal > 0) fEcal = iEcal; 
  fHcal = 0; if(iHcal > 0) fHcal = iHcal; 
  double lFrac = fEcal/(fHcal+fEcal); 
  int    lIFrac=int(lFrac*10.);
  fIEta = iEta+30;
  fFrac = lFrac;
  if(lIFrac > fNFrac-1) lIFrac = 0; 
  fPtCorr = fGraph[iEta+30][lIFrac]->Eval(fHcal+fEcal);
  return fPtCorr;
}
double corrector::ecalFrac() { 
  if(fEcal == 0) return 0;
  double lEcalCorr = fGraph[fIEta][fNFrac-1]->Eval(fEcal);
  return lEcalCorr/fPtCorr;
}
double corrector::hcalFrac() { 
  double lEcalCorr = fGraph[fIEta][10]->Eval(fEcal);
  return (1.-lEcalCorr/fPtCorr);
}
