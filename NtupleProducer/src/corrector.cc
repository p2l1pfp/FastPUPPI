#include "../interface/corrector.hh"
#include <iostream>
#include <sstream>
#include <TFile.h>

corrector::corrector(const std::string iFile,int iNFrac) {
  TFile *lFile = new TFile(iFile.c_str());
  fNEta  = 61;
  fNPhi  = 72;
  fNFrac = iNFrac;
  fGraph = new TGraph***[fNEta];
  lFile->Print();
  for(int i0 = 0; i0 < fNEta; i0++) { 
    fGraph[i0] = new TGraph**[fNPhi];
    for(int i1 = 0; i1 < fNPhi; i1++) { 
      fGraph[i0][i1] = new TGraph*[fNFrac];
      for(int i2 = 0; i2 < fNFrac; i2++) { 
	std::stringstream pSS; 
	pSS << "eta_";
	int pEta = i0;
	if(iNFrac > 1) pEta = i0-30;
	pSS << pEta;
	if(iNFrac > 1) pSS << "phi_" << (i1+1);
	pSS << "_frac_" << i2;
	fGraph[i0][i1][i2] = (TGraph*) lFile->FindObjectAny(pSS.str().c_str());
      }
    }
  }
}
double corrector::correct(double iHcal,double iEcal,int iEta,int iPhi) { 
  if(fNFrac == 1  && iEcal       < 1. ) return 0;
  if(fNFrac == 11 && iEcal+iHcal < 2. ) return 0; 
  fEcal = 0; if(iEcal > 0) fEcal = iEcal; 
  fHcal = 0; if(iHcal > 0) fHcal = iHcal; 
  double lFrac = fEcal/(fHcal+fEcal); 
  int    lIFrac=int(lFrac*10.);
  fIEta = iEta+30;
  fIPhi = iPhi-1;
  fFrac = lFrac;
  if(lIFrac > fNFrac-1) lIFrac = 0; 
  fPtCorr = (fGraph[fIEta][fIPhi][lIFrac])->Eval(fHcal+fEcal);
  if(fPtCorr < 1.) fPtCorr = 0;
  return fPtCorr;
}
double corrector::ecalFrac() { 
  if(fEcal == 0) return 0;
  double lEcalCorr = fGraph[fIEta][fIPhi][fNFrac-1]->Eval(fEcal);
  return lEcalCorr/fPtCorr;
}
double corrector::hcalFrac() { 
  double lEcalCorr = fGraph[fIEta][fIPhi][10]->Eval(fEcal);
  return (1.-lEcalCorr/fPtCorr);
}
