#include "../interface/corrector.hh"
#include "../interface/L1TPFUtils.h"
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <TFile.h>
#include <TKey.h>
#include "FWCore/Utilities/interface/CPUTimer.h"
#include <unordered_map>

corrector::corrector(const std::string iFile, int iNFrac, int debug) {
  TFile *lFile = new TFile(iFile.c_str());
  fNEta  = l1tpf::towerNEta();
  fNPhi  = l1tpf::towerNPhi(1);
  fNFrac = iNFrac;
  fGraph = new TGraph***[fNEta];
  if (debug) lFile->Print();
  edm::CPUTimer timer;
  timer.start();
  char buff[1023];
  std::unordered_map<std::string,TGraph *> graphs;
  TKey *key;
  TIter nextkey(lFile->GetListOfKeys());
  while ((key = (TKey *) nextkey())) {
      if (strncmp(key->GetName(), "eta_", 4) == 0) {
          graphs[key->GetName()] = (TGraph*) key->ReadObj();
      }
  }
  unsigned int ngraphs = 0;
  for(int i0 = 0; i0 < fNEta; i0++) { 
    fGraph[i0] = new TGraph**[fNPhi];
    for(int i1 = 0; i1 < fNPhi; i1++) { 
      fGraph[i0][i1] = new TGraph*[fNFrac];
      for(int i2 = 0; i2 < fNFrac; i2++) { 
        if (iNFrac > 1) snprintf(buff, 1022, "eta_%dphi_%d_frac_%d", i0-fNEta/2, i1+1, i2);
        else            snprintf(buff, 1022, "eta_%d_frac_%d", i0, i2);
        fGraph[i0][i1][i2] = graphs[buff]; 
      }
    }
  }
  timer.stop();
  if (debug) std::cout << "Read " << ngraphs << " graphs from " << iFile << " in " << timer.realTime() << " s"  << std::endl; 
}
double corrector::correct(double iHcal,double iEcal,int iEta,int iPhi) { 
  if(fNFrac == 1  && iEcal       < 0.5 ) return 0;
  if(fNFrac == 11 && iEcal+iHcal < 1.0 ) return 0; 
  fEcal = 0; if(iEcal > 0) fEcal = iEcal; 
  fHcal = 0; if(iHcal > 0) fHcal = iHcal; 
  if(iEcal < 0) fEcal = 0;
  if(abs(iEta) > fNEta/2-1 || iPhi < 1 || iPhi > fNPhi-1) return fHcal+fEcal;
  double lFrac = fEcal/(fHcal+fEcal); 
  int    lIFrac=int(lFrac*10.);
  fIEta = iEta+fNEta/2;
  fIPhi = iPhi-1;
  fFrac = lFrac;
  if(lIFrac > fNFrac-1) lIFrac = 0; 
  fPtCorr = (fGraph[fIEta][fIPhi][lIFrac])->Eval(fHcal+fEcal);
  fPtCorr = std::min(10.*(fHcal+fEcal),fPtCorr);
  if(fPtCorr < 1.) fPtCorr = 0;
  if(fabs(l1tpf::towerEta(iEta)) > 3.0) fPtCorr *= 0.66;//
  if(fabs(l1tpf::towerEta(iEta)) > 2.853 && fabs(l1tpf::towerEta(iEta)) < 3.1 && iHcal > 4.) fPtCorr = 0;
  return fPtCorr;
}
//Double correction
double corrector::correct(double iCorr,double iHcal,double iEcal,int iEta,int iPhi) { 
  if(fNFrac == 1  && iCorr   < 0.1 ) return 0;
  if(fNFrac == 11 && iCorr   < 0.1 ) return 0; 
  fEcal = 0; if(iEcal > 0) fEcal = iEcal; 
  fHcal = 0; if(iHcal > 0) fHcal = iHcal; 
  if(iEcal < 0) fEcal = 0;
  if(abs(iEta) > fNEta/2-1 || iPhi < 1 || iPhi > fNPhi-1) return iCorr;
  double lFrac = fEcal/(fHcal+fEcal); 
  int    lIFrac=int(lFrac*10.);
  fIEta = iEta+fNEta/2;
  fIPhi = iPhi-1;
  fFrac = lFrac;
  if(lIFrac > fNFrac-1) lIFrac = 0; 
  fPtCorr = (fGraph[fIEta][fIPhi][lIFrac])->Eval(iCorr);
  fPtCorr = std::min(3.*fPtCorr,fPtCorr);
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
