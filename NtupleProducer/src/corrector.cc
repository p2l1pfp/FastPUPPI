#include "../interface/corrector.hh"
#include "../interface/L1TPFUtils.h"
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <TFile.h>
#include <TKey.h>
#include <TH1.h>
#include <TAxis.h>
#include "FWCore/Utilities/interface/CPUTimer.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <unordered_map>
#include <cassert>

corrector::corrector(const std::string iFile, int iNFrac, double iFracMax, int debug) :
    fNEta(l1tpf::towerNEta()),
    fNPhi(l1tpf::towerNPhi(1)),
    fNFrac(iNFrac),
    fFracMax(iFracMax)
{
  TFile *lFile = TFile::Open(iFile.c_str());
  assert(lFile);
  assert(fNFrac >= 1); // other cases are not supported yet
  fGraph = new TGraph***[fNEta];
  if (debug) lFile->Print();
  TH1 *index = (TH1*) lFile->Get("INDEX");
  TAxis *etaAxis = nullptr, *emfAxis = nullptr;
  if (index) {
    etaAxis = index->GetXaxis();
    if (iNFrac > 1) emfAxis = index->GetYaxis();
  }
  
  edm::CPUTimer timer;
  timer.start();
  char buff[1023];
  std::unordered_map<std::string,TGraph *> graphs;
  TKey *key;
  TIter nextkey(lFile->GetListOfKeys());
  while ((key = (TKey *) nextkey())) {
      if (strncmp(key->GetName(), "eta_", 4) == 0) {
          TGraph *gr = (TGraph*) key->ReadObj();
          if (!gr->TestBit(TGraph::kIsSortedX)) gr->Sort();
          graphs[key->GetName()] = gr;
      }
  }
  unsigned int ngraphs = 0;
  for(int i0 = 0; i0 < fNEta; i0++) { 
    fGraph[i0] = new TGraph**[fNPhi];
    for(int i1 = 0; i1 < fNPhi; i1++) { 
      fGraph[i0][i1] = new TGraph*[fNFrac];
      for(int i2 = 0; i2 < fNFrac; i2++) { 
        if (etaAxis) {
            float eta = std::min(std::abs(l1tpf::towerEta(i0-fNEta/2)), 4.999f);
            int   etaBin = std::min(std::max(etaAxis->FindBin(eta), 1), index->GetNbinsX());
            if (emfAxis) {
                float emf = i2/float(fNFrac-1);
                int   emfBin = std::min(std::max(emfAxis->FindBin(emf), 1), index->GetNbinsY());
                if (eta > 3.0) emfBin = 1; // no EMF bins in HF
                snprintf(buff, 1022, "eta_bin%d_emf_bin%d", etaBin, emfBin);
            } else {
                snprintf(buff, 1022, "eta_bin%d", etaBin);
            }
        } else if (iNFrac > 1) {
            snprintf(buff, 1022, "eta_%dphi_%d_frac_%d", i0-fNEta/2, i1+1, i2);
        } else {
            snprintf(buff, 1022, "eta_%d_frac_%d", i0, i2);
        }
        fGraph[i0][i1][i2] = graphs[buff]; 
        if (fGraph[i0][i1][i2] != nullptr) ngraphs++;
        //if (i1 == 0) printf("fGraph[%2d/%2d][%2d/%2d][%2d/%2d] is %s @ %p\n", i0, fNEta, i1, fNPhi, i2, fNFrac, buff, (const void *)(fGraph[i0][i1][i2]));
      }
    }
  }
  timer.stop();
  if (debug) std::cout << "Read " << ngraphs << " graphs from " << iFile << " in " << timer.realTime() << " s"  << std::endl; 
}
double corrector::correct(double iTotal,double iEcal,int iEta,int iPhi) { 
  if(fNFrac == 1  && iEcal       < 0.5 ) return 0;
  if(fNFrac >  1  && iTotal      < 1.0 ) return 0; 
  double fEcal  = 0; if(iEcal  > 0) fEcal  = iEcal; 
  double fTotal = 0; if(iTotal > 0) fTotal = iTotal; 
  if(iEcal < 0) fEcal = 0;
  if(abs(iEta) > fNEta/2-1 || iPhi < 1 || iPhi > fNPhi-1) return fTotal; //overflow
  if(fEcal > fTotal) fTotal = fEcal;
  double lFrac = fEcal/(fTotal); 
  int    lIFrac=int(floor(lFrac*(fNFrac-1)));
  int lIEta = iEta+fNEta/2;
  int lIPhi = iPhi-1;
  if(lIFrac > fNFrac-1) lIFrac = 0; 
  assert(lIEta >= 0 && lIEta < fNEta);
  assert(lIPhi >= 0 && lIPhi < fNPhi);
  assert(lIFrac >= 0 && lIFrac < fNFrac);
  if (!fGraph[lIEta][lIPhi][lIFrac]) {
    throw cms::Exception("RuntimeError") << "Error trying to read calibration [" << lIEta << "][" << lIPhi << "][" << lIFrac << "] for iTotal = " << iTotal << " iEcal = " << iEcal << " iEta " << iEta << " iPhi " << iPhi << std::endl;
  }
  double fPtCorr = (fGraph[lIEta][lIPhi][lIFrac])->Eval(fTotal);
  fPtCorr = std::min(4.*(fTotal),fPtCorr); // Just in case there is a bug don't go overboard
  if(fPtCorr < (iTotal > 0 ? 1. : 0.5)) fPtCorr = 0;
  ///======> Tuned parameters for optimal MET resolution
  if(iTotal > 0 && lFrac > fFracMax) return fTotal;   //Use just Ecal correction for Hi Ecal composition
  //if(fabs(l1tpf::towerEta(iEta)) > 3.0) fPtCorr *= 0.66;//
  //if(fabs(l1tpf::towerEta(iEta)) > 2.853 && fabs(l1tpf::towerEta(iEta)) < 3.1 && iHcal > 4.) fPtCorr = 0;
  return fPtCorr;
}

