#ifndef FastPUPPI_NtupleProducer_L1TPFUtils_h
#define FastPUPPI_NtupleProducer_L1TPFUtils_h
#include <vector>
#include <DataFormats/Math/interface/LorentzVector.h>

namespace l1tpf {
   void  propagate(int iOption,std::vector<double> &iVars,const math::XYZTLorentzVector& iMom,const math::XYZTLorentzVector& iVtx,double iCharge,double iBField) ;
   std::pair<float,float> towerEtaBounds(int ieta);
   float towerEtaSize(int ieta);
   float towerPhiSize(int ieta);
   int   towerNPhi(int ieta);
   int   translateIEta(float eta);
   int   translateIPhi(float phi,float eta);
   float towerEta(int ieta);
   float towerPhi(int ieta, int iphi);
   //Map to array
   int   translateAEta(int ieta,bool iInvert=false);
   int   translateAPhi(int iphi,bool iInvert=false);
}

#endif

