#ifndef FastPUPPI_NtupleProducer_L1TPFUtils_h
#define FastPUPPI_NtupleProducer_L1TPFUtils_h
#include <vector>
#include <DataFormats/Math/interface/LorentzVector.h>

namespace l1tpf {
   void propagate(int iOption,std::vector<double> &iVars,const math::XYZTLorentzVector& iMom,const math::XYZTLorentzVector& iVtx,double iCharge,double iBField) ;
}

#endif

