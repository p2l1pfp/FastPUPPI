#ifndef FASTPUPPI_NTUPLERPRODUCER_ALTERNATIVEPF_H
#define FASTPUPPI_NTUPLERPRODUCER_ALTERNATIVEPF_H

// more PF algos
#include "FastPUPPI/NtupleProducer/interface/DiscretePF.h"

namespace l1tpf_int { 
class PFAlgo3 : public PFAlgo {
    public:
        PFAlgo3( const edm::ParameterSet& ) ;
        virtual void runPF(Region &r) const override;
    protected:
        float drMatchEm_, drMatchEmHad_;
        bool caloReLinkStep_; float caloReLinkDr_, caloReLinkThreshold_;
        bool sumTkCaloErr2_, ecalPriority_;
  };

} // end namespace

#endif
