#ifndef FASTPUPPI_NTUPLERPRODUCER_BITWISEPF_H
#define FASTPUPPI_NTUPLERPRODUCER_BITWISEPF_H

#include "FastPUPPI/NtupleProducer/interface/DiscretePF.h"

namespace l1tpf_int { 
class BitwisePF : public PFAlgo {
    public:
        BitwisePF( const edm::ParameterSet& ) ;
        virtual void runPF(Region &r) const override;

  };

} // end namespace

#endif
