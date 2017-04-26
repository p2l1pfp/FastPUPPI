#include "FastPUPPI/NtupleProducer/interface/L1TPFParticle.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "Rtypes.h" 
#include <vector>

namespace FastPUPPI_NtupleProducer {
    struct dictionary {
        edm::Wrapper<std::vector<l1tpf::Particle>> wvl1tpf;
    };
}
