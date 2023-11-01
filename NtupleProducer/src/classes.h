#include "DataFormats/Common/interface/Wrapper.h"
#include "Rtypes.h" 
#include <vector>
#include "TLorentzVector.h"

namespace FastPUPPI_NtupleProducer {
    struct dictionary {
        std::vector<TLorentzVector> vlv;
        std::vector<std::pair<TLorentzVector,int>> vplvi;
        std::vector<std::pair<float,int>> vpfi;
        std::vector<std::tuple<float,float,int>> vtffi;
        std::vector<std::tuple<float,float,float,int>> vtfffi;
        std::vector<std::map<std::pair<int,int>,float>> vmpiif;
    };
}
