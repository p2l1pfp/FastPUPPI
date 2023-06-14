#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"
#include "DataFormats/L1Trigger/interface/VertexWord.h"

typedef SimpleFlatTableProducer<l1t::VertexWord> VertexWordFlatTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(VertexWordFlatTableProducer);