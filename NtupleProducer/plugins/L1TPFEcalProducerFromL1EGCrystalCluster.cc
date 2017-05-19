#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "FastPUPPI/NtupleProducer/interface/L1TPFParticle.h"
#include "DataFormats/Phase2L1CaloTrig/interface/L1EGCrystalCluster.h"

namespace l1tpf {
    class EcalProducerFromL1EGCrystalCluster : public edm::stream::EDProducer<> {
        public:
            explicit EcalProducerFromL1EGCrystalCluster(const edm::ParameterSet&) ;
            ~EcalProducerFromL1EGCrystalCluster() {}

        private:
            edm::EDGetTokenT<l1slhc::L1EGCrystalClusterCollection> src_;
            double etCut_;

            virtual void produce(edm::Event&, const edm::EventSetup&) override;
    }; // class
} // namespace

l1tpf::EcalProducerFromL1EGCrystalCluster::EcalProducerFromL1EGCrystalCluster(const edm::ParameterSet & iConfig) :
    src_(consumes<l1slhc::L1EGCrystalClusterCollection>(iConfig.getParameter<edm::InputTag>("src"))),
    etCut_(iConfig.getParameter<double>("etMin"))
{
    produces<std::vector<l1tpf::Particle>>();
}


void 
l1tpf::EcalProducerFromL1EGCrystalCluster::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) 
{
    std::unique_ptr<std::vector<l1tpf::Particle>> out(new std::vector<l1tpf::Particle>());
    edm::Handle<l1slhc::L1EGCrystalClusterCollection> clusters;
    iEvent.getByToken(src_, clusters);

    for(auto it = clusters->begin(), ed = clusters->end(); it != ed; ++it) {
        if (it->pt() <= etCut_) continue;
        out->emplace_back(it->pt(), it->eta(), it->phi(), 0., 0);
        out->back().setHOverE(it->hovere());
    }

    iEvent.put(std::move(out));
}
using l1tpf::EcalProducerFromL1EGCrystalCluster;
DEFINE_FWK_MODULE(EcalProducerFromL1EGCrystalCluster);
