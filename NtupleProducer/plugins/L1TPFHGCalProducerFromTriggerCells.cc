#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "FastPUPPI/NtupleProducer/interface/L1TPFParticle.h"

namespace l1tpf {
    class HGCalProducerFromTriggerCells : public edm::stream::EDProducer<> {
        public:
            explicit HGCalProducerFromTriggerCells(const edm::ParameterSet&) ;
            ~HGCalProducerFromTriggerCells() {}

        private:
            edm::EDGetTokenT<l1t::HGCalTriggerCellBxCollection> src_;
            double etCut_;

            virtual void produce(edm::Event&, const edm::EventSetup&) override;
    }; // class
} // namespace

l1tpf::HGCalProducerFromTriggerCells::HGCalProducerFromTriggerCells(const edm::ParameterSet & iConfig) :
    src_(consumes<l1t::HGCalTriggerCellBxCollection>(iConfig.getParameter<edm::InputTag>("src"))),
    etCut_(iConfig.getParameter<double>("etMin"))
{
    produces<std::vector<l1tpf::Particle>>();
}


void 
l1tpf::HGCalProducerFromTriggerCells::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) 
{
    std::unique_ptr<std::vector<l1tpf::Particle>> out(new std::vector<l1tpf::Particle>());
    edm::Handle<l1t::HGCalTriggerCellBxCollection> multiclusters;
    iEvent.getByToken(src_, multiclusters);

    for(auto it = multiclusters->begin(0), ed = multiclusters->end(0); it != ed; ++it) {
        if (it->pt() <= etCut_) continue;
        out->emplace_back(it->pt(), it->eta(), it->phi(), 0., 0., it->layer(), 0, it->eta(), it->phi());
    }

    iEvent.put(std::move(out));
}
using l1tpf::HGCalProducerFromTriggerCells;
DEFINE_FWK_MODULE(HGCalProducerFromTriggerCells);
