#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "FastPUPPI/NtupleProducer/interface/L1TPFParticle.h"

namespace l1tpf {
    class HGCalProducerFrom3DTPs : public edm::stream::EDProducer<> {
        public:
            explicit HGCalProducerFrom3DTPs(const edm::ParameterSet&) ;
            ~HGCalProducerFrom3DTPs() {}

        private:
            edm::EDGetTokenT<l1t::HGCalMulticlusterBxCollection> src_;
            double etCut_;

            virtual void produce(edm::Event&, const edm::EventSetup&) override;
    }; // class
} // namespace

l1tpf::HGCalProducerFrom3DTPs::HGCalProducerFrom3DTPs(const edm::ParameterSet & iConfig) :
    src_(consumes<l1t::HGCalMulticlusterBxCollection>(iConfig.getParameter<edm::InputTag>("src"))),
    etCut_(iConfig.getParameter<double>("etMin"))
{
    produces<std::vector<l1tpf::Particle>>();
}


void 
l1tpf::HGCalProducerFrom3DTPs::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) 
{
    std::unique_ptr<std::vector<l1tpf::Particle>> out(new std::vector<l1tpf::Particle>());
    edm::Handle<l1t::HGCalMulticlusterBxCollection> multiclusters;
    iEvent.getByToken(src_, multiclusters);

    for(auto it = multiclusters->begin(0), ed = multiclusters->end(0); it != ed; ++it) {
        if (it->pt() <= etCut_) continue;
        out->emplace_back(it->pt(), it->eta(), it->phi(), 0, 0, 0, 0, it->eta(), it->phi());
        out->back().setHOverE(it->hOverE());
    }

    iEvent.put(std::move(out));
}
using l1tpf::HGCalProducerFrom3DTPs;
DEFINE_FWK_MODULE(HGCalProducerFrom3DTPs);
