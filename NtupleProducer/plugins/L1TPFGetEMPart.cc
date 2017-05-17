#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FastPUPPI/NtupleProducer/interface/L1TPFParticle.h"

namespace l1tpf {
    class GetEMPart : public edm::stream::EDProducer<> {
        public:
            explicit GetEMPart(const edm::ParameterSet&) ;
            ~GetEMPart() {}

        private:
            edm::EDGetTokenT<std::vector<l1tpf::Particle>> src_;

            virtual void produce(edm::Event&, const edm::EventSetup&) override;
    }; // class
} // namespace

l1tpf::GetEMPart::GetEMPart(const edm::ParameterSet & iConfig) :
    src_(consumes<std::vector<l1tpf::Particle>>(iConfig.getParameter<edm::InputTag>("src")))
{
    produces<std::vector<l1tpf::Particle>>();
}


void 
l1tpf::GetEMPart::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) 
{
    std::unique_ptr<std::vector<l1tpf::Particle>> out(new std::vector<l1tpf::Particle>());
    edm::Handle<std::vector<l1tpf::Particle>> src;
    iEvent.getByToken(src_, src);
    for (const l1tpf::Particle & p : *src) {
        if (p.emEt() > 0) {
            out->push_back(p);
            out->back().setPt(p.emEt());
        }
    }
    iEvent.put(std::move(out));
}
using l1tpf::GetEMPart;
DEFINE_FWK_MODULE(GetEMPart);
