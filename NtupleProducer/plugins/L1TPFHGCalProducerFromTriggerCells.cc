#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FastPUPPI/NtupleProducer/interface/L1TPFParticle.h"

namespace l1tpf {
    class HGCalProducerFromTriggerCells : public edm::stream::EDProducer<> {
        public:
            explicit HGCalProducerFromTriggerCells(const edm::ParameterSet&) ;
            ~HGCalProducerFromTriggerCells() {}

        private:
            edm::EDGetTokenT<l1t::HGCalTriggerCellBxCollection> src_;
            double etCut_;

            struct SimpleHit {
                float et, eta, phi;
                SimpleHit(float aet, float aeta, float aphi) : et(aet), eta(aeta), phi(aphi) {}
            };

            virtual void produce(edm::Event&, const edm::EventSetup&) override;
    }; // class
} // namespace

l1tpf::HGCalProducerFromTriggerCells::HGCalProducerFromTriggerCells(const edm::ParameterSet & iConfig) :
    src_(consumes<l1t::HGCalTriggerCellBxCollection>(iConfig.getParameter<edm::InputTag>("src"))),
    etCut_(iConfig.getParameter<double>("etMin"))
{
    produces<std::vector<l1tpf::Particle>>();
    produces<std::vector<l1tpf::Particle>>("towersEE");
    produces<std::vector<l1tpf::Particle>>("towersFHBH");
}


void 
l1tpf::HGCalProducerFromTriggerCells::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) 
{
    std::unique_ptr<std::vector<l1tpf::Particle>> out(new std::vector<l1tpf::Particle>());
    edm::Handle<l1t::HGCalTriggerCellBxCollection> multiclusters;
    iEvent.getByToken(src_, multiclusters);

    std::map<std::pair<int,int>,std::vector<SimpleHit>> towersEE, towersFHBH;

    for(auto it = multiclusters->begin(0), ed = multiclusters->end(0); it != ed; ++it) {
        if (it->pt() <= etCut_) continue;
        HGCalDetId id(it->detId());
        bool em = id.subdetId() == HGCEE;
        (em ? towersEE : towersFHBH)[make_pair(id.zside(), id.wafer())].emplace_back(it->pt(), it->eta(), it->phi());
        out->emplace_back(it->pt(), it->eta(), it->phi(), 0., em ? l1tpf::Particle::GAMMA : l1tpf::Particle::NH);
    }

    iEvent.put(std::move(out));

    for (int em = 0; em <= 1; ++em) {
        std::unique_ptr<std::vector<l1tpf::Particle>> out_tower(new std::vector<l1tpf::Particle>());
        for (const auto & pair : (em ? towersEE : towersFHBH)) {
            double etsum = 0., etaetsum = 0., phietsum = 0.;
            double reta = pair.second.front().eta, rphi = pair.second.front().phi;
            for (const SimpleHit & hit : pair.second) {
                etsum += hit.et;
                etaetsum += (hit.eta - reta) * hit.et;
                phietsum += reco::deltaPhi(hit.phi, rphi) * hit.et;
            }
            etaetsum /= etsum;
            phietsum /= etsum;
            float eta = etaetsum + reta, phi = reco::deltaPhi(phietsum + rphi, 0.);
            out_tower->emplace_back(etsum, eta, phi, 0.,  em ? l1tpf::Particle::GAMMA : l1tpf::Particle::NH);
        }
        iEvent.put(std::move(out_tower), em ? "towersEE" : "towersFHBH");
    }

}
using l1tpf::HGCalProducerFromTriggerCells;
DEFINE_FWK_MODULE(HGCalProducerFromTriggerCells);
