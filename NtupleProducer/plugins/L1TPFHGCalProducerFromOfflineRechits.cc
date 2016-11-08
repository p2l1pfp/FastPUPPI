#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "FastPUPPI/NtupleProducer/interface/L1TPFParticle.h"


namespace l1tpf {
    class HGCalProducerFromOfflineRechits : public edm::stream::EDProducer<> {
        public:
            explicit HGCalProducerFromOfflineRechits(const edm::ParameterSet&) ;
            ~HGCalProducerFromOfflineRechits() {}

        private:
            edm::EDGetTokenT<HGCRecHitCollection> src_;
            double etCut_, eCut_;
            edm::ESHandle<CaloGeometry> pG;

            virtual void produce(edm::Event&, const edm::EventSetup&) override;
            
            struct SimpleHit { 
                float et, eta, phi; 
                SimpleHit(float aet, float aeta, float aphi) : et(aet), eta(aeta), phi(aphi) {}
            };

    }; // class
} // namespace

l1tpf::HGCalProducerFromOfflineRechits::HGCalProducerFromOfflineRechits(const edm::ParameterSet & iConfig) :
    src_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("src"))),
    etCut_(iConfig.getParameter<double>("etMin")),
    eCut_(iConfig.getParameter<double>("eMin")) 
{
    //produces<std::vector<l1tpf::Particle>>("crystals");
    produces<std::vector<l1tpf::Particle>>("towers");
}


void 
l1tpf::HGCalProducerFromOfflineRechits::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) 
{
    //std::unique_ptr<std::vector<l1tpf::Particle>> out_crystal(new std::vector<l1tpf::Particle>());
    std::unique_ptr<std::vector<l1tpf::Particle>> out_tower(new std::vector<l1tpf::Particle>());

    std::map<std::pair<int,int>,std::vector<SimpleHit>> towers;

    //Get Calo Geometry
    iSetup.get<CaloGeometryRecord>().get(pG);
    const CaloGeometry* caloGeom = pG.product();
    const HGCalGeometry* geomEE = (static_cast<const HGCalGeometry*>(caloGeom->getSubdetectorGeometry(HGCalDetId(HGCEE,1,1,0,0,0))));
    const HGCalGeometry* geomFH = (static_cast<const HGCalGeometry*>(caloGeom->getSubdetectorGeometry(HGCalDetId(HGCHEF,1,1,0,0,0))));

    edm::Handle<HGCRecHitCollection> src;
    iEvent.getByToken(src_, src);
    for (const HGCRecHit & hit : *src) {
        assert(hit.detid().det() == DetId::Forward && (3 <= hit.detid().subdetId() && hit.detid().subdetId() <= 4));
        if (hit.energy() < eCut_) continue;
        const GlobalPoint & pos = (hit.detid().subdetId() == 3 ? geomEE : geomFH)->getPosition(hit.detid());
        double et = pos.perp()/pos.mag() * hit.energy();
        if (et < etCut_) continue; 
        HGCalDetId detid(hit.detid());
        //out_crystal->emplace_back(et, pos.eta(), pos.phi(), 0, 0, 0, 0, pos.eta(), pos.phi());
        towers[make_pair(detid.zside(), detid.wafer())].emplace_back(et, pos.eta(), pos.phi());
    }
    for (const auto & pair : towers) {
        double etsum = 0., etaetsum = 0., phietsum = 0.; 
        double reta = pair.second.front().eta, rphi = pair.second.front().phi;
        for (const SimpleHit & hit : pair.second) {
            etsum += hit.et;
            etaetsum += (hit.eta - reta) * hit.et;
            phietsum += reco::deltaPhi(hit.phi, rphi) * hit.et;
        }
        etaetsum /= etsum;
        phietsum /= etsum;
        out_tower->emplace_back(etsum, etaetsum + reta, reco::deltaPhi(phietsum + rphi, 0.), 0, 0,0,0);
        //out_tower->back().setIEtaIPhi(pair.first.first, pair.first.second);
        //out_tower->back().setCaloEtaPhi( l1t::CaloTools::towerEta(pair.first.first), 
        //                                 l1t::CaloTools::towerPhi(pair.first.first, pair.first.second) );
    }

    //iEvent.put(std::move(out_crystal), "crystals");
    iEvent.put(std::move(out_tower), "towers");
}
using l1tpf::HGCalProducerFromOfflineRechits;
DEFINE_FWK_MODULE(HGCalProducerFromOfflineRechits);
