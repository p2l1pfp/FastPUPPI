#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "FastPUPPI/NtupleProducer/interface/L1TPFParticle.h"
#include "FastPUPPI/NtupleProducer/interface/L1TPFUtils.h"


namespace l1tpf {
    class EcalProducerFromOfflineRechits : public edm::stream::EDProducer<> {
        public:
            explicit EcalProducerFromOfflineRechits(const edm::ParameterSet&) ;
            ~EcalProducerFromOfflineRechits() {}

        private:
            edm::EDGetTokenT<EcalRecHitCollection> src_;
            double etCut_, eCut_;
            edm::ESHandle<CaloGeometry> pG;

            virtual void produce(edm::Event&, const edm::EventSetup&) override;
            
            struct SimpleHit { 
                float et, eta, phi; 
                SimpleHit(float aet, float aeta, float aphi) : et(aet), eta(aeta), phi(aphi) {}
            };

    }; // class
} // namespace

l1tpf::EcalProducerFromOfflineRechits::EcalProducerFromOfflineRechits(const edm::ParameterSet & iConfig) :
    src_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("src"))),
    etCut_(iConfig.getParameter<double>("etMin")),
    eCut_(iConfig.getParameter<double>("eMin")) 
{
    produces<std::vector<l1tpf::Particle>>("crystals");
    produces<std::vector<l1tpf::Particle>>("towers");
}


void 
l1tpf::EcalProducerFromOfflineRechits::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) 
{
    std::unique_ptr<std::vector<l1tpf::Particle>> out_crystal(new std::vector<l1tpf::Particle>());
    std::unique_ptr<std::vector<l1tpf::Particle>> out_tower(new std::vector<l1tpf::Particle>());

    std::map<std::pair<int,int>,std::vector<SimpleHit>> towers;

    //Get Calo Geometry
    iSetup.get<CaloGeometryRecord>().get(pG);
    const CaloGeometry* caloGeom = pG.product();

    edm::Handle<EcalRecHitCollection> src;
    iEvent.getByToken(src_, src);
    for (const EcalRecHit & hit : *src) {
        if (hit.energy() <= eCut_) continue;
        if (hit.checkFlag(EcalRecHit::kOutOfTime) || hit.checkFlag(EcalRecHit::kL1SpikeFlag)) continue;
        const GlobalPoint & pos = caloGeom->getPosition(hit.detid());
        double et = pos.perp()/pos.mag() * hit.energy();
        if (et < etCut_) continue; 
        EBDetId id(hit.detid());
        out_crystal->emplace_back(et, pos.eta(), pos.phi(), 0., 0);
	//Using Local iEta,iPhi conventions since all the others are driving me nuts
	towers[make_pair(l1tpf::translateIEta(pos.eta()),l1tpf::translateIPhi(pos.phi(),pos.eta()))].emplace_back(et, pos.eta(), pos.phi());
	//std::cout << "check ieta " << id.tower_ieta() << " -- " << id.tower_iphi() << " -- " << l1tpf::translateIEta(pos.eta()) << " -- " << l1tpf::translateIPhi(pos.phi(),pos.eta()) << std::endl;
    }
    for (const auto & pair : towers) {
        double etsum = 0., etaetsum = 0., phietsum = 0.; 
        double reta = pair.second.front().eta, rphi = pair.second.front().phi;
        for (const SimpleHit & hit : pair.second) {
	  etsum += hit.et;
	  etaetsum += (hit.eta - reta) * hit.et;
	  phietsum += reco::deltaPhi(hit.phi, rphi) * hit.et;
        }
	if(etsum > 0)  {
	  etaetsum /= etsum;
	  phietsum /= etsum;
	}
	out_tower->emplace_back(etsum, etaetsum + reta, reco::deltaPhi(phietsum + rphi, 0.), 0., 0);
    }

  iEvent.put(std::move(out_crystal), "crystals");
  iEvent.put(std::move(out_tower), "towers");
}
using l1tpf::EcalProducerFromOfflineRechits;
DEFINE_FWK_MODULE(EcalProducerFromOfflineRechits);
