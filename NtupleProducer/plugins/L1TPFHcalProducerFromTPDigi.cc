#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"
#include "FastPUPPI/NtupleProducer/interface/L1TPFParticle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"

const double PI = 3.1415926535897;

namespace l1tpf {
    class HcalProducerFromTPDigi : public edm::EDProducer {
        public:
            explicit HcalProducerFromTPDigi(const edm::ParameterSet&) ;
            ~HcalProducerFromTPDigi() {}

        private:
            edm::EDGetTokenT<HcalTrigPrimDigiCollection> HcalTPTag_;
            std::unique_ptr<L1CaloHcalScale> hcalScale_;
            edm::ESHandle<CaloGeometry> calo;
            edm::ESHandle<HcalTrigTowerGeometry> theTrigTowerGeometry;


            virtual void produce(edm::Event&, const edm::EventSetup&) override;

    }; // class
} // namespace

l1tpf::HcalProducerFromTPDigi::HcalProducerFromTPDigi(const edm::ParameterSet & iConfig) :
    HcalTPTag_(consumes<HcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("HcalTPTag"))),
    hcalScale_(new L1CaloHcalScale(0.5)) // HCAL LSB = 0.5
{
    produces<std::vector<l1tpf::Particle>>();
}


void 
l1tpf::HcalProducerFromTPDigi::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) 
{
  std::unique_ptr<std::vector<l1tpf::Particle>> out(new std::vector<l1tpf::Particle>());

  iSetup.get<CaloGeometryRecord>().get(calo);
  const CaloGeometry* geo = calo.product(); 

  iSetup.get<CaloGeometryRecord>().get(theTrigTowerGeometry);
  const HcalTrigTowerGeometry* geoTrig = theTrigTowerGeometry.product();

  // / ----------------HCAL INFO-------------------
  // / Stealing some other code!
  // / https://github.com/cms-sw/cmssw/blob/0397259dd747cee94b68928f17976224c037057a/L1Trigger/L1TNtuples/src/L1AnalysisCaloTP.cc#L40
  // / 72 in phi and 56 in eta = 4032
  // / 144 TP is HF: 4 (in eta) * 18 (in phi) * 2 (sides)
  // / Hcal TPs
  edm::Handle< HcalTrigPrimDigiCollection > hcalTPs;
  iEvent.getByToken(HcalTPTag_, hcalTPs);
  for (HcalTrigPrimDigiCollection::const_iterator it = hcalTPs->begin(); it != hcalTPs->end(); ++it) {
      double et = hcalScale_->et( it->SOI_compressedEt(),it->id().ietaAbs(), it->id().zside() );

      // convert ieta-iphi to eta,phi 
      std::vector<HcalDetId> detids = geoTrig->detIds(it->id());
      unsigned int ndetsPerTower = detids.size();
      float towerEta = 0.;
      float towerPhi = 0.;
      float towerR   = 0.;
      float curphimax = -9999.;
      float curphimin = 9999.;
      for (unsigned int i = 0; i < ndetsPerTower; i++){
          const CaloCellGeometry *cell = geo->getGeometry( detids[i] );       
          // std::cout << detids[i].det() << "," << detids[i].subdetId() << "," << detids[i].iphi() << "," << detids[i].ieta() << "," << cell->etaPos() << "," << cell->phiPos() << std::endl;
          towerEta += cell->etaPos();
          towerR   += cell->getPosition().mag();
          if (curphimax < cell->phiPos()) curphimax = cell->phiPos();
          if (curphimin > cell->phiPos()) curphimin = cell->phiPos();
      }
      // std::cout << "curphimax = " << curphimax << ", curphimin = " << curphimin << std::endl;
      if (curphimax - curphimin < PI){ towerPhi = (curphimax + curphimin)/2.; }
      else{ towerPhi = -3.14159;} // special case
      // else{ towerPhi = (curphimax - curphimin + 2*PI)/2.; }
      // std::cout << "towerPhi = " << towerPhi << std::endl;
      // if (fabs(towerPhi) > PI && towerPhi > 0) towerPhi -= PI;
      // if (fabs(towerPhi) > PI && towerPhi < 0) towerPhi += PI;
      towerPhi *= 1.;
      towerEta /= float(ndetsPerTower);
      towerR /= float(ndetsPerTower);
      out->emplace_back( et, towerEta, towerPhi, 0,  0,0,0, towerEta, towerPhi, 0 );
      out->back().setIEtaIPhi( it->id().iphi(), it->id().ieta() );
  }

  iEvent.put(std::move(out));
}
using l1tpf::HcalProducerFromTPDigi;
DEFINE_FWK_MODULE(HcalProducerFromTPDigi);
