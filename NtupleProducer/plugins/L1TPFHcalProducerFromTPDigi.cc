#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"
#include "FastPUPPI/NtupleProducer/interface/L1TPFParticle.h"
//#include "FastPUPPI/NtupleProducer/interface/L1TPFUtils.h"
#include "FWCore/Framework/interface/ESHandle.h"
//#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
//#include "Geometry/Records/interface/CaloGeometryRecord.h"
//include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
//#include "DataFormats/DetId/interface/DetId.h"
//#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
//#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"


namespace l1tpf {
    class HcalProducerFromTPDigi : public edm::EDProducer {
        public:
            explicit HcalProducerFromTPDigi(const edm::ParameterSet&) ;
            ~HcalProducerFromTPDigi() {}

        private:
            edm::EDGetTokenT<HcalTrigPrimDigiCollection> HcalTPTag_;
            edm::ESHandle<CaloTPGTranscoder> decoder_;
            //edm::ESHandle<CaloGeometry> calo;
            //edm::ESHandle<HcalTrigTowerGeometry> theTrigTowerGeometry;


            virtual void produce(edm::Event&, const edm::EventSetup&) override;

    }; // class
} // namespace

l1tpf::HcalProducerFromTPDigi::HcalProducerFromTPDigi(const edm::ParameterSet & iConfig) :
    HcalTPTag_(consumes<HcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("HcalTPTag")))
{
    produces<std::vector<l1tpf::Particle>>();
}


// https://github.com/cms-l1t-offline/cmssw/blob/phase2-l1t-integration-CMSSW_9_1_0_pre2/L1Trigger/L1TNtuples/plugins/L1CaloTowerTreeProducer.cc
void 
l1tpf::HcalProducerFromTPDigi::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) 
{
  std::unique_ptr<std::vector<l1tpf::Particle>> out(new std::vector<l1tpf::Particle>());

  iSetup.get<CaloTPGRecord>().get(decoder_);

  //iSetup.get<CaloGeometryRecord>().get(calo);
  //const CaloGeometry* geo = calo.product(); 
  //iSetup.get<CaloGeometryRecord>().get(theTrigTowerGeometry);
  //const HcalTrigTowerGeometry* geoTrig = theTrigTowerGeometry.product();

  edm::Handle< HcalTrigPrimDigiCollection > hcalTPs;
  iEvent.getByToken(HcalTPTag_, hcalTPs);
  for (const auto & itr : *hcalTPs) {
      HcalTrigTowerDetId id = itr.id();
      double et = decoder_->hcaletValue(itr.id(), itr.t0());
      if (et <= 0) continue;
      float towerEta = l1t::CaloTools::towerEta(id.ieta());
      float towerPhi = l1t::CaloTools::towerPhi(id.ieta(), id.iphi());
      out->emplace_back( et, towerEta, towerPhi, 0,  0,0,0, towerEta, towerPhi, 0 );
  }

  iEvent.put(std::move(out));
}
using l1tpf::HcalProducerFromTPDigi;
DEFINE_FWK_MODULE(HcalProducerFromTPDigi);
