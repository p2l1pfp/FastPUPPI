#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"
#include "CondFormats/L1TObjects/interface/L1CaloEcalScale.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "FastPUPPI/NtupleProducer/interface/L1TPFParticle.h"
#include "FastPUPPI/NtupleProducer/interface/L1TPFUtils.h"

namespace l1tpf {
    class EcalProducerFromTPDigi : public edm::EDProducer {
        public:
            explicit EcalProducerFromTPDigi(const edm::ParameterSet&) ;
            ~EcalProducerFromTPDigi() {}

        private:
            edm::EDGetTokenT<EcalTrigPrimDigiCollection> EcalTPTag_;
            std::unique_ptr<L1CaloEcalScale> ecalScale_;

            virtual void produce(edm::Event&, const edm::EventSetup&) override;

    }; // class
} // namespace

l1tpf::EcalProducerFromTPDigi::EcalProducerFromTPDigi(const edm::ParameterSet & iConfig) :
    EcalTPTag_(consumes<EcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("EcalTPTag"))),
    ecalScale_(new L1CaloEcalScale(0.5)) // ECAL LSB = 0.5
{
    produces<std::vector<l1tpf::Particle>>();
}


void 
l1tpf::EcalProducerFromTPDigi::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) 
{
  std::unique_ptr<std::vector<l1tpf::Particle>> out(new std::vector<l1tpf::Particle>());

  /// ----------------ECAL INFO-------------------
  /// 72 in phi and 56 in eta
  /// Ecal TPs
  edm::Handle< EcalTrigPrimDigiCollection > ecalTPs;
  iEvent.getByToken(EcalTPTag_, ecalTPs);
  for (EcalTrigPrimDigiCollection::const_iterator it = ecalTPs->begin(); it != ecalTPs->end(); ++it) {
      short sign = 1; if(it->id().ieta() < 0) sign=-1;
      double et = ecalScale_->et( it->compressedEt(),(unsigned short)abs(it->id().ieta()), sign);
      if(et < 0.1) continue;
      //!!!!Check ME
      float curTowerEta = l1tpf::towerEta(it->id().ieta());
      float curTowerPhi = l1tpf::towerPhi(it->id().ieta(),it->id().iphi());
      out->emplace_back( et, curTowerEta, curTowerPhi, 0,  0,0,0, curTowerEta, curTowerPhi, 0 );
  }

  iEvent.put(std::move(out));
}
using l1tpf::EcalProducerFromTPDigi;
DEFINE_FWK_MODULE(EcalProducerFromTPDigi);
