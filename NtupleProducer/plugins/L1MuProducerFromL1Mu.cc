#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FastPUPPI/NtupleProducer/interface/L1TPFParticle.h"
#include "FastPUPPI/NtupleProducer/interface/L1TPFUtils.h"

namespace l1tpf {
    class MuProducerFromL1Mu : public edm::EDProducer {
        public:
            explicit MuProducerFromL1Mu(const edm::ParameterSet&) ;
            ~MuProducerFromL1Mu() {}

        private:
            edm::EDGetTokenT<l1t::MuonBxCollection> MuonTag_;
            bool muonGunVeto_;
            edm::EDGetTokenT<reco::CandidateView> GenTagForVeto_;

            virtual void produce(edm::Event&, const edm::EventSetup&) override;

    }; // class
} // namespace

l1tpf::MuProducerFromL1Mu::MuProducerFromL1Mu(const edm::ParameterSet&iConfig) :
    MuonTag_(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("MuonTag"))),
    muonGunVeto_(iConfig.getParameter<bool>("MuonGunVeto"))
{
    if (muonGunVeto_) GenTagForVeto_ = consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("GenTagForVeto"));
    produces<std::vector<l1tpf::Particle>>();
}


void 
l1tpf::MuProducerFromL1Mu::produce(edm::Event & iEvent, const edm::EventSetup &) 
{
  std::unique_ptr<std::vector<l1tpf::Particle>> out(new std::vector<l1tpf::Particle>());
  std::vector<const reco::Candidate *> vetos;
  if (muonGunVeto_) {
      edm::Handle<reco::CandidateView> gens;
      iEvent.getByToken(GenTagForVeto_, gens);
      for (const reco::Candidate & c : *gens) vetos.push_back(&c);
  }
     
  edm::Handle<l1t::MuonBxCollection> muon;    
  iEvent.getByToken(MuonTag_, muon);
  if (muon.isValid()){ 
      for (int ibx = muon->getFirstBX(); ibx <= muon->getLastBX(); ++ibx) {
          if (ibx != 0) continue; // only the first bunch crossing

          for (auto it=muon->begin(ibx); it!=muon->end(ibx); it++){      
              if (it->et() == 0) continue; // if you don't care about L1T candidates with zero ET.

              l1tpf::Particle mu( it->pt(), it->eta(), ::deltaPhi(it->phi(), 0.), 0.105, 4 ); // the deltaPhi is to get the wrapping correct, 4 is the muon ID code
              mu.setCharge(it->charge());
              mu.setQuality(it->hwQual());

              bool passveto = true;
              for (const reco::Candidate * veto : vetos) {
                  if (::deltaR2(*veto, mu) < 0.5) {
                      passveto = false; break;
                  }
              }
              if (!passveto) continue;

              out->push_back(mu);
          }
      }
  } else {
      edm::LogWarning("MissingProduct") << "L1Upgrade muon bx collection not found." << std::endl;
  }
  iEvent.put(std::move(out));
}
using l1tpf::MuProducerFromL1Mu;
DEFINE_FWK_MODULE(MuProducerFromL1Mu);
