#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FastPUPPI/NtupleProducer/interface/L1TPFParticle.h"
#include "FastPUPPI/NtupleProducer/interface/L1TPFUtils.h"

#define HASL1TK
#ifdef HASL1TK
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"     
#endif


namespace l1tpf {
    class TkProducerFromL1Tracks : public edm::EDProducer {
        public:
            explicit TkProducerFromL1Tracks(const edm::ParameterSet&) ;
            ~TkProducerFromL1Tracks() {}

        private:
            edm::EDGetToken TrackTag_;
            bool muonGunVeto_;
            edm::EDGetTokenT<reco::CandidateView> GenTagForVeto_;
            float fBz_;

            virtual void produce(edm::Event&, const edm::EventSetup&) override;

            virtual void beginRun(edm::Run const&, edm::EventSetup const& iSetup) override {
                edm::ESHandle<MagneticField> magneticField;
                iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
                fBz_ = magneticField->inTesla(GlobalPoint(0,0,0)).z();
            }

    }; // class
} // namespace

l1tpf::TkProducerFromL1Tracks::TkProducerFromL1Tracks(const edm::ParameterSet&iConfig) :
    muonGunVeto_(iConfig.getParameter<bool>("MuonGunVeto"))
{
#ifdef HASL1TK
    TrackTag_ = consumes<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>>(iConfig.getParameter<edm::InputTag>("L1TrackTag"));
#endif
    if (muonGunVeto_) GenTagForVeto_ = consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("GenTagForVeto"));
    produces<std::vector<l1tpf::Particle>>();
}


void 
l1tpf::TkProducerFromL1Tracks::produce(edm::Event & iEvent, const edm::EventSetup &) 
{
  std::unique_ptr<std::vector<l1tpf::Particle>> out(new std::vector<l1tpf::Particle>());
  std::vector<const reco::Candidate *> vetos;
  if (muonGunVeto_) {
      edm::Handle<reco::CandidateView> gens;
      iEvent.getByToken(GenTagForVeto_, gens);
      for (const reco::Candidate & c : *gens) vetos.push_back(&c);
  }

#ifdef HASL1TK
 // https://github.com/skinnari/cmssw/blob/80c19f1b721325c3a02ee0482f72fb974a4c3bf7/L1Trigger/TrackFindingTracklet/test/L1TrackNtupleMaker.cc
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > tracks;
  iEvent.getByToken(TrackTag_, tracks);

  unsigned L1Tk_nPar = 4;
  for (const auto &tk : *tracks) {

      float pt   = tk.getMomentum(L1Tk_nPar).perp();
      float eta  = tk.getMomentum(L1Tk_nPar).eta();
      float phi  = tk.getMomentum(L1Tk_nPar).phi();
      float z0   = tk.getPOCA(L1Tk_nPar).z(); //cm
      int charge = tk.getRInv() > 0 ? +1 : -1;

      bool passveto = true;
      for (const reco::Candidate * veto : vetos) {
        if (::deltaR2(veto->eta(), veto->phi(), eta, phi) < 0.01 && pt > 0.5 * std::min(veto->pt(),100.)) {
          passveto = false; break;
        }
      }
      if (!passveto) continue;

      out->emplace_back(pt, eta, phi, 0.137, 0, 0.0, z0);
      l1tpf::Particle & me = out->back();
      // More info
      me.setQuality(tk.getStubRefs().size());
      me.setNormalizedChi2(tk.getChi2Red(L1Tk_nPar));
      me.setVertex(reco::Particle::Point(0,0,z0));
      me.setCharge(charge);

      // Calo propagation
      math::XYZTLorentzVector tMom(me.px(), me.py(), me.pz(), me.energy());
      math::XYZTLorentzVector tVtx(0,0,z0,0.);
      std::vector<double> lVars;
      l1tpf::propagate(1,lVars,tMom,tVtx,charge,fBz_);
      me.setCaloEtaPhi(lVars[4], lVars[5]); 
  }
#endif
  iEvent.put(std::move(out));
}
using l1tpf::TkProducerFromL1Tracks;
DEFINE_FWK_MODULE(TkProducerFromL1Tracks);
