#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FastPUPPI/NtupleProducer/interface/L1TPFParticle.h"
#include "FastPUPPI/NtupleProducer/interface/L1TPFUtils.h"

namespace l1tpf {
    class TkProducerFromOfflineTracks : public edm::EDProducer {
        public:
            explicit TkProducerFromOfflineTracks(const edm::ParameterSet&) ;
            ~TkProducerFromOfflineTracks() {}

        private:
            edm::EDGetTokenT<reco::TrackCollection> TrackTag_;
            StringCutObjectSelector<reco::Track>    TrackSel_;
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

l1tpf::TkProducerFromOfflineTracks::TkProducerFromOfflineTracks(const edm::ParameterSet&iConfig) :
    TrackTag_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("TrackTag"))),
    TrackSel_(iConfig.getParameter<std::string>("TrackSel")),
    muonGunVeto_(iConfig.getParameter<bool>("MuonGunVeto"))
{
    if (muonGunVeto_) GenTagForVeto_ = consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("GenTagForVeto"));
    produces<std::vector<l1tpf::Particle>>();
}


void 
l1tpf::TkProducerFromOfflineTracks::produce(edm::Event & iEvent, const edm::EventSetup &) 
{
  std::unique_ptr<std::vector<l1tpf::Particle>> out(new std::vector<l1tpf::Particle>());
  std::vector<const reco::Candidate *> vetos;
  if (muonGunVeto_) {
      edm::Handle<reco::CandidateView> gens;
      iEvent.getByToken(GenTagForVeto_, gens);
      for (const reco::Candidate & c : *gens) vetos.push_back(&c);
  }
     
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(TrackTag_, tracks);
  for (const reco::Track &tk : *tracks) {
      if (!TrackSel_(tk)) continue;

      bool passveto = true;
      for (const reco::Candidate * veto : vetos) {
        if (::deltaR2(*veto, tk) < 0.01 && tk.pt() > 0.5 * std::min(veto->pt(),100.)) {
          passveto = false; break;
        }
      }
      if (!passveto) continue;

      const reco::Track::Vector & momentum = tk.momentum();
      const reco::Track::Point  & poca     = tk.referencePoint();  // point of closest approach
      bool iIsEle=false;
      double mass  = iIsEle ? 0.0005 : 0.137; 
      double energy = hypot(mass, tk.p());
      const math::XYZTLorentzVector tMom(momentum.X(),momentum.Y(),momentum.Z(),energy);
      const math::XYZTLorentzVector tVtx(poca.X()     ,poca.Y()     ,poca.Z(),0.);

      //Check below
      std::vector<double> lVars;
      l1tpf::propagate(1,lVars,tMom,tVtx,tk.charge(),fBz_);
      float trkEcalEta = lVars[4];
      float trkEcalPhi = lVars[5];
 
      int pQuality = 0; // note: these are made-up numbers, as I don't have yet a proper offline equivalent of stubs
      if (tk.hitPattern().stripLayersWithMeasurement() > 4) pQuality++;   
      if (tk.hitPattern().stripLayersWithMeasurement() > 6) pQuality += 2;
      if (tk.normalizedChi2() < 6.) pQuality += 4;  

      out->emplace_back(tk.pt(), tk.eta(), tk.phi(), 0.137, 
                        0, 0.0, poca.Z(), trkEcalEta, trkEcalPhi, tk.charge(), pQuality);
      out->back().setVertex(reco::Particle::Point(poca.x(),poca.y(),poca.z()));
  }
  iEvent.put(std::move(out));
}
using l1tpf::TkProducerFromOfflineTracks;
DEFINE_FWK_MODULE(TkProducerFromOfflineTracks);
