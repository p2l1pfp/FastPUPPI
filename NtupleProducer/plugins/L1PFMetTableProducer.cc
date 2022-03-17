// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

class L1PFMetTableProducer : public edm::global::EDProducer<>  {
   public:
      explicit L1PFMetTableProducer(const edm::ParameterSet&);
      ~L1PFMetTableProducer();

   private:
       virtual void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;


      struct MetBranch {
          std::string name;
          edm::EDGetTokenT<reco::CandidateView> token;
          
          MetBranch(const std::string &n, const edm::EDGetTokenT<reco::CandidateView> &t) : 
              name(n), token(t) {}

      };
      std::vector<MetBranch> mets_;

      edm::EDGetTokenT<reco::CandidateView> genMetToken_;
      std::string flavour_;
};

L1PFMetTableProducer::L1PFMetTableProducer(const edm::ParameterSet& iConfig) :
    genMetToken_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("genMet"))),
    flavour_(iConfig.getParameter<std::string>("flavour"))
{
    produces<nanoaod::FlatTable>("genMet");
    edm::ParameterSet mets = iConfig.getParameter<edm::ParameterSet>("mets");
    for (const std::string & name : mets.getParameterNamesForType<edm::InputTag>()) {
        mets_.emplace_back(name, consumes<reco::CandidateView>(mets.getParameter<edm::InputTag>(name)));
        produces<nanoaod::FlatTable>(name+"Met");
    }
}

L1PFMetTableProducer::~L1PFMetTableProducer() { }


// ------------ method called for each event  ------------
void
L1PFMetTableProducer::produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const 
{
    edm::Handle<reco::CandidateView> hgenMet, hmet;
    
    iEvent.getByToken(genMetToken_, hgenMet);

    const reco::Candidate & gen = (*hgenMet)[0];

    auto outg = std::make_unique<nanoaod::FlatTable>(1, "genMet"+flavour_, true); 
    outg->addColumnValue<float>("pt",  gen.pt(), "genMet pt");
    outg->addColumnValue<float>("phi", gen.phi(), "genMet phi");
    iEvent.put(std::move(outg), "genMet");

    for (auto & m : mets_) {
        auto out = std::make_unique<nanoaod::FlatTable>(1, m.name+"Met"+flavour_, true); 
        iEvent.getByToken(m.token, hmet);
        const reco::Candidate & rec = (*hmet)[0];
        out->addColumnValue<float>("pt", rec.pt(), m.name + "Met pt");
        out->addColumnValue<float>("phi", rec.phi(), m.name + "Met phi");
        out->addColumnValue<float>("para", rec.pt() * std::cos(rec.phi() - gen.phi()), m.name + "Met component parallel to gen");
        out->addColumnValue<float>("perp", rec.pt() * std::sin(rec.phi() - gen.phi()), m.name + "Met component perpendicular to gen");
        iEvent.put(std::move(out), m.name+"Met");
    }

}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PFMetTableProducer);
