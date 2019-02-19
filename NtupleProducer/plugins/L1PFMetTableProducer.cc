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
          bool central;
          
          MetBranch(const std::string &n, const edm::EDGetTokenT<reco::CandidateView> &t) : 
              name(n), token(t), central(name.find("Central") != std::string::npos) {}

      };
      std::vector<MetBranch> mets_;

      edm::EDGetTokenT<reco::CandidateView> genMetToken_, genMetCentralToken_;
};

L1PFMetTableProducer::L1PFMetTableProducer(const edm::ParameterSet& iConfig) :
    genMetToken_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("genMet"))),
    genMetCentralToken_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("genMetCentral")))
{
    produces<nanoaod::FlatTable>("genMet");
    produces<nanoaod::FlatTable>("genCentralMet");
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
    edm::Handle<reco::CandidateView> hgenMet, hgenMetCentral, hmet;
    
    iEvent.getByToken(genMetToken_, hgenMet);
    iEvent.getByToken(genMetCentralToken_, hgenMetCentral);

    const reco::Candidate & genMet = (*hgenMet)[0];
    const reco::Candidate & genMetCentral = (*hgenMetCentral)[0];

    auto outg = std::make_unique<nanoaod::FlatTable>(1, "genMet", true); 
    auto outc = std::make_unique<nanoaod::FlatTable>(1, "genCentralMet", true);
    outg->addColumnValue<float>("pt",  genMet.pt(), "genMet pt", nanoaod::FlatTable::FloatColumn);
    outg->addColumnValue<float>("phi", genMet.phi(), "genMet phi", nanoaod::FlatTable::FloatColumn);
    outc->addColumnValue<float>("pt",  genMetCentral.pt(), "genCentralMet pt", nanoaod::FlatTable::FloatColumn);
    outc->addColumnValue<float>("phi", genMetCentral.phi(), "genCentralMet phi", nanoaod::FlatTable::FloatColumn);
    iEvent.put(std::move(outg), "genMet");
    iEvent.put(std::move(outc), "genCentralMet");

    for (auto & m : mets_) {
        auto out = std::make_unique<nanoaod::FlatTable>(1, m.name+"Met", true); 
        iEvent.getByToken(m.token, hmet);
        const reco::Candidate & rec = (*hmet)[0];
        const reco::Candidate & gen = m.central ? genMetCentral : genMet;
        out->addColumnValue<float>("pt", rec.pt(), m.name + "Met pt", nanoaod::FlatTable::FloatColumn);
        out->addColumnValue<float>("phi", rec.phi(), m.name + "Met phi", nanoaod::FlatTable::FloatColumn);
        out->addColumnValue<float>("para", rec.pt() * std::cos(rec.phi() - gen.phi()), m.name + "Met component parallel to gen", nanoaod::FlatTable::FloatColumn);
        out->addColumnValue<float>("perp", rec.pt() * std::sin(rec.phi() - gen.phi()), m.name + "Met component perpendicular to gen", nanoaod::FlatTable::FloatColumn);
        iEvent.put(std::move(out), m.name+"Met");
    }

}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PFMetTableProducer);
