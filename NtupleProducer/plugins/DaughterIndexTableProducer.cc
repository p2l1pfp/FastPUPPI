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

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include <algorithm>

class DaughterIndexTableProducer : public edm::global::EDProducer<> {
public:
  explicit DaughterIndexTableProducer(const edm::ParameterSet&);
  ~DaughterIndexTableProducer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

  edm::EDGetTokenT<reco::CandidateView> srcJets_;
  edm::EDGetTokenT<reco::CandidateView> srcConstituents_;
  std::string tableName_, varName_, doc_;
};

DaughterIndexTableProducer::DaughterIndexTableProducer(const edm::ParameterSet& iConfig)
    : srcJets_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("jets"))),
      srcConstituents_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("constituents"))),
      tableName_(iConfig.getParameter<std::string>("tableName")),
      varName_(iConfig.getParameter<std::string>("varName")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  if (doc_.empty()) {
    doc_ = "Index of this element in " + iConfig.getParameter<edm::InputTag>("jets").encode() +
           " collection (-1 if not found)";
  }
  produces<nanoaod::FlatTable>();
}

DaughterIndexTableProducer::~DaughterIndexTableProducer() {}

// ------------ method called for each event  ------------
void DaughterIndexTableProducer::produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  edm::Handle<reco::CandidateView> srcJets;
  edm::Handle<reco::CandidateView> srcConstituents;
  iEvent.getByToken(srcJets_, srcJets);
  iEvent.getByToken(srcConstituents_, srcConstituents);
  auto ncands = srcConstituents->size();
  std::vector<int> indices(srcConstituents->size(), -1);
  for (size_t ijet = 0; ijet < srcJets->size(); ++ijet) {
    const auto& jet = (*srcJets)[ijet];
    for (size_t idau = 0, ndau = jet.numberOfSourceCandidatePtrs(); idau < ndau; ++idau) {
      const auto dau = jet.sourceCandidatePtr(idau);
      if (dau.isNull())
        throw cms::Exception("CorruptData") << "Null daughter found for jet index " << ijet << "\n";
      if ((dau.id() != srcConstituents.id() || dau.key() > ncands))
        throw cms::Exception("CorruptData") << "Daughter with invalid reference found for jet index " << ijet <<
            " id " << dau.id() << " (expected " << srcConstituents.id() << ")" <<
            " key " << dau.key() << "(expected <= " << ncands << ")\n";
      indices[dau.key()] = ijet;
    }
  }

  auto out = std::make_unique<nanoaod::FlatTable>(ncands, tableName_, /*singleton=*/false, /*extension=*/true);
  out->addColumn<int>(varName_, indices, doc_);

  iEvent.put(std::move(out));
}

void DaughterIndexTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription psetDesc;
  psetDesc.add<edm::InputTag>("jets");
  psetDesc.add<edm::InputTag>("constituents");
  psetDesc.add<std::string>("tableName");
  psetDesc.add<std::string>("varName", "jetIdx");
  psetDesc.add<std::string>("doc", "");
  descriptions.addWithDefaultLabel(psetDesc);
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DaughterIndexTableProducer);