// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/L1TParticleFlow/interface/PFTrack.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TPFUtils.h"

#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"

#include <algorithm>

class L1DecTkTruthTableProducer : public edm::global::EDProducer<> {
public:
  explicit L1DecTkTruthTableProducer(const edm::ParameterSet&);
  ~L1DecTkTruthTableProducer();

private:
  virtual void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;
  //std::vector<l1t::PFTrack>

  std::string name_;
  edm::EDGetTokenT<std::vector<l1t::PFTrack>> decTks_;
  edm::EDGetTokenT<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>> ttTrackMCTruthToken_;
};

L1DecTkTruthTableProducer::L1DecTkTruthTableProducer(const edm::ParameterSet& iConfig)
    : name_(iConfig.getParameter<std::string>("name")),
      decTks_(consumes<std::vector<l1t::PFTrack>>(iConfig.getParameter<edm::InputTag>("src"))),
      ttTrackMCTruthToken_(consumes<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>>(
          iConfig.getParameter<edm::InputTag>("MCTruthTrackInputTag")))
{
  produces<nanoaod::FlatTable>();
}

L1DecTkTruthTableProducer::~L1DecTkTruthTableProducer() {}

// ------------ method called for each event  ------------
void L1DecTkTruthTableProducer::produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  edm::Handle<std::vector<l1t::PFTrack>> decTks;
  edm::Handle<TTTrackAssociationMap<Ref_Phase2TrackerDigi_>> MCTruthTTTrackHandle;
  iEvent.getByToken(decTks_, decTks);
  iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);

  // create the table
  unsigned int ncands = decTks->size();
  auto out = std::make_unique<nanoaod::FlatTable>(ncands, name_, false, true);

  std::vector<int> isGenuine(ncands), isLooselyGenuine(ncands), isUnknown(ncands), isCombinatoric(ncands), isReal(ncands);
  std::vector<unsigned int> tpmatch(ncands);

  for (unsigned int i = 0; i < ncands; ++i) {
    const auto tkPtr = edm::refToPtr((*decTks)[i].track());
    isGenuine[i] = MCTruthTTTrackHandle->isGenuine(tkPtr);
    isLooselyGenuine[i] = MCTruthTTTrackHandle->isLooselyGenuine(tkPtr);
    isUnknown[i] = MCTruthTTTrackHandle->isUnknown(tkPtr);
    isCombinatoric[i] = MCTruthTTTrackHandle->isCombinatoric(tkPtr);
    tpmatch[i] = (isCombinatoric[i]) + (isUnknown[i] << 1) + (isLooselyGenuine[i] << 2) + (isGenuine[i] << 3);

    edm::Ptr<TrackingParticle> my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(tkPtr);
    int tmp_isReal;
    if (my_tp.isNull()) {
      tmp_isReal = 0;
    } else {
      int tmp = (my_tp->eventId().event());
      if (tmp > 0) {
        tmp_isReal = 2;
      } else {
        tmp_isReal = 1;
      }
    }
    isReal[i] = tmp_isReal;
  }

  out->addColumn<int>("isGenuine", isGenuine, "");
  out->addColumn<int>("isLooselyGenuine", isLooselyGenuine, "");
  out->addColumn<int>("isUnknown", isUnknown, "");
  out->addColumn<int>("isCombinatoric", isCombinatoric, "");
  out->addColumn<int>("isReal", isReal, "");
  out->addColumn<int>("tpmatch", tpmatch, "");

  // save to the event branches
  iEvent.put(std::move(out));
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1DecTkTruthTableProducer);
