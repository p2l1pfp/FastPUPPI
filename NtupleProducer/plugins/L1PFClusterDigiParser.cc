// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/L1TParticleFlow/interface/PFCluster.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TPFUtils.h"


#include <algorithm>


// iso flag: two bits, least significant bit is the standalone WP (true or false), second bit is the looseTk WP (true or false)
// e.g. 0b01 : standalone iso flag passed, loose Tk iso flag did not pass
uint isoFlags(uint64_t digiWord) { return ((digiWord >> 38) & 0x3); }  // (two 1's) 0b11 = 0x3
bool passes_iso(uint64_t digiWord) { return (isoFlags(digiWord) & 0x1); }               // standalone iso WP
bool passes_looseTkiso(uint64_t digiWord) { return (isoFlags(digiWord) & 0x2); }        // loose Tk iso WP

// shower shape shape flag: two bits, least significant bit is the standalone WP, second bit is the looseTk WP
// e.g. 0b01 : standalone shower shape flag passed, loose Tk shower shape flag did not pass
uint shapeFlags(uint64_t digiWord) { return ((digiWord >> 51) & 0x3); }

bool passes_ss(uint64_t digiWord) { return (shapeFlags(digiWord) & 0x1); }         // standalone shower shape WP
bool passes_looseTkss(uint64_t digiWord) { return (shapeFlags(digiWord) & 0x2); }  // loose Tk shower shape WP

// brems: not saved in the current emulator
uint get_brems(uint64_t digiWord) { return ((digiWord >> 53) & 0x3); }

bool passes_standaloneWP(uint64_t digiWord)  { return (passes_iso(digiWord) && passes_ss(digiWord)); }
bool passes_looseL1TkMatchWP(uint64_t digiWord)  { return (passes_looseTkiso(digiWord) && passes_looseTkss(digiWord)); }

class L1PFClusterDigiParser : public edm::global::EDProducer<> {
public:
  explicit L1PFClusterDigiParser(const edm::ParameterSet&);
  ~L1PFClusterDigiParser();

private:
  virtual void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

  std::string name_;
  edm::EDGetTokenT<std::vector<l1t::PFCluster>> pfClusters_;

};

L1PFClusterDigiParser::L1PFClusterDigiParser(const edm::ParameterSet& iConfig)
    : name_(iConfig.getParameter<std::string>("name")),
      pfClusters_(consumes<std::vector<l1t::PFCluster>>(iConfig.getParameter<edm::InputTag>("src"))){
  produces<nanoaod::FlatTable>();
}

L1PFClusterDigiParser::~L1PFClusterDigiParser() {}

// ------------ method called for each event  ------------
void L1PFClusterDigiParser::produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  edm::Handle<std::vector<l1t::PFCluster>> pfClusters;
  iEvent.getByToken(pfClusters_, pfClusters);


  // create the table
  unsigned int ncands = pfClusters->size();
  auto out = std::make_unique<nanoaod::FlatTable>(ncands, name_, false, true);

  std::vector<uint> is_ss(ncands), is_iso(ncands), is_looseTkiso(ncands), is_looseTkss(ncands),brems(ncands), standaloneWP(ncands), looseL1TkMatchWP(ncands);


  for (unsigned int i = 0; i < ncands; ++i) {
    const auto pfCluster = (*pfClusters)[i];
    uint64_t digiWord = pfCluster.digiWord();
    is_iso[i] = passes_iso(digiWord);
    is_ss[i] = passes_ss(digiWord);
    is_looseTkiso[i] = passes_looseTkiso(digiWord);
    is_looseTkss[i] = passes_looseTkss(digiWord);
    brems[i] = get_brems(digiWord);
    standaloneWP[i] = passes_standaloneWP(digiWord);
    looseL1TkMatchWP[i] = passes_looseL1TkMatchWP(digiWord);



  }

  out->addColumn<uint>("isIso", is_iso, "");
  out->addColumn<uint>("isSS", is_ss, "");
  out->addColumn<uint>("isLooseTkIso", is_looseTkiso, "");
  out->addColumn<uint>("isLooseTkSS", is_looseTkss, "");
  out->addColumn<uint>("brems", brems, "");
  out->addColumn<uint>("standaloneWP", standaloneWP, "");
  out->addColumn<uint>("looseL1TkMatchWP", looseL1TkMatchWP, "");





  // save to the event branches
  iEvent.put(std::move(out));
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PFClusterDigiParser);
