// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFJet.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "DataFormats/L1TCorrelator/interface/TkElectron.h"
#include "DataFormats/L1TCorrelator/interface/TkEm.h"

#include <cstdio>
#include <cstdint>
#include <fstream>

class PuppiDumperHelper {
public:
  PuppiDumperHelper(const edm::ParameterSet &cfg, edm::ConsumesCollector cc)
      : src_(cc.consumes<edm::View<l1t::PFCandidate>>(cfg.getParameter<edm::InputTag>("src"))) {}
  void dump(const edm::Event &iEvent, std::fstream &out) {
    edm::Handle<edm::View<l1t::PFCandidate>> src;
    iEvent.getByToken(src_, src);
    uint64_t nobj = src->size();
    nobj |= (0b10llu << 62);  // event header
    out.write(reinterpret_cast<const char *>(&nobj), sizeof(uint64_t));
    for (auto &c : *src) {
      uint64_t packed = c.encodedPuppi64();
      out.write(reinterpret_cast<const char *>(&packed), sizeof(uint64_t));
    }
  }

private:
  edm::EDGetTokenT<edm::View<l1t::PFCandidate>> src_;
};

class JetDumperHelper {
public:
  JetDumperHelper(const edm::ParameterSet &cfg, edm::ConsumesCollector cc)
      : src_(cc.consumes<edm::View<l1t::PFJet>>(cfg.getParameter<edm::InputTag>("src"))),
        ptMin_(cfg.getParameter<double>("ptMin")) {}
  void dump(const edm::Event &iEvent, std::fstream &out) {
    edm::Handle<edm::View<l1t::PFJet>> src;
    iEvent.getByToken(src_, src);
    std::vector<uint64_t> data(1, 0u);  // leave one empty word at the beginning
    for (l1t::PFJet c : *src) {
      if (c.pt() > ptMin_) {
        const std::array<uint64_t, 2> &enc = c.encodedJet();
        data.push_back(enc[0]);
        data.push_back(enc[1]);
      }
    }
    // now fill the size
    data[0] = (data.size() - 1);
    data[0] |= (0b10llu << 62);  // event header
    out.write(reinterpret_cast<const char *>(&data[0]), data.size() * sizeof(uint64_t));
  }

private:
  edm::EDGetTokenT<edm::View<l1t::PFJet>> src_;
  float ptMin_;
};

class TrackerMuonDumperHelper {
public:
  TrackerMuonDumperHelper(const edm::ParameterSet &cfg, edm::ConsumesCollector cc)
      : src_(cc.consumes<edm::View<l1t::TrackerMuon>>(cfg.getParameter<edm::InputTag>("src"))) {}
  void dump(const edm::Event &iEvent, std::fstream &out) {
    edm::Handle<edm::View<l1t::TrackerMuon>> src;
    iEvent.getByToken(src_, src);
    std::vector<uint64_t> data(1, 0u);  // leave one empty word at the beginning
    bool even = true;
    for (auto &c : *src) {
      const std::array<uint64_t, 2> &enc = c.word();
      // unfortunately currently in CMSSW the first word is missing the valid bit and is shifted
      // we fix it on the fly
      uint64_t enc0fixed = (enc[0] << 1) | 1;
      if (even) {
        data.push_back(enc0fixed);
        data.push_back(enc[1]);
      } else {
        data.back() |= (enc[1] << 32);
        data.push_back(enc0fixed);
      }
      even = !even;
    }
    data[0] = (data.size() - 1);
    data[0] |= (0b10llu << 62);  // event header
    out.write(reinterpret_cast<const char *>(&data[0]), data.size() * sizeof(uint64_t));
  }

private:
  edm::EDGetTokenT<edm::View<l1t::TrackerMuon>> src_;
};

class CTL2EgammaDumperHelper {
public:
  CTL2EgammaDumperHelper(const edm::ParameterSet &cfg, edm::ConsumesCollector cc)
      : srcEle_(cc.consumes<edm::View<l1t::TkElectron>>(cfg.getParameter<edm::InputTag>("srcEle"))),
        srcEm_(cc.consumes<edm::View<l1t::TkEm>>(cfg.getParameter<edm::InputTag>("srcEm"))),
        interleave_(cfg.getParameter<bool>("interleaveOutputs")) {}
  void dump(const edm::Event &iEvent, std::fstream &out) {
    edm::Handle<edm::View<l1t::TkElectron>> srcEle;
    iEvent.getByToken(srcEle_, srcEle);
    edm::Handle<edm::View<l1t::TkEm>> srcEm;
    iEvent.getByToken(srcEm_, srcEm);
    std::vector<uint64_t> data(1, 0u);  // leave one empty word at the beginning
    if (interleave_) {
      for (unsigned int i = 0, nEle = srcEle->size(), nEm = srcEm->size(), nMax = std::max(nEle, nEm); i < nMax; ++i) {
        // photons first
        ap_uint<96> wEle = (i < nEle ? (*srcEle)[i].egBinaryWord<96>() : ap_uint<96>(0));
        ap_uint<96> wEm = (i < nEm ? (*srcEm)[i].egBinaryWord<96>() : ap_uint<96>(0));
        ap_uint<64> word;
        word = wEm(63, 0);
        data.push_back(word.to_uint64());
        word(31, 0) = wEm(95, 64);
        word(63, 32) = wEle(95, 64);
        data.push_back(word.to_uint64());
        word = wEle(63, 0);
        data.push_back(word.to_uint64());
      }
    } else {
      // photons first, zero-padded to 12
      ap_uint<64> word = 0;
      for (unsigned int i = 0, nEm = srcEm->size(); i < 12; ++i) {
        ap_uint<96> wEm = (i < nEm ? (*srcEm)[i].egBinaryWord<96>() : ap_uint<96>(0));
        if (i % 2 == 0) {
          word = wEm(63, 0);
          data.push_back(word.to_uint64());
          word(31, 0) = wEm(95, 64);
        } else {
          word(63, 32) = wEm(95, 64);
          data.push_back(word.to_uint64());
          word = wEm(63, 0);
          data.push_back(word.to_uint64());
        }
      }
      // then we do electrons
      word = 0;
      for (unsigned int i = 0, nEle = srcEle->size(); i < nEle; ++i) {
        ap_uint<96> wEle = (*srcEle)[i].egBinaryWord<96>();
        if (i % 2 == 0) {
          word = wEle(63, 0);
          data.push_back(word.to_uint64());
          word(31, 0) = wEle(95, 64);
        } else {
          word(63, 32) = wEle(95, 64);
          data.push_back(word.to_uint64());
          word = wEle(63, 0);
          data.push_back(word.to_uint64());
          word = 0;
        }
      }
      if (word != 0)
        data.push_back(word.to_uint64());
    }
    data[0] = (data.size() - 1);
    data[0] |= (0b10llu << 62);  // event header
    out.write(reinterpret_cast<const char *>(&data[0]), data.size() * sizeof(uint64_t));
  }

private:
  edm::EDGetTokenT<edm::View<l1t::TkElectron>> srcEle_;
  edm::EDGetTokenT<edm::View<l1t::TkEm>> srcEm_;
  bool interleave_;
};

template <typename Helper>
class BinaryDumper : public edm::one::EDAnalyzer<> {
public:
  explicit BinaryDumper(const edm::ParameterSet &iConfig)
      : helper_(iConfig, consumesCollector()), outName_(iConfig.getParameter<std::string>("outName")), out_() {}

private:
  Helper helper_;
  std::string outName_;
  std::fstream out_;

  void beginJob() override { out_.open(outName_, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc); }
  void analyze(const edm::Event &iEvent, const edm::EventSetup &) override { helper_.dump(iEvent, out_); }
  void endJob() override { out_.close(); }
};

typedef BinaryDumper<PuppiDumperHelper> L1PuppiBinaryDumper;
typedef BinaryDumper<JetDumperHelper> L1JetBinaryDumper;
typedef BinaryDumper<TrackerMuonDumperHelper> L1TrackerMuonBinaryDumper;
typedef BinaryDumper<CTL2EgammaDumperHelper> L1CTL2EgammaBinaryDumper;
//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PuppiBinaryDumper);
DEFINE_FWK_MODULE(L1JetBinaryDumper);
DEFINE_FWK_MODULE(L1TrackerMuonBinaryDumper);
DEFINE_FWK_MODULE(L1CTL2EgammaBinaryDumper);
