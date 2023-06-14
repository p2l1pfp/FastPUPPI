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

#include <cstdio>
#include <cstdint>
#include <fstream>


class PuppiDumperHelper {
    public:
        PuppiDumperHelper(const edm::ParameterSet &cfg, edm::ConsumesCollector cc) :
            src_(cc.consumes<edm::View<l1t::PFCandidate>>(cfg.getParameter<edm::InputTag>("src"))) {}
        void dump(const edm::Event & iEvent, std::fstream & out) {
            edm::Handle<edm::View<l1t::PFCandidate>> src;
            iEvent.getByToken(src_, src);
            uint64_t nobj = src->size();
            out.write(reinterpret_cast<const char*>(&nobj), sizeof(uint64_t));
            for (auto & c : *src) {
                uint64_t packed = c.encodedPuppi64();
                out.write(reinterpret_cast<const char*>(&packed), sizeof(uint64_t));
            }
        }
    private:
        edm::EDGetTokenT<edm::View<l1t::PFCandidate>> src_;
};

class JetDumperHelper {
    public:
        JetDumperHelper(const edm::ParameterSet &cfg, edm::ConsumesCollector cc) :
            src_(cc.consumes<edm::View<l1t::PFJet>>(cfg.getParameter<edm::InputTag>("src"))),
            ptMin_(cfg.getParameter<double>("ptMin")) {}
        void dump(const edm::Event & iEvent, std::fstream & out) {
            edm::Handle<edm::View<l1t::PFJet>> src;
            iEvent.getByToken(src_, src);
            std::vector<uint64_t> data(1, 0u); // leave one empty word at the beginning
            for (l1t::PFJet c : *src) {
                if (c.pt() > ptMin_) {
                    const std::array<uint64_t, 2> & enc = c.encodedJet();
                    data.push_back(enc[0]);
                    data.push_back(enc[1]);
                }
            }
            // now fill the size
            data[0] = (data.size()-1)/2;
            out.write(reinterpret_cast<const char*>(&data[0]), data.size()*sizeof(uint64_t));
        }
    private:
        edm::EDGetTokenT<edm::View<l1t::PFJet>> src_;
        float ptMin_;
};

template<typename Helper>
class BinaryDumper: public edm::one::EDAnalyzer<>  {
    public:
        explicit BinaryDumper(const edm::ParameterSet &iConfig) :
            helper_(iConfig, consumesCollector()),
            outName_(iConfig.getParameter<std::string>("outName")),
            out_() {}
    private:
        Helper helper_;
        std::string outName_;
        std::fstream out_;

        void beginJob() override {
            out_.open(outName_, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        }
        void analyze(const edm::Event & iEvent, const edm::EventSetup&) override {
            helper_.dump(iEvent, out_);
        }
        void endJob() override {
            out_.close();
        }

};

typedef BinaryDumper<PuppiDumperHelper> L1PuppiBinaryDumper;
typedef BinaryDumper<JetDumperHelper> L1JetBinaryDumper;
//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PuppiBinaryDumper);
DEFINE_FWK_MODULE(L1JetBinaryDumper);
