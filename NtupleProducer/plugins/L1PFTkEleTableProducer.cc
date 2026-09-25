// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/L1TCorrelator/interface/TkElectron.h"
#include "DataFormats/L1TCorrelator/interface/TkElectronFwd.h"
#include "DataFormats/L1TParticleFlow/interface/egamma.h"
#include "DataFormats/L1TParticleFlow/interface/PFTrack.h"
#include "DataFormats/L1TParticleFlow/interface/PFCluster.h"

#include <algorithm>

class L1PFTkEleTableProducer : public edm::global::EDProducer<> {
public:
    explicit L1PFTkEleTableProducer(const edm::ParameterSet&);
    ~L1PFTkEleTableProducer() override = default;

private:
    void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

    std::string name_;
    edm::EDGetTokenT<edm::View<l1t::TkElectron>> src_;
    StringCutObjectSelector<l1t::TkElectron> sel_;

    // decoded track/calo collections to cross-reference against, labelled by region (e.g. "HGCal", "Barrel")
    std::vector<std::pair<std::string, edm::EDGetTokenT<l1t::PFTrackCollection>>> decTkColls_;
    std::vector<std::pair<std::string, edm::EDGetTokenT<l1t::PFClusterCollection>>> decCaloColls_;
};

namespace {
    // suffix of the decoded table name, e.g. "l1tLayer1HGCal:DecodedTK" -> "HGCal"
    std::string decRegionSuffix(const edm::InputTag& tag) {
        const std::string prefix = "l1tLayer1";
        const std::string& label = tag.label();
        if (label.compare(0, prefix.size(), prefix) == 0)
            return label.substr(prefix.size());
        return label;
    }
}  // namespace

L1PFTkEleTableProducer::L1PFTkEleTableProducer(const edm::ParameterSet& iConfig)
    : name_(iConfig.getParameter<std::string>("name")),
      src_(consumes<edm::View<l1t::TkElectron>>(iConfig.getParameter<edm::InputTag>("src"))),
      sel_(iConfig.getParameter<std::string>("cut"), true) {
    for (const auto& tag : iConfig.getParameter<std::vector<edm::InputTag>>("decTkRegions"))
        decTkColls_.emplace_back(decRegionSuffix(tag), consumes<l1t::PFTrackCollection>(tag));
    for (const auto& tag : iConfig.getParameter<std::vector<edm::InputTag>>("decCaloRegions"))
        decCaloColls_.emplace_back(decRegionSuffix(tag), consumes<l1t::PFClusterCollection>(tag));
    produces<nanoaod::FlatTable>();
}

void L1PFTkEleTableProducer::produce(edm::StreamID id,
                                     edm::Event& iEvent,
                                     const edm::EventSetup& iSetup) const {
    edm::Handle<edm::View<l1t::TkElectron>> src;
    iEvent.getByToken(src_, src);

    std::vector<const l1t::TkElectron*> selected;
    for (const auto& ele : *src)
        if (sel_(ele))
            selected.push_back(&ele);

    unsigned int ncands = selected.size();
    auto out = std::make_unique<nanoaod::FlatTable>(ncands, name_, false, true);

    std::vector<float> hwPt(ncands, 0), hwEta(ncands, 0), hwPhi(ncands, 0),
        hwQual(ncands, 0), hwIso(ncands, 0),
        hwDEta(ncands, 0), hwDPhi(ncands, 0), hwZ0(ncands, 0),
        hwCharge(ncands, 0), hwIDScore(ncands, 0),
        hwTkRedChi2RPhi(ncands, 0), hwTkCaloDphi(ncands, 0),
        hwCaloShowerShape(ncands, 0), hwCaloTkPtRatio(ncands, 0);

    for (unsigned int i = 0; i < ncands; ++i) {
        const auto& ele = *selected[i];
        if (ele.encoding() != l1t::TkEm::HWEncoding::CT)
            continue;

        auto word = ele.egBinaryWord<l1ct::EGIsoEleObj::BITWIDTH>();
        l1ct::EGIsoEleObj ct;
        ct.initFromBits(word);

        hwPt[i]              = ct.hwPt.to_float();
        hwEta[i]             = ct.hwEta.to_float();
        hwPhi[i]             = ct.hwPhi.to_float();
        hwQual[i]            = float(ct.hwQual);
        hwIso[i]             = ct.hwIso.to_float();
        hwDEta[i]            = ct.hwDEta.to_float();
        hwDPhi[i]            = ct.hwDPhi.to_float();
        hwZ0[i]              = ct.hwZ0.to_float();
        hwCharge[i]          = float(ct.hwCharge);
        hwIDScore[i]         = ct.hwIDScore.to_float();
        hwTkRedChi2RPhi[i]   = ct.hwTkRedChi2RPhi.to_float();
        hwTkCaloDphi[i]      = ct.hwTkCaloDphi.to_float();
        hwCaloShowerShape[i] = ct.hwCaloShowerShape.to_float();
        hwCaloTkPtRatio[i]   = ct.hwCaloTkPtRatio.to_float();
    }

    out->addColumn<float>("hwPt",              hwPt,              "CT hardware pT [LSB = 0.25 GeV]");
    out->addColumn<float>("hwEta",             hwEta,             "CT hardware eta at calo face [LSB = pi/720]");
    out->addColumn<float>("hwPhi",             hwPhi,             "CT hardware phi at calo face [LSB = pi/720]");
    out->addColumn<float>("hwQual",            hwQual,            "CT hardware quality flags");
    out->addColumn<float>("hwIso",             hwIso,             "CT hardware isolation pT [LSB = 0.25 GeV]");
    out->addColumn<float>("hwDEta",            hwDEta,            "CT hardware track-calo delta eta [LSB = pi/720]");
    out->addColumn<float>("hwDPhi",            hwDPhi,            "CT hardware track-calo delta phi [LSB = pi/720]");
    out->addColumn<float>("hwZ0",              hwZ0,              "CT hardware track z0 [LSB = 0.05 cm]");
    out->addColumn<float>("hwCharge",          hwCharge,          "CT hardware charge (1=positive, 0=negative)");
    out->addColumn<float>("hwIDScore",         hwIDScore,         "CT hardware ID score [-1, 1]");
    out->addColumn<float>("hwTkRedChi2RPhi",   hwTkRedChi2RPhi,   "CT hardware track reduced chi2 in R-phi (4-bit bin)");
    out->addColumn<float>("hwTkCaloDphi",      hwTkCaloDphi,      "CT hardware track-calo delta phi (7 bits)");
    out->addColumn<float>("hwCaloShowerShape", hwCaloShowerShape, "CT hardware calo shower shape (6 bits)");
    out->addColumn<float>("hwCaloTkPtRatio",   hwCaloTkPtRatio,   "CT hardware calo/track pT ratio (10 bits)");

    // cross-references to the decoded track/cluster tables, matched by pointer identity
    for (unsigned int k = 0; k < decTkColls_.size(); ++k) {
        edm::Handle<l1t::PFTrackCollection> decTks;
        if (!iEvent.getByToken(decTkColls_[k].second, decTks))
            continue;
        std::vector<int> idx(ncands, -1);
        for (unsigned int i = 0; i < ncands; ++i) {
            const auto& ele = *selected[i];
            if (ele.trkPtr().isNull())
                continue;
            const auto* trk = ele.trkPtr().get();
            for (unsigned int j = 0; j < decTks->size(); ++j) {
                if ((*decTks)[j].track().isNonnull() && (*decTks)[j].track().get() == trk) {
                    idx[i] = j;
                    break;
                }
            }
        }
        out->addColumn<int>("decTk" + decTkColls_[k].first + "Idx",
                            idx,
                            "index in the DecTk" + decTkColls_[k].first +
                                " table of the decoded track used to build this electron (-1 if not found)");
    }
    for (unsigned int k = 0; k < decCaloColls_.size(); ++k) {
        edm::Handle<l1t::PFClusterCollection> decClus;
        if (!iEvent.getByToken(decCaloColls_[k].second, decClus))
            continue;
        std::vector<int> idx(ncands, -1);
        for (unsigned int i = 0; i < ncands; ++i) {
            const auto& ele = *selected[i];
            if (ele.egCaloPtr().isNull())
                continue;
            const auto* clu = ele.egCaloPtr().get();
            for (unsigned int j = 0; j < decClus->size(); ++j) {
                const auto& cons = (*decClus)[j].constituentsAndFractions();
                if (!cons.empty() && cons.front().first.isNonnull() && cons.front().first.get() == clu) {
                    idx[i] = j;
                    break;
                }
            }
        }
        out->addColumn<int>("decEmCalo" + decCaloColls_[k].first + "Idx",
                            idx,
                            "index in the DecEmCalo" + decCaloColls_[k].first +
                                " table of the decoded cluster used to build this electron (-1 if not found)");
    }

    iEvent.put(std::move(out));
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PFTkEleTableProducer);
