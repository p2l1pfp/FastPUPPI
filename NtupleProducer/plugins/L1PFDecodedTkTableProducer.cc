// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/L1TParticleFlow/interface/PFTrack.h"

#include <algorithm>

class L1PFDecodedTkTableProducer : public edm::global::EDProducer<> {
public:
    explicit L1PFDecodedTkTableProducer(const edm::ParameterSet&);
    ~L1PFDecodedTkTableProducer() override = default;

private:
    void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

    std::string name_;
    edm::EDGetTokenT<l1t::PFTrackCollection> tracks_;
    StringCutObjectSelector<l1t::PFTrack> sel_;
};

L1PFDecodedTkTableProducer::L1PFDecodedTkTableProducer(const edm::ParameterSet& iConfig)
    : name_(iConfig.getParameter<std::string>("name")),
      tracks_(consumes<l1t::PFTrackCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      sel_(iConfig.getParameter<std::string>("cut"), true) {
    produces<nanoaod::FlatTable>();
}

void L1PFDecodedTkTableProducer::produce(edm::StreamID id,
                                         edm::Event& iEvent,
                                         const edm::EventSetup& iSetup) const {
    edm::Handle<l1t::PFTrackCollection> tracks;
    iEvent.getByToken(tracks_, tracks);

    std::vector<const l1t::PFTrack*> selected;
    for (const auto& track : *tracks)
        if (sel_(track))
            selected.push_back(&track);

    // create the table
    unsigned int ncands = selected.size();
    auto out = std::make_unique<nanoaod::FlatTable>(ncands, name_, false, true);

    std::vector<float> hwPt(ncands, 0), hwEta(ncands, 0), hwPhi(ncands, 0),
        hwQual(ncands, 0), charge(ncands, 0);

    std::vector<float> twValid(ncands, 0),
        twRinvBits(ncands, 0), twPhi0Bits(ncands, 0), twTanlBits(ncands, 0),
        twZ0Bits(ncands, 0), twD0Bits(ncands, 0),
        twChi2RPhiBits(ncands, 0), twChi2RZBits(ncands, 0), twBendChi2Bits(ncands, 0),
        twHitPatternBits(ncands, 0), twMVAQualityBits(ncands, 0), twMVAOtherBits(ncands, 0),
        twRinv(ncands, 0), twPhi0(ncands, 0), twTanl(ncands, 0),
        twZ0(ncands, 0), twD0(ncands, 0), twMVAOther(ncands, 0);

    std::vector<float> srcPt(ncands, 0), srcEta(ncands, 0), srcPhi(ncands, 0),
        srcZ0(ncands, 0), srcD0(ncands, 0),
        srcChi2Red(ncands, 0), srcChi2ZRed(ncands, 0), srcBendChi2(ncands, 0),
        srcNStubs(ncands, 0);

    for (unsigned int i = 0; i < ncands; ++i) {
        const l1t::PFTrack* cand = selected[i];

        // hardware values of the track decoded by the regionizer and stored in the PFTrack
        hwPt[i]   = cand->hwPt();
        hwEta[i]  = cand->hwEta();
        hwPhi[i]  = cand->hwPhi();
        hwQual[i] = cand->hwQual();
        charge[i] = cand->charge();

        // raw and dequantized fields of the 96-bit track word of the source TTTrack
        const TTTrack_TrackWord& tw = cand->trackWord();
        twValid[i]         = tw.getValid();
        twRinvBits[i]      = tw.getRinvBits();
        twPhi0Bits[i]      = tw.getPhiBits();
        twTanlBits[i]      = tw.getTanlBits();
        twZ0Bits[i]        = tw.getZ0Bits();
        twD0Bits[i]        = tw.getD0Bits();
        twChi2RPhiBits[i]  = tw.getChi2RPhiBits();
        twChi2RZBits[i]    = tw.getChi2RZBits();
        twBendChi2Bits[i]  = tw.getBendChi2Bits();
        twHitPatternBits[i]= tw.getHitPatternBits();
        twMVAQualityBits[i]= tw.getMVAQualityBits();
        twMVAOtherBits[i]  = tw.getMVAOtherBits();

        twRinv[i]    = tw.getRinv();
        twPhi0[i]    = tw.getPhi();
        twTanl[i]    = tw.getTanl();
        twZ0[i]      = tw.getZ0();
        twD0[i]      = tw.getD0();
        twMVAOther[i]= tw.getMVAOther();

        // full-precision values of the source TTTrack
        if (cand->track().isNonnull()) {
            srcPt[i]      = cand->track()->momentum().perp();
            srcEta[i]     = cand->track()->momentum().eta();
            srcPhi[i]     = cand->track()->momentum().phi();
            srcZ0[i]      = cand->track()->z0();
            srcD0[i]      = cand->track()->d0();
            srcChi2Red[i] = cand->track()->chi2Red();
            srcChi2ZRed[i]= cand->track()->chi2ZRed();
            srcBendChi2[i]= cand->track()->chi2Bend();
            srcNStubs[i]  = cand->track()->getStubRefs().size();
        }
    }

    out->addColumn<float>("hwPt",   hwPt,   "decoded regionizer track pT as integer [LSB = 0.25 GeV]");
    out->addColumn<float>("hwEta",  hwEta,  "decoded regionizer track eta at calo face, relative to region center [LSB = pi/720 rad]");
    out->addColumn<float>("hwPhi",  hwPhi,  "decoded regionizer track phi, relative to region center [LSB = pi/720 rad]");
    out->addColumn<float>("hwQual", hwQual, "decoded regionizer track quality flags");
    out->addColumn<float>("charge", charge, "track charge (+1/-1)");

    out->addColumn<float>("twValid",          twValid,          "valid bit of the source track word");
    out->addColumn<float>("twRinvBits",       twRinvBits,       "raw Rinv field of the source track word (15 bits)");
    out->addColumn<float>("twPhi0Bits",       twPhi0Bits,       "raw phi0 field of the source track word (12 bits)");
    out->addColumn<float>("twTanlBits",       twTanlBits,       "raw tan(lambda) field of the source track word (16 bits)");
    out->addColumn<float>("twZ0Bits",         twZ0Bits,         "raw z0 field of the source track word (12 bits)");
    out->addColumn<float>("twD0Bits",         twD0Bits,         "raw d0 field of the source track word (13 bits)");
    out->addColumn<float>("twChi2RPhiBits",   twChi2RPhiBits,   "raw chi2/dof r-phi field of the source track word (4 bits)");
    out->addColumn<float>("twChi2RZBits",     twChi2RZBits,     "raw chi2/dof r-z field of the source track word (4 bits)");
    out->addColumn<float>("twBendChi2Bits",   twBendChi2Bits,   "raw bend-chi2 field of the source track word (3 bits)");
    out->addColumn<float>("twHitPatternBits", twHitPatternBits, "raw stub hit pattern of the source track word (7 bits)");
    out->addColumn<float>("twMVAQualityBits", twMVAQualityBits, "raw track quality MVA field of the source track word (3 bits)");
    out->addColumn<float>("twMVAOtherBits",   twMVAOtherBits,   "raw specialized MVA field of the source track word (6 bits)");

    out->addColumn<float>("twRinv",     twRinv,     "dequantized Rinv from the source track word [cm^-1]");
    out->addColumn<float>("twPhi0",     twPhi0,     "dequantized sector-local phi0 from the source track word [rad]");
    out->addColumn<float>("twTanl",     twTanl,     "dequantized tan(lambda) from the source track word");
    out->addColumn<float>("twZ0",       twZ0,       "dequantized z0 from the source track word [cm]");
    out->addColumn<float>("twD0",       twD0,       "dequantized d0 from the source track word [cm]");
    out->addColumn<float>("twMVAOther", twMVAOther, "specialized MVA bits from the source track word");

    out->addColumn<float>("srcPt",       srcPt,       "full-precision pT of the source TTTrack [GeV]");
    out->addColumn<float>("srcEta",      srcEta,      "full-precision eta of the source TTTrack");
    out->addColumn<float>("srcPhi",      srcPhi,      "full-precision phi of the source TTTrack");
    out->addColumn<float>("srcZ0",       srcZ0,       "full-precision z0 of the source TTTrack [cm]");
    out->addColumn<float>("srcD0",       srcD0,       "full-precision d0 of the source TTTrack [cm]");
    out->addColumn<float>("srcChi2Red",  srcChi2Red,  "full-precision reduced chi2 r-phi of the source TTTrack");
    out->addColumn<float>("srcChi2ZRed", srcChi2ZRed, "full-precision reduced chi2 r-z of the source TTTrack");
    out->addColumn<float>("srcBendChi2", srcBendChi2, "full-precision bend-chi2 of the source TTTrack");
    out->addColumn<float>("srcNStubs",   srcNStubs,   "number of stubs of the source TTTrack");

    // save to the event branches
    iEvent.put(std::move(out));
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PFDecodedTkTableProducer);
