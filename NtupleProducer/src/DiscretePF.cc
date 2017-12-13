#include "../interface/DiscretePF.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Math/ProbFunc.h"
#include <TH1F.h>

namespace { 
    std::vector<float> vd2vf(const std::vector<double> & vd) {
        std::vector<float> ret;
        ret.insert(ret.end(), vd.begin(), vd.end());
        return ret;
    }
}
using namespace l1tpf_int;

const char * Region::inputTypeName(int type) {
    switch(InputType(type)) {
        case calo_type: return "Calo";
        case emcalo_type: return "EmCalo";
        case track_type: return "TK";
        case l1mu_type: return "Mu";
        case n_input_types: assert(false);
    }
    return "NO_SUCH_INPUT_TYPE";
}
const char * Region::outputTypeName(int type) {
    switch(OutputType(type)) {
        case any_type: return "";
        case charged_type: return "Charged";
        case neutral_type: return "Neutral";
        case electron_type: return "Electron";
        case pfmuon_type: return "Muon";
        case charged_hadron_type: return "ChargedHadron";
        case neutral_hadron_type: return "NeutralHadron";
        case photon_type: return "Photon";
        case n_output_types: assert(false);
    }
    return "NO_SUCH_OUTPUT_TYPE";
}

unsigned int Region::nInput(InputType type) const {
    switch(type) {
        case calo_type: return calo.size();
        case emcalo_type: return emcalo.size();
        case track_type: return track.size();
        case l1mu_type: return muon.size();
        case n_input_types: assert(false);
    }
    return 9999;
}
unsigned int Region::nOutput(OutputType type, bool usePuppi) const {
    unsigned int ret = 0;
    for (const auto &p : (usePuppi ? puppi : pf)) {
        if (p.hwPt <= 0) continue;
        if (!fiducialLocal(p.floatEta(), p.floatPhi())) continue;
        switch(type) {
            case any_type: ret++; break;
            case charged_type: if (p.intCharge() != 0) ret++; break;
            case neutral_type: if (p.intCharge() == 0) ret++; break;
            case electron_type: if (p.hwId == l1tpf::Particle::EL) ret++; break;
            case pfmuon_type: if (p.hwId == l1tpf::Particle::MU) ret++; break;
            case charged_hadron_type: if (p.hwId == l1tpf::Particle::CH) ret++; break;
            case neutral_hadron_type: if (p.hwId == l1tpf::Particle::NH) ret++; break;
            case photon_type: if (p.hwId == l1tpf::Particle::GAMMA) ret++; break;
            case n_output_types: assert(false);
        }
    }
    return ret;
}

RegionMapper::RegionMapper( const edm::ParameterSet& iConfig )  :
    useRelativeRegionalCoordinates_(false)
{
    if (iConfig.existsAs<std::vector<edm::ParameterSet>>("regions")) {
        useRelativeRegionalCoordinates_ = iConfig.getParameter<bool>("useRelativeRegionalCoordinates");
        for (const edm::ParameterSet & preg : iConfig.getParameter<std::vector<edm::ParameterSet>>("regions")) {
            std::vector<double> etaBoundaries = preg.getParameter<std::vector<double>>("etaBoundaries");
            unsigned int phiSlices = preg.getParameter<uint32_t>("phiSlices");
            float etaExtra = preg.getParameter<double>("etaExtra");
            float phiExtra = preg.getParameter<double>("phiExtra");
            float phiWidth = 2*M_PI/phiSlices;
            unsigned int ncalomax = 0, nemcalomax = 0, ntrackmax = 0, nmuonmax = 0, npfmax = 0, npuppimax = 0;
            if (preg.existsAs<uint32_t>("caloNMax")) ncalomax = preg.getParameter<uint32_t>("caloNMax");
            if (preg.existsAs<uint32_t>("emcaloNMax")) nemcalomax = preg.getParameter<uint32_t>("emcaloNMax");
            if (preg.existsAs<uint32_t>("trackNMax")) ntrackmax = preg.getParameter<uint32_t>("trackNMax");
            if (preg.existsAs<uint32_t>("muonNMax")) nmuonmax = preg.getParameter<uint32_t>("muonNMax");
            if (preg.existsAs<uint32_t>("pfNMax")) npfmax = preg.getParameter<uint32_t>("pfNMax");
            if (preg.existsAs<uint32_t>("puppiNMax")) npuppimax = preg.getParameter<uint32_t>("puppiNMax");
            for (unsigned int ieta = 0, neta = etaBoundaries.size()-1; ieta < neta; ++ieta) {
                for (unsigned int iphi = 0; iphi < phiSlices; ++iphi) {
                    float phiCenter = (iphi+0.5)*phiWidth-M_PI;
                    regions_.push_back(Region(
                            etaBoundaries[ieta], etaBoundaries[ieta+1], phiCenter, phiWidth, 
                            phiExtra, etaExtra, useRelativeRegionalCoordinates_,
                            ncalomax, nemcalomax, ntrackmax, nmuonmax, npfmax, npuppimax)); 
                }
            }
        }
        std::string trackRegionMode = "any";
        if (iConfig.existsAs<std::string>("trackRegionMode")) trackRegionMode = iConfig.getParameter<std::string>("trackRegionMode");
        if (trackRegionMode == "atVertex")    trackRegionMode_ = atVertex;
        else if (trackRegionMode == "atCalo") trackRegionMode_ = atCalo;
        else if (trackRegionMode == "any")    trackRegionMode_ = any;
        else throw cms::Exception("Configuration", "Unsupported value for trackRegionMode: " + trackRegionMode+" (allowed are 'atVertex', 'atCalo', 'any')");
        std::cout << "L1 RegionMapper: made " << regions_.size() << " regions" << std::endl;
    } else {
        // start off with a dummy region
        unsigned int ncalomax = 0, nemcalomax = 0, ntrackmax = 0, nmuonmax = 0, npfmax = 0, npuppimax = 0;
        regions_.push_back(Region(-5.5,5.5, 0,2*M_PI, 0.5, 0.5, useRelativeRegionalCoordinates_,
                                 ncalomax, nemcalomax, ntrackmax, nmuonmax, npfmax, npuppimax));
    }
}

void RegionMapper::addTrack( const l1tpf::Particle & t ) {
    // now let's be optimistic and make things very simple
    // we propagate in floating point the track to the calo
    // we add the track to the region corresponding to its vertex (eta,phi) coordinates AND its (eta,phi) calo coordinates
   for (Region &r : regions_) {
        bool inside = true;
        switch (trackRegionMode_) {
            case atVertex: inside = r.contains(t.eta(), t.phi()); break;
            case atCalo  : inside = r.contains(t.caloEta(), t.caloPhi()); break;
            case any     : inside = r.contains(t.eta(), t.phi()) || r.contains(t.caloEta(), t.caloPhi()); break;
        }
        if (inside) {
            PropagatedTrack prop;
            prop.fillInput(t.pt(), r.localEta(t.eta()), r.localPhi(t.phi()), t.charge(), t.dz(), 0);
            prop.fillPropagated(t.pt(), t.sigma(), t.caloSigma(), r.localEta(t.caloEta()), r.localPhi(t.caloPhi()), 0);
            float ndf = 2*t.quality()-4;
            prop.hwChi2  = round(t.normalizedChi2()*ndf*10);
            prop.hwStubs = round(t.quality());
            r.track.push_back(prop);
        }
    } 
}

void RegionMapper::addMuon( const l1tpf::Particle &mu ) {
    // now let's be optimistic and make things very simple
    // we don't propagate anything
    for (Region &r : regions_) {
        if (r.contains(mu.eta(), mu.phi())) {
            Muon prop;
            prop.fill(mu.pt(), r.localEta(mu.eta()), r.localPhi(mu.phi()), mu.charge(), mu.quality());
            r.muon.push_back(prop);
        }
    } 
}

void RegionMapper::addCalo( const l1tpf::Particle &p ) { 
    if (p.pt() == 0) return;
    for (Region &r : regions_) {
        if (r.contains(p.eta(), p.phi())) {
            CaloCluster calo;
            calo.fill(p.pt(), p.rawEmEt(), p.sigma(), r.localEta(p.eta()), r.localPhi(p.phi()), p.pdgId() == l1tpf::Particle::GAMMA, 0);
            r.calo.push_back(calo);
        }
    } 
}
void RegionMapper::addEmCalo( const l1tpf::Particle &p ) { 
    if (p.pt() == 0) return;
    for (Region &r : regions_) {
        if (r.contains(p.eta(), p.phi())) {
            CaloCluster calo;
            calo.fill(p.pt(), p.rawEmEt(), p.sigma(), r.localEta(p.eta()), r.localPhi(p.phi()), p.pdgId() == l1tpf::Particle::GAMMA, 0);
            r.emcalo.push_back(calo);
        }
    } 
}


std::vector<l1tpf::Particle> RegionMapper::fetch(bool puppi, float ptMin, bool discarded) const {
    std::vector<l1tpf::Particle> ret;
    for (const Region &r : regions_) {
        for (const PFParticle & p : (puppi ? r.puppi : (discarded ? r.pfdiscarded : r.pf ))) {
            if (regions_.size() > 1) {
                bool inside = true;
                switch (trackRegionMode_) {
                    case atVertex: inside = r.fiducialLocal(p.floatVtxEta(), p.floatVtxPhi()); break;
                    case atCalo  : inside = r.fiducialLocal(p.floatEta(), p.floatPhi()); break;
                    case any     : inside = r.fiducialLocal(p.floatVtxEta(), p.floatVtxPhi()); break; // WARNING: this may not be the best choice
                }
                if (!inside) continue;
            }
            if (p.floatPt() > ptMin) {
                ret.emplace_back( p.floatPt(), r.globalEta(p.floatVtxEta()), r.globalPhi(p.floatVtxPhi()), 0.13f, p.hwId, 0.f, p.floatDZ(), r.globalEta(p.floatEta()), r.globalPhi(p.floatPhi()), p.intCharge()  );
                ret.back().setStatus(p.hwStatus);
                ret.back().setPuppiWeight(p.floatPuppiW());
            }
        }
    }
    return ret;
}

std::vector<l1tpf::Particle> RegionMapper::fetchCalo(float ptMin, bool emcalo) const {
    std::vector<l1tpf::Particle> ret;
    for (const Region &r : regions_) {
        for (const CaloCluster & p : (emcalo ? r.emcalo : r.calo)) {
            if (regions_.size() > 1) {
                if (!r.fiducialLocal(p.floatEta(), p.floatPhi())) continue;
            }
            if (p.floatPt() > ptMin) {
                ret.emplace_back( p.floatPt(), r.globalEta(p.floatEta()), r.globalPhi(p.floatPhi()), 0.13f, (p.isEM || emcalo) ? l1tpf::Particle::GAMMA : l1tpf::Particle::NH );
            }
        }
    }
    return ret;
}

std::vector<l1tpf::Particle> RegionMapper::fetchTracks(float ptMin, bool fromPV) const {
    std::vector<l1tpf::Particle> ret;
    for (const Region &r : regions_) {
        for (const PropagatedTrack & p : r.track) {
            if (fromPV && !p.fromPV) continue;
            if (regions_.size() > 1) {
                if (!r.fiducialLocal(p.floatVtxEta(), p.floatVtxPhi())) continue;
            }
            if (p.floatPt() > ptMin) {
                ret.emplace_back( p.floatVtxPt(), r.globalEta(p.floatVtxEta()), r.globalPhi(p.floatVtxPhi()), 0.13f, p.muonLink ? l1tpf::Particle::MU : l1tpf::Particle::CH, 0.f, p.floatDZ(), r.globalEta(p.floatEta()), r.globalPhi(p.floatPhi()), p.intCharge() );
            }
        }
    }
    return ret;
}

std::pair<unsigned,unsigned> RegionMapper::totAndMaxInput(int type) const {
    unsigned ntot = 0, nmax = 0;
    for (const auto & r : regions_) {
        unsigned int ni = r.nInput(Region::InputType(type));
        ntot += ni;
        nmax = std::max(nmax, ni);
    }
    return std::make_pair(ntot,nmax);
}
std::pair<unsigned,unsigned> RegionMapper::totAndMaxOutput(int type, bool puppi) const {
    unsigned ntot = 0, nmax = 0;
    for (const auto & r : regions_) {
        unsigned int ni = r.nOutput(Region::OutputType(type),puppi);
        ntot += ni;
        nmax = std::max(nmax, ni);
    }
    return std::make_pair(ntot,nmax);
}




PFAlgo::PFAlgo( const edm::ParameterSet & iConfig ) :
    skipMuons_(iConfig.getParameter<bool>("metRate")),
    etaCharged_(iConfig.getParameter<double>("etaCharged")),
    puppiDr_(iConfig.getParameter<double>("puppiDr")),
    puppiEtaCuts_(vd2vf(iConfig.getParameter<std::vector<double>>("puppiEtaCuts"))),
    puppiPtCuts_(vd2vf(iConfig.getParameter<std::vector<double>>("puppiPtCuts"))),
    puppiPtCutsPhotons_(vd2vf(iConfig.getParameter<std::vector<double>>("puppiPtCutsPhotons"))),
    vtxRes_(iConfig.getParameter<double>("vtxRes")),
    vtxAdaptiveCut_(iConfig.getParameter<bool>("vtxAdaptiveCut")),
    drMatch_(0.2), ptMatchLow_(2.0), ptMatchHigh_(2.0), maxInvisiblePt_(20.0),
    useTrackCaloSigma_(false), rescaleUnmatchedTrack_(false),
    debug_(iConfig.getUntrackedParameter<int>("debug",0))
{
    edm::ParameterSet linkcfg = iConfig.getParameter<edm::ParameterSet>("linking");
    drMatch_ = linkcfg.getParameter<double>("trackCaloDR");
    ptMatchLow_ = linkcfg.getParameter<double>("trackCaloNSigmaLow");
    ptMatchHigh_ = linkcfg.getParameter<double>("trackCaloNSigmaHigh");
    useTrackCaloSigma_ = linkcfg.getParameter<bool>("useTrackCaloSigma");
    rescaleUnmatchedTrack_ = linkcfg.getParameter<bool>("rescaleUnmatchedTrack");
    maxInvisiblePt_  = linkcfg.getParameter<double>("maxInvisiblePt");
    if (rescaleUnmatchedTrack_) std::cout << "WARNING: rescaleUnmatchedTrack not yet implemented for integer code" << std::endl;

    intDrMatchBox_ = std::ceil(drMatch_ * CaloCluster::ETAPHI_SCALE * std::sqrt(M_PI/4));
    intPtMatchLowX4_ = std::ceil(ptMatchLow_ * 4);
    intPtMatchHighX4_ = std::ceil(ptMatchHigh_ * 4);
    intMaxInvisiblePt_ = std::round(maxInvisiblePt_ * CaloCluster::PT_SCALE);

    intDrMuonMatchBox_ = std::ceil(0.20 * CaloCluster::ETAPHI_SCALE * std::sqrt(M_PI/4));

    if (puppiEtaCuts_.size() != puppiPtCuts_.size() || puppiPtCuts_.size() != puppiPtCutsPhotons_.size()) {
        throw cms::Exception("Configuration", "Bad PUPPI config");
    }
    for (unsigned int i = 0, n = puppiEtaCuts_.size(); i < n; ++i) {
        intPuppiEtaCuts_.push_back( std::round(puppiEtaCuts_[i] * CaloCluster::ETAPHI_SCALE) );
        intPuppiPtCuts_.push_back( std::round(puppiPtCuts_[i] * CaloCluster::PT_SCALE) );
        intPuppiPtCutsPhotons_.push_back( std::round(puppiPtCutsPhotons_[i] * CaloCluster::PT_SCALE) );
    }
}

void PFAlgo::runPF(Region &r) const {
    initRegion(r);
    muonTrackLink(r);
    caloTrackLink(r);
    // then do muons
    if(!skipMuons_) { 
        for (int itk = 0, ntk = r.track.size(); itk < ntk; ++itk) {
            if (r.track[itk].muonLink) addTrackToPF(r, r.track[itk]);
        }
    }
}

void PFAlgo::initRegion(Region &r) const {
    r.inputSort();
    r.pf.clear(); r.puppi.clear();
    r.pfdiscarded.clear();
    for (auto & c : r.calo) c.used = false;
    for (auto & c : r.emcalo) c.used = false;
    for (auto & t : r.track) { t.used = false; t.muonLink = false; }
}

void PFAlgo::muonTrackLink(Region &r) const {
    // do a rectangular match for the moment; make a box of the same are as a 0.2 cone
    
    // first do the muon/tracking matching
    if (debug_) printf("INTMU Trying to link. I have %d tracks, %d muons\n", int(r.track.size()), int(r.muon.size()));
    for (int imu = 0, nmu = r.muon.size(); imu < nmu; ++imu) {
        if (debug_>1) printf("INTMU \t muon %d (pt %7.2f)\n", imu, r.muon[imu].floatPt());
        // FIXME the current pt matching is probably a bit tricky to do in integers
        float minPtDiff = 4; int imatch = -1;
        for (int itk = 0, ntk = r.track.size(); itk < ntk; ++itk) {
            if (debug_>1) printf("INTMU \t\t track %d (pt %7.2f): ", itk, r.track[itk].floatPt());
            if (std::abs(r.muon[imu].hwEta - r.track[itk].hwEta) >= intDrMuonMatchBox_) {
                if (debug_>1) printf("outside deta.\n");
                continue; }
            if (std::abs((r.muon[imu].hwPhi - r.track[itk].hwPhi) % CaloCluster::PHI_WRAP) >= intDrMuonMatchBox_) {
                if (debug_>1) printf("outside dphi.\n");  // phi wrapping would play nice if we had a power of 2 bins in phi
                continue; }
            // FIXME for the moment, we do the floating point matching in pt
            float dpt = (r.muon[imu].hwPt > r.track[itk].hwPt ? r.muon[imu].floatPt()/r.track[itk].floatPt() : r.track[itk].floatPt()/r.muon[imu].floatPt());
            if (dpt < minPtDiff) {
                if (debug_>1) printf("new best match.\n");
                if (imatch >= 0) r.track[imatch].muonLink = false;
                r.track[itk].muonLink = true;
                minPtDiff = dpt; imatch = itk;
            } else if (debug_>1) printf("new match, but %g worse than %d (%g).\n", dpt, imatch, minPtDiff);
        }
    }
}

void PFAlgo::caloTrackLink(Region &r) const {
    // then do track + calo linking
    if (debug_) printf("INT Trying to link. I have %d tracks, %d calo\n", int(r.track.size()), int(r.calo.size()));
    for (int itk = 0, ntk = r.track.size(); itk < ntk; ++itk) {
        if (r.track[itk].muonLink) continue;
        int iMatch = -1; int16_t dptMatch = 0;
        if (debug_>1) printf("INT \t track %d (pt %7.2f)\n", itk, r.track[itk].floatPt());
        for (int ic = 0, nc = r.calo.size(); ic < nc; ++ic) {
            int16_t hwCaloErr = useTrackCaloSigma_ ? r.track[itk].hwCaloPtErr : r.calo[ic].hwPtErr ;
            if (debug_>1) printf("INT \t\t calo %d (pt %7.2f): ", ic, r.calo[ic].floatPt());
            if (r.calo[ic].used) { 
                if (debug_>1) printf("used.\n"); 
                continue; }
            if (std::abs(r.calo[ic].hwEta - r.track[itk].hwEta) >= intDrMatchBox_)  { 
                if (debug_>1) printf("outside deta.\n"); 
                continue; }
            if (std::abs((r.calo[ic].hwPhi - r.track[itk].hwPhi) % CaloCluster::PHI_WRAP) >= intDrMatchBox_) {
                if (debug_>1) printf("outside dphi.\n"); 
                continue; } // phi wrapping would play nice if we had a power of 2 bins in phi
            if (r.calo[ic].hwPt + ((intPtMatchLowX4_ * hwCaloErr)/4) < r.track[itk].hwPt) { 
                if (debug_>1) printf("calo pt is too low.\n"); 
                continue; }
            int16_t dpt = std::abs(r.calo[ic].hwPt - r.track[itk].hwPt);
            if (iMatch == -1 || dpt < dptMatch) {
                if (debug_>1) printf("new best match.\n");
                dptMatch = dpt; iMatch = ic;
            } else if (debug_>1) printf("new match, but worse than %d.\n", iMatch);
        }
        if (iMatch > -1) {
            if (debug_) printf("INT Linking track of pt %7.2f eta %+5.2f phi %+5.2f (calo eta %+5.2f phi %+5.2f) to calo of pt %7.2f eta %+5.2f phi %+5.2f (dR = %.3f, deta = %.3f, dphi = %.3f)\n",
                            r.track[itk].floatVtxPt(), r.track[itk].floatVtxEta(), r.track[itk].floatVtxPhi(), r.track[itk].floatEta(), r.track[itk].floatPhi(),
                            r.calo[iMatch].floatPt(), r.calo[iMatch].floatEta(), r.calo[iMatch].floatPhi(),
                            ::deltaR(r.track[itk].floatEta(),r.track[itk].floatPhi(),r.calo[iMatch].floatEta(),r.calo[iMatch].floatPhi()), 
                            (r.track[itk].floatEta()-r.calo[iMatch].floatEta()), ::deltaPhi(r.track[itk].floatPhi(),r.calo[iMatch].floatPhi()));
            mergeTkCalo(r, r.track[itk], r.calo[iMatch]); 
        } else if (r.track[itk].hwPt < intMaxInvisiblePt_) {
            addTrackToPF(r, r.track[itk]);
        }
    }
    // then do all the un-linked calo (neutrals)
    for (int ic = 0, nc = r.calo.size(); ic < nc; ++ic) {
        if (!r.calo[ic].used) addCaloToPF(r, r.calo[ic]);
    }
    
}
PFParticle & PFAlgo::addTrackToPF(std::vector<PFParticle> &pfs, const PropagatedTrack &tk) const {
    PFParticle pf;
    pf.hwPt = tk.hwPt;
    pf.hwEta = tk.hwEta;
    pf.hwPhi = tk.hwPhi;
    pf.hwVtxEta = tk.hwEta; // FIXME: get from the track
    pf.hwVtxPhi = tk.hwPhi; // before propagation
    pf.track = tk;
    pf.cluster.hwPt = 0;
    pf.hwId = (tk.muonLink ? l1tpf::Particle::MU : l1tpf::Particle::CH);
    pf.hwStatus = 0;
    pfs.push_back(pf);
    return pfs.back();
}
PFParticle & PFAlgo::addCaloToPF(std::vector<PFParticle> &pfs, const CaloCluster &calo) const {
    PFParticle pf;
    pf.hwPt = calo.hwPt;
    pf.hwEta = calo.hwEta;
    pf.hwPhi = calo.hwPhi;
    pf.hwVtxEta = calo.hwEta; 
    pf.hwVtxPhi = calo.hwPhi; 
    pf.track.hwPt = 0;
    pf.cluster = calo;
    pf.hwId = (calo.isEM ? l1tpf::Particle::GAMMA : l1tpf::Particle::NH);
    pf.hwStatus = 0;
    pfs.push_back(pf);
    return pfs.back();
}

void PFAlgo::mergeTkCalo(Region &r, const PropagatedTrack &tk, CaloCluster & calo) const {
    int16_t hwCaloErr = useTrackCaloSigma_ ? tk.hwCaloPtErr : calo.hwPtErr;
    int16_t totErr2 = (intPtMatchHighX4_*(tk.hwPtErr + hwCaloErr))/4; // should be in quadrature
    if (calo.hwPt - tk.hwPt > totErr2) { // calo > track + 2 sigma
        calo.hwPt -= tk.hwPt;
        if (debug_) printf("INT   case 1: add track, reduce calo pt to %7.2f\n", calo.floatPt());
        addTrackToPF(r, tk);
    } else { // |calo - track|  < 2 sigma
        // FIXME for the moment, doing the full weighted average
        // in float; to test alternatives in the future.
        float wcalo = 0, wtk = 1;
        int16_t ptAvg = tk.hwPt;
        if (tk.hwPtErr > 0) { // it can easily be zero due to underflows
            wcalo = 1.0/std::pow(float(hwCaloErr),2);
            wtk   = 1.0/std::pow(float(tk.hwPtErr),2);
            ptAvg = std::round( ( float(calo.hwPt)*wcalo + float(tk.hwPt)*wtk ) / (wcalo + wtk) );
        }
        PFParticle & pf = addTrackToPF(r, tk);
        pf.hwPt = ptAvg;
        pf.cluster = calo;
        pf.hwId = (calo.isEM ? l1tpf::Particle::EL : l1tpf::Particle::CH);
        calo.used = true; // used, not to be added to the PF again
        if (debug_) printf("INT   case 2: merge, avg pt %7.2f from calo (%7.2f +- %.2f), track (%7.2f +- %.2f) weights %.3f %.3f\n", 
                                pf.floatPt(), calo.floatPt(), calo.floatPtErr(), tk.floatPt(), tk.floatPtErr(), wcalo/(wcalo+wtk), wtk/(wcalo+wtk));
    } 
}

void PFAlgo::runPuppi(Region &r, float npu, float alphaCMed, float alphaCRms, float alphaFMed, float alphaFRms) const {
    computePuppiWeights(r, alphaCMed, alphaCRms, alphaFMed, alphaFRms);
    fillPuppi(r);
}

void PFAlgo::runChargedPV(Region &r, float z0) const {
    int16_t iZ0 = round(z0 * InputTrack::Z0_SCALE);
    int16_t iDZ  = round(1.5 * vtxRes_ * InputTrack::Z0_SCALE);
    int16_t iDZ2 = vtxAdaptiveCut_ ? round(4.0 * vtxRes_ * InputTrack::Z0_SCALE) : iDZ;
    for (PFParticle & p : r.pf) {
        bool barrel = std::abs(p.track.hwVtxEta) < InputTrack::VTX_ETA_1p3;
        if (r.relativeCoordinates) barrel = (std::abs(r.globalAbsEta(p.track.floatVtxEta())) < 1.3); // FIXME could make a better integer implementation
        p.chargedPV = (p.hwId <= 1 && std::abs(p.track.hwZ0 - iZ0) < (barrel ? iDZ : iDZ2));
    }
}

void PFAlgo::computePuppiWeights(Region &r, float alphaCMed, float alphaCRms, float alphaFMed, float alphaFRms) const {
    int16_t ietacut = std::round(etaCharged_ * CaloCluster::ETAPHI_SCALE);
    // FIXME floats for now
    float puppiDr2 = std::pow(puppiDr_,2);
    for (PFParticle & p : r.pf) {
        // charged
        if (p.hwId <= 1) {
            p.setPuppiW(p.chargedPV ? 1.0 : 0); 
            if (debug_) printf("PUPPI \t charged id %1d pt %7.2f eta %+5.2f phi %+5.2f  alpha %+6.2f x2 %+6.2f --> puppi weight %.3f   puppi pt %7.2f \n", p.hwId, p.floatPt(), p.floatEta(), p.floatPhi(), 0., 0., p.floatPuppiW(), p.floatPt()*p.floatPuppiW());
            continue;
        }
        // neutral
        float alphaC = 0, alphaF = 0;
        for (const PFParticle & p2 : r.pf) {
            float dr2 = ::deltaR2(p.floatEta(), p.floatPhi(), p2.floatEta(), p2.floatPhi());
            if (dr2 > 0 && dr2 < puppiDr2) {
                float w = std::pow(p2.floatPt(),2) / dr2;
                alphaF += w;
                if (p2.chargedPV) alphaC += w;
            }
        }
        float alpha = -99, x2 = -99;
        bool central = std::abs(p.hwEta) < ietacut;
        if (r.relativeCoordinates) central = (std::abs(r.globalAbsEta(p.floatEta())) < etaCharged_); // FIXME could make a better integer implementation
        if (central) {
            if (alphaC > 0) {
                alpha = std::log(alphaC);
                x2 = (alpha - alphaCMed) * std::abs(alpha - alphaCMed) / std::pow(alphaCRms,2);
                p.setPuppiW( ROOT::Math::chisquared_cdf(x2,1) );
            } else {
                p.setPuppiW(0);
            }
        } else {
            if (alphaF > 0) {
                alpha = std::log(alphaF);
                x2 = (alpha - alphaFMed) * std::abs(alpha - alphaFMed) / std::pow(alphaFRms,2);
                p.setPuppiW( ROOT::Math::chisquared_cdf(x2,1) );
            } else {
                p.setPuppiW(0);
            }
        }
        if (debug_) printf("PUPPI \t neutral id %1d pt %7.2f eta %+5.2f phi %+5.2f  alpha %+6.2f x2 %+7.2f --> puppi weight %.3f   puppi pt %7.2f \n", p.hwId, p.floatPt(), p.floatEta(), p.floatPhi(), alpha, x2, p.floatPuppiW(), p.floatPt()*p.floatPuppiW());
    }
}

void PFAlgo::doVertexing(std::vector<Region> &rs, VertexAlgo algo, float &pvdz) const {
    int lNBins = int(40./vtxRes_);
    if (algo == TPVtxAlgo) lNBins *= 3;
    std::unique_ptr<TH1F> h_dz(new TH1F("h_dz","h_dz",lNBins,-20,20));
    for (const Region & r : rs) {
        for (const PropagatedTrack & p : r.track) {
            if (rs.size() > 1) {
                if (!r.fiducialLocal(p.floatVtxEta(), p.floatVtxPhi())) continue; // skip duplicates
            }
            h_dz->Fill( p.floatDZ(), std::min(p.floatPt(), 50.f) );
        }
    }
    switch(algo) {
        case OldVtxAlgo: {
                             int imaxbin = h_dz->GetMaximumBin();
                             pvdz = h_dz->GetXaxis()->GetBinCenter(imaxbin);
                         }; break;
        case TPVtxAlgo: {
                            float max = 0; int bmax = -1;
                            for (int b = 1; b <= lNBins; ++b) {
                                float sum3 = h_dz->GetBinContent(b) + h_dz->GetBinContent(b+1) + h_dz->GetBinContent(b-1);
                                if (bmax == -1 || sum3 > max) { max = sum3; bmax = b; }
                            }
                            pvdz = h_dz->GetXaxis()->GetBinCenter(bmax); 
                        }; break;
    }
    int16_t iZ0 = round(pvdz * InputTrack::Z0_SCALE);
    int16_t iDZ  = round(1.5 * vtxRes_ * InputTrack::Z0_SCALE);
    int16_t iDZ2 = vtxAdaptiveCut_ ? round(4.0 * vtxRes_ * InputTrack::Z0_SCALE) : iDZ;
    for (Region & r : rs) {
        for (PropagatedTrack & p : r.track) {
            bool central = std::abs(p.hwVtxEta) < InputTrack::VTX_ETA_1p3;
            if (r.relativeCoordinates) central = (std::abs(r.globalAbsEta(p.floatVtxEta())) < 1.3); // FIXME could make a better integer implementation
            p.fromPV = (std::abs(p.hwZ0 - iZ0) < (central ? iDZ : iDZ2));
        }
    }

}

void PFAlgo::computePuppiMedRMS(const std::vector<Region> &rs, float &alphaCMed, float &alphaCRms, float &alphaFMed, float &alphaFRms) const {
    std::vector<float> alphaFs;
    std::vector<float> alphaCs;
    int16_t ietacut = std::round(etaCharged_ * CaloCluster::ETAPHI_SCALE);
    float puppiDr2 = std::pow(puppiDr_,2);
    for (const Region & r : rs) {
        for (const PFParticle & p : r.pf) {
            bool central = std::abs(p.hwEta) < ietacut;
            if (r.relativeCoordinates) central = (r.globalAbsEta(p.floatEta()) < etaCharged_); // FIXME could make a better integer implementation
            if (central) {
                if (p.hwId > 1 || p.chargedPV) continue;
            }
            float alphaC = 0, alphaF = 0;
            for (const PFParticle & p2 : r.pf) {
                float dr2 = ::deltaR2(p.floatEta(), p.floatPhi(), p2.floatEta(), p2.floatPhi());
                if (dr2 > 0 && dr2 < puppiDr2) {
                    float w = std::pow(p2.floatPt(),2) / std::max<float>(0.01f, dr2);
                    alphaF += w;
                    if (p2.chargedPV) alphaC += w;
                }
            }
            if (central) {
                if (alphaC > 0) alphaCs.push_back(std::log(alphaC));
            } else {
                if (alphaF > 0) alphaFs.push_back(std::log(alphaF));
            }
        }
    }
  std::sort(alphaCs.begin(),alphaCs.end());
  std::sort(alphaFs.begin(),alphaFs.end());

  if (alphaCs.size() > 1){
      alphaCMed = alphaCs[alphaCs.size()/2+1];
      double sum = 0.0;
      for (float alpha : alphaCs) sum += std::pow(alpha-alphaCMed,2);
      alphaCRms = std::sqrt(float(sum)/alphaCs.size());
  } else {
      alphaCMed = 8.; alphaCRms = 8.;
  }

  if (alphaFs.size() > 1){
      alphaFMed = alphaFs[alphaFs.size()/2+1];
      double sum = 0.0;
      for (float alpha : alphaFs) sum += std::pow(alpha-alphaFMed,2);
      alphaFRms = std::sqrt(float(sum)/alphaFs.size());
  } else {
      alphaFMed = 6.; alphaFRms = 6.;
  }
  if (debug_) printf("PUPPI \t alphaC = %+6.2f +- %6.2f (%4lu), alphaF = %+6.2f +- %6.2f (%4lu)\n", alphaCMed, alphaCRms, alphaCs.size(), alphaFMed, alphaFRms, alphaFs.size());
}

void PFAlgo::fillPuppi(Region &r) const {
    constexpr uint16_t PUPPIW_0p01 = std::round(0.01 * PFParticle::PUPPI_SCALE);
    r.puppi.clear();
    for (PFParticle & p : r.pf) {
        if (p.hwId == l1tpf::Particle::MU) {
            r.puppi.push_back(p);
        } else if (p.hwId <= 1) { // charged
            if (p.hwPuppiWeight > 0) {
                r.puppi.push_back(p);
            }
        } else { // neutral
            if (p.hwPuppiWeight > PUPPIW_0p01) {
                // FIXME would work better with PUPPI_SCALE being a power of two, to do the shift
                // FIXME done with floats
                int16_t hwPt = ( float(p.hwPt) * float(p.hwPuppiWeight) / float(PFParticle::PUPPI_SCALE) );
                int16_t hwPtCut = 0, hwAbsEta = r.relativeCoordinates ? round(r.globalAbsEta(p.floatEta()) * CaloCluster::ETAPHI_SCALE) : std::abs(p.hwEta);
                for (unsigned int ietaBin = 0, nBins = intPuppiEtaCuts_.size(); ietaBin < nBins; ++ietaBin) {
                    if (hwAbsEta < intPuppiEtaCuts_[ietaBin]) {
                        hwPtCut = (p.hwId == l1tpf::Particle::GAMMA ? intPuppiPtCutsPhotons_[ietaBin] : intPuppiPtCuts_[ietaBin]);
                        break;
                    }
                }
                if (hwPt > hwPtCut) {
                    r.puppi.push_back(p);
                    r.puppi.back().hwPt = hwPt;
                }
            }
        }
    }
}

