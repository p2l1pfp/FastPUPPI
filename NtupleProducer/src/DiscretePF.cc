#include "../interface/DiscretePF.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Math/ProbFunc.h"

using namespace l1tpf_int;

RegionMapper::RegionMapper( const edm::ParameterSet& iConfig) 
{
    if (iConfig.existsAs<std::vector<edm::ParameterSet>>("regions")) {
        for (const edm::ParameterSet & preg : iConfig.getParameter<std::vector<edm::ParameterSet>>("regions")) {
            std::vector<double> etaBoundaries = preg.getParameter<std::vector<double>>("etaBoundaries");
            unsigned int phiSlices = preg.getParameter<uint32_t>("phiSlices");
            float etaExtra = preg.getParameter<double>("etaExtra");
            float phiExtra = preg.getParameter<double>("phiExtra");
            float phiWidth = 2*M_PI/phiSlices;
            for (unsigned int ieta = 0, neta = etaBoundaries.size()-1; ieta < neta; ++ieta) {
                for (unsigned int iphi = 0; iphi < phiSlices; ++iphi) {
                    float phiCenter = (iphi+0.5)*phiWidth-M_PI;
                    regions_.push_back(Region(
                            etaBoundaries[ieta], etaBoundaries[ieta+1], phiCenter, phiWidth, 
                            phiExtra, etaExtra,
                            9999,9999,999,9,9999,9999)); 
                }
            }
        }
        std::cout << "L1 RegionMapper: made " << regions_.size() << " regions" << std::endl;
    } else {
        // start off with a dummy region
        regions_.push_back(Region(-5.5,5.5, 0,2*M_PI, 0.5, 0.5,
                                 9999,9999,999,9,9999,9999));
    }
}

void RegionMapper::addTrack( const l1tpf::Particle & t ) {
    // now let's be optimistic and make things very simple
    // we propagate in floating point the track to the calo
    // we add the track to the region corresponding to its vertex (eta,phi) coordinates AND its (eta,phi) calo coordinates
    PropagatedTrack prop;
    prop.fillInput(t.pt(), t.eta(), t.phi(), t.charge(), t.dz(), 0);
    prop.fillPropagated(t.pt(), t.sigma(), t.caloSigma(), t.caloEta(), t.caloPhi(), 0);
    for (Region &r : regions_) {
        if (r.contains(t.eta(), t.phi()) || r.contains(t.caloEta(), t.caloPhi())) {
            r.track.push_back(prop);
        }
    } 
}

void RegionMapper::addMuon( const l1tpf::Particle &mu ) {
    // now let's be optimistic and make things very simple
    // we don't propagate anything
    Muon prop;
    prop.fill(mu.pt(), mu.eta(), mu.phi(), mu.charge(), mu.quality());
    for (Region &r : regions_) {
        if (r.contains(mu.eta(), mu.phi())) {
            r.muon.push_back(prop);
        }
    } 
}

void RegionMapper::addCalo( const l1tpf::Particle &p ) { 
    CaloCluster calo;
    calo.fill(p.pt(), p.sigma(), p.eta(), p.phi(), p.pdgId() == l1tpf::Particle::GAMMA, 0);
    for (Region &r : regions_) {
        if (r.contains(p.eta(), p.phi())) {
            r.calo.push_back(calo);
        }
    } 
}
void RegionMapper::addEmCalo( const l1tpf::Particle &p ) { 
    CaloCluster calo;
    calo.fill(p.pt(), p.sigma(), p.eta(), p.phi(), p.pdgId() == l1tpf::Particle::GAMMA, 0);
    for (Region &r : regions_) {
        if (r.contains(p.eta(), p.phi())) {
            r.emcalo.push_back(calo);
        }
    } 
}


std::vector<l1tpf::Particle> RegionMapper::fetch(bool puppi, float ptMin) const {
    std::vector<l1tpf::Particle> ret;
    for (const Region &r : regions_) {
        for (const PFParticle & p : (puppi ? r.puppi : r.pf)) {
            if (regions_.size() > 1) {
                if (!r.fiducial(p.floatVtxEta(), p.floatVtxPhi())) continue;
            }
            if (p.floatPt() > ptMin) {
                ret.emplace_back( p.floatPt(), p.floatVtxEta(), p.floatVtxPhi(), 0.13f, p.hwId, 0.f, p.floatDZ(), p.floatEta(), p.floatPhi(), p.intCharge()  );
            }
        }
    }
    return ret;
}

std::vector<l1tpf::Particle> RegionMapper::fetchCalo(float ptMin) const {
    std::vector<l1tpf::Particle> ret;
    for (const Region &r : regions_) {
        for (const CaloCluster & p : r.calo) {
            if (regions_.size() > 1) {
                if (!r.fiducial(p.floatEta(), p.floatPhi())) continue;
            }
            if (p.floatPt() > ptMin) {
                ret.emplace_back( p.floatPt(), p.floatEta(), p.floatPhi(), 0.13f, p.isEM ? l1tpf::Particle::GAMMA : l1tpf::Particle::NH );
            }
        }
    }
    return ret;
}

std::vector<l1tpf::Particle> RegionMapper::fetchTracks(float ptMin) const {
    std::vector<l1tpf::Particle> ret;
    for (const Region &r : regions_) {
        for (const PropagatedTrack & p : r.track) {
            if (regions_.size() > 1) {
                if (!r.fiducial(p.floatVtxEta(), p.floatVtxPhi())) continue;
            }
            if (p.floatPt() > ptMin) {
                ret.emplace_back( p.floatVtxPt(), p.floatVtxEta(), p.floatVtxPhi(), 0.13f, p.muonLink ? l1tpf::Particle::MU : l1tpf::Particle::CH, 0.f, p.floatDZ(), p.floatEta(), p.floatPhi(), p.intCharge() );
            }
        }
    }
    return ret;
}


PFAlgo::PFAlgo( const edm::ParameterSet & iConfig ) :
    skipMuons_(iConfig.getParameter<bool>("metRate")),
    etaCharged_(iConfig.getParameter<double>("etaCharged")),
    puppiDr_(iConfig.getParameter<double>("puppiDr")),
    puppiPtCutC_(1*iConfig.getParameter<double>("puppiPtCut")),
    puppiPtCutF_(2*iConfig.getParameter<double>("puppiPtCut")),
    vtxCut_(1.5*iConfig.getParameter<double>("vtxRes")),
    drMatch_(0.2), ptMatchLow_(2.0), ptMatchHigh_(2.0), maxInvisiblePt_(20.0),
    useTrackCaloSigma_(false), rescaleUnmatchedTrack_(false),
    debug_(iConfig.getUntrackedParameter<int>("debug",0))
{
    if (iConfig.existsAs<edm::ParameterSet>("linking")) {
        edm::ParameterSet linkcfg = iConfig.getParameter<edm::ParameterSet>("linking");
        drMatch_ = linkcfg.getParameter<double>("trackCaloDR");
        ptMatchLow_ = linkcfg.getParameter<double>("trackCaloNSigmaLow");
        ptMatchHigh_ = linkcfg.getParameter<double>("trackCaloNSigmaHigh");
        useTrackCaloSigma_ = linkcfg.getParameter<bool>("useTrackCaloSigma");
        rescaleUnmatchedTrack_ = linkcfg.getParameter<bool>("rescaleUnmatchedTrack");
        maxInvisiblePt_  = linkcfg.getParameter<double>("maxInvisiblePt");
        if (rescaleUnmatchedTrack_) std::cout << "WARNING: rescaleUnmatchedTrack not yet implemented for integer code" << std::endl;
    }
    intDrMatchBox_ = std::ceil(drMatch_ * CaloCluster::ETAPHI_SCALE * std::sqrt(M_PI/4));
    intPtMatchLowX4_ = std::ceil(ptMatchLow_ * 4);
    intPtMatchHighX4_ = std::ceil(ptMatchHigh_ * 4);
    intMaxInvisiblePt_ = std::round(maxInvisiblePt_ * CaloCluster::PT_SCALE);

    intDrMuonMatchBox_ = std::ceil(0.20 * CaloCluster::ETAPHI_SCALE * std::sqrt(M_PI/4));
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
PFParticle & PFAlgo::addTrackToPF(Region &r, const PropagatedTrack &tk) const {
    PFParticle pf;
    pf.hwPt = tk.hwPt;
    pf.hwEta = tk.hwEta;
    pf.hwPhi = tk.hwPhi;
    pf.hwVtxEta = tk.hwEta; // FIXME: get from the track
    pf.hwVtxPhi = tk.hwPhi; // before propagation
    pf.track = tk;
    pf.cluster.hwPt = 0;
    pf.hwId = (tk.muonLink ? l1tpf::Particle::MU : l1tpf::Particle::CH);
    r.pf.push_back(pf);
    return r.pf.back();
}
PFParticle & PFAlgo::addCaloToPF(Region &r, const CaloCluster &calo) const {
    PFParticle pf;
    pf.hwPt = calo.hwPt;
    pf.hwEta = calo.hwEta;
    pf.hwPhi = calo.hwPhi;
    pf.hwVtxEta = calo.hwEta; 
    pf.hwVtxPhi = calo.hwPhi; 
    pf.track.hwPt = 0;
    pf.cluster = calo;
    pf.hwId = (calo.isEM ? l1tpf::Particle::GAMMA : l1tpf::Particle::NH);
    r.pf.push_back(pf);
    return r.pf.back();
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

void PFAlgo::runPuppi(Region &r, float z0, float npu, float alphaCMed, float alphaCRms, float alphaFMed, float alphaFRms) const {
    makeChargedPV(r, z0);
    computePuppiWeights(r, alphaCMed, alphaCRms, alphaFMed, alphaFRms);
    fillPuppi(r);
}
void PFAlgo::makeChargedPV(Region &r, float z0) const {
    int16_t iZ0 = round(z0 * InputTrack::Z0_SCALE);
    int16_t iDZ = round(vtxCut_ * InputTrack::Z0_SCALE);
    for (PFParticle & p : r.pf) {
        p.chargedPV = (p.hwId <= 1 && std::abs(p.track.hwZ0 - iZ0) < iDZ);
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
        if (std::abs(p.hwEta) < ietacut) {
            if (alphaC > 0) {
                float alpha = std::log(alphaC);
                float x2 = (alpha - alphaCMed) * std::abs(alpha - alphaCMed) / std::pow(alphaCRms,2);
                p.setPuppiW( ROOT::Math::chisquared_cdf(x2,1) );
            } else {
                p.setPuppiW(0);
            }
        } else {
            if (alphaF > 0) {
                float alpha = std::log(alphaF);
                float x2 = (alpha - alphaFMed) * std::abs(alpha - alphaFMed) / std::pow(alphaFRms,2);
                p.setPuppiW( ROOT::Math::chisquared_cdf(x2,1) );
            } else {
                p.setPuppiW(0);
            }
        }
    }
}

void PFAlgo::fillPuppi(Region &r) const {
    int16_t ietacut = std::round(etaCharged_ * CaloCluster::ETAPHI_SCALE);
    int16_t iptcutC = std::round(puppiPtCutC_ * CaloCluster::PT_SCALE);
    int16_t iptcutF = std::round(puppiPtCutF_ * CaloCluster::PT_SCALE);
    constexpr uint16_t PUPPIW_0p01 = std::round(0.01 * PFParticle::PUPPI_SCALE);
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
                if (hwPt > (std::abs(p.hwEta) < ietacut ? iptcutC : iptcutF)) {
                    r.puppi.push_back(p);
                    r.puppi.back().hwPt = hwPt;
                }
            }
        }
    }
}

