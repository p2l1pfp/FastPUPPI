#ifndef FASTPUPPI_NTUPLERPRODUCER_DISCRETEPF_H
#define FASTPUPPI_NTUPLERPRODUCER_DISCRETEPF_H
/** 
 * Classes for integer-level L1 PF
 * 
 * */

// fwd declarations
namespace edm { class ParameterSet; }

// real includes
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include "DiscretePFInputs.h"
#include "L1TPFParticle.h"
#include "DataFormats/Math/interface/deltaPhi.h"

namespace l1tpf_int { 


  struct Region : public InputRegion {
    std::vector<PFParticle>       pf;
    std::vector<PFParticle>       puppi;
    std::vector<PFParticle>       pfdiscarded; // for debugging
    unsigned int caloOverflow, emcaloOverflow, trackOverflow, muonOverflow, pfOverflow, puppiOverflow;

    const bool relativeCoordinates; // whether the eta,phi in each region are global or relative to the region center
    const unsigned int ncaloMax, nemcaloMax, ntrackMax, nmuonMax, npfMax, npuppiMax;
    Region(float etamin, float etamax, float phicenter, float phiwidth, float etaextra, float phiextra, bool useRelativeCoordinates,
           unsigned int ncalomax, unsigned int nemcalomax, unsigned int ntrackmax, unsigned int nmuonmax, unsigned int npfmax, unsigned int npuppimax) :
        InputRegion(0.5*(etamin+etamax), etamin, etamax, phicenter, 0.5*phiwidth, etaextra, phiextra), relativeCoordinates(useRelativeCoordinates),
        ncaloMax(ncalomax), nemcaloMax(nemcalomax), ntrackMax(ntrackmax), nmuonMax(nmuonmax), npfMax(npfmax), npuppiMax(npuppimax) {}

    enum InputType { calo_type=0, emcalo_type=1, track_type=2, l1mu_type=3, n_input_types=4 };
    static const char * inputTypeName(int inputType) ;

    enum OutputType { any_type=0, charged_type=1, neutral_type=2, electron_type=3, pfmuon_type=4, charged_hadron_type=5, neutral_hadron_type=6, photon_type=7, n_output_types=8 };
    static const char * outputTypeName(int outputType) ;

    unsigned int nInput(InputType type) const ;
    unsigned int nOutput(OutputType type, bool puppi) const ;

    // global coordinates
    bool contains(float eta, float phi) const { 
        float dphi = deltaPhi(phiCenter,phi);
        return (etaMin-etaExtra < eta && eta <= etaMax+etaExtra && 
                -phiHalfWidth-phiExtra < dphi && dphi <= phiHalfWidth+phiExtra); 
    }
    // global coordinates
    bool fiducial(float eta, float phi) const { 
        float dphi = deltaPhi(phiCenter,phi);
        return (etaMin < eta && eta <= etaMax && 
                -phiHalfWidth < dphi && dphi <= phiHalfWidth); 
    }
    // possibly local coordinates
    bool fiducialLocal(float localEta, float localPhi) const { 
        if (relativeCoordinates) {
            float dphi = deltaPhi(0.f,localPhi);
            return (etaMin < localEta+etaCenter && localEta+etaCenter <= etaMax && 
                -phiHalfWidth < dphi && dphi <= phiHalfWidth); 
        }
        float dphi = deltaPhi(phiCenter,localPhi);
        return (etaMin < localEta && localEta <= etaMax && 
                -phiHalfWidth < dphi && dphi <= phiHalfWidth); 
    }
    float regionAbsEta() const { return std::abs(etaCenter); }
    float globalAbsEta(float localEta) const { return std::abs(relativeCoordinates ? localEta + etaCenter : localEta); }
    float globalEta(float localEta) const { return relativeCoordinates ? localEta + etaCenter : localEta; }
    float globalPhi(float localPhi) const { return relativeCoordinates ? localPhi + phiCenter : localPhi; }
    float localEta(float globalEta) const { return relativeCoordinates ? globalEta - etaCenter : globalEta; }
    float localPhi(float globalPhi) const { return relativeCoordinates ? deltaPhi(globalPhi,phiCenter) : globalPhi; }

    void zero() {
        calo.clear(); emcalo.clear(); track.clear(); muon.clear(); pf.clear(); puppi.clear();
        caloOverflow = 0; emcaloOverflow = 0; trackOverflow = 0; muonOverflow = 0; pfOverflow = 0; puppiOverflow = 0;
    }
    void inputSort() {
        std::sort(calo.begin(),  calo.end());
        std::sort(emcalo.begin(),  emcalo.end());
        std::sort(track.begin(), track.end());
        std::sort(muon.begin(),  muon.end());
        if (ncaloMax > 0 && calo.size() > ncaloMax) { caloOverflow = calo.size() - ncaloMax; calo.resize(ncaloMax); }
        if (nemcaloMax > 0 && emcalo.size() > nemcaloMax) { emcaloOverflow = emcalo.size() - nemcaloMax; emcalo.resize(nemcaloMax); }
        if (ntrackMax > 0 && track.size() > ntrackMax) { trackOverflow = track.size() - ntrackMax; track.resize(ntrackMax); }
    }
    void outputSort() {
        std::sort(puppi.begin(), puppi.end());
    }

  };

  class RegionMapper {
    // This does the input and filling of regions. 
    public:
        RegionMapper( const edm::ParameterSet& ) ;
        void addTrack( const l1tpf::Particle & t ) ;
        void addMuon( const l1tpf::Particle & t );
        void addCalo( const l1tpf::Particle & t ); 
        void addEmCalo( const l1tpf::Particle & t ); 

        void clear() { for (Region & r : regions_) r.zero(); }
        std::vector<Region> & regions() { return regions_; }

        std::vector<l1tpf::Particle> fetch(bool puppi=true, float ptMin=0.01, bool discarded = false) const ;
        std::vector<l1tpf::Particle> fetchCalo(float ptMin=0.01, bool emcalo=false) const ;
        std::vector<l1tpf::Particle> fetchTracks(float ptMin=0.01, bool fromPV=false) const ;

        std::pair<unsigned,unsigned> totAndMaxInput(/*Region::InputType*/int type) const ;
        std::pair<unsigned,unsigned> totAndMaxOutput(/*Region::OutputType*/int type, bool puppi) const ;
    protected:
        std::vector<Region> regions_;
        bool useRelativeRegionalCoordinates_; // whether the eta,phi in each region are global or relative to the region center
        enum TrackAssoMode { atVertex, atCalo, any=999 } trackRegionMode_;
  };

  class PFAlgo {
    public:
        PFAlgo( const edm::ParameterSet& ) ;
        virtual void runPF(Region &r) const ;
        virtual void runChargedPV(Region &r, float z0) const ;
        virtual void runPuppi(Region &r, float npu, float alphaCMed, float alphaCRms, float alphaFMed, float alphaFRms) const ;
        /// global operations
        enum VertexAlgo { OldVtxAlgo, TPVtxAlgo };
        virtual void doVertexing(std::vector<Region> &rs, VertexAlgo algo, float &vz) const ; // region is not const since it sets the fromPV bit of the tracks
        virtual void computePuppiMedRMS(const std::vector<Region> &rs, float &alphaCMed, float &alphaCRms, float &alphaFMed, float &alphaFRms) const ;
    protected:
        bool skipMuons_;
        float etaCharged_, puppiDr_; 
        std::vector<float> puppiEtaCuts_, puppiPtCuts_, puppiPtCutsPhotons_;
        std::vector<int16_t> intPuppiEtaCuts_, intPuppiPtCuts_, intPuppiPtCutsPhotons_;
        float vtxRes_;
        bool vtxAdaptiveCut_; 
        float drMatch_, ptMatchLow_, ptMatchHigh_, maxInvisiblePt_;
        int16_t intDrMuonMatchBox_, intDrMatchBox_, intPtMatchLowX4_, intPtMatchHighX4_, intMaxInvisiblePt_;
        bool useTrackCaloSigma_, rescaleUnmatchedTrack_;
        int debug_;
        void initRegion(Region &r) const ;
        void muonTrackLink(Region &r) const ;
        void caloTrackLink(Region &r) const ;
        void mergeTkCalo(Region &r, const PropagatedTrack &tk, CaloCluster & calo) const ;
        void computePuppiWeights(Region &r, float alphaCMed, float alphaCRms, float alphaFMed, float alphaFRms) const ;
        void fillPuppi(Region &r) const ;
        PFParticle & addTrackToPF(Region &r, const PropagatedTrack &tk) const { return addTrackToPF(r.pf, tk); }
        PFParticle & addCaloToPF(Region &r, const CaloCluster &calo) const { return addCaloToPF(r.pf, calo); }
        PFParticle & discardTrack(Region &r, const PropagatedTrack &tk, int status) const { 
            PFParticle & ret = addTrackToPF(r.pfdiscarded, tk); 
            ret.hwStatus = status;
            return ret;
        }
        PFParticle & discardCalo(Region &r, const CaloCluster &calo, int status) const { 
            PFParticle & ret = addCaloToPF(r.pfdiscarded, calo); 
            ret.hwStatus = status;
            return ret;
        }
        PFParticle & addTrackToPF(std::vector<PFParticle> &pfs, const PropagatedTrack &tk) const ;
        PFParticle & addCaloToPF(std::vector<PFParticle> &pfs, const CaloCluster &calo) const ;
  };

} // end namespace

#endif
