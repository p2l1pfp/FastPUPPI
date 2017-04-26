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
#include <cmath>
#include <vector>
#include <algorithm>
#include "L1TPFParticle.h"

namespace l1tpf_int { 

  struct CaloCluster {
      int16_t hwPt;   
      int16_t  hwPtErr;   
      int16_t  hwEta;   
      int16_t  hwPhi;   
      uint16_t hwFlags;
      bool     isEM, used;
      static constexpr float PT_SCALE = 4.0;     // quantize in units of 0.25 GeV (can be changed)
      static constexpr float ETAPHI_FACTOR = 4;  // size of an ecal crystal in phi in integer units (our choice)
      static constexpr float ETAPHI_SCALE = ETAPHI_FACTOR*(180./M_PI);  // M_PI/180 is the size of an ECal crystal; we make a grid that is 4 times that size
      static constexpr int16_t PHI_WRAP = 360*ETAPHI_FACTOR;            // what is 3.14 in integer
      // sorting
      bool operator<(const CaloCluster &other) const { return hwPt > other.hwPt; }
      // filling from floating point
      void fill(float pt, float ptErr, float eta, float phi, bool em, unsigned int flags) {
          hwPt  = round(pt  * CaloCluster::PT_SCALE);
          hwPtErr = round(ptErr  * CaloCluster::PT_SCALE);
          hwEta = round(eta * CaloCluster::ETAPHI_SCALE);
          hwPhi = int16_t(round(phi * CaloCluster::ETAPHI_SCALE)) % CaloCluster::PHI_WRAP;
          isEM  = em;
          used = false;
          hwFlags = flags;
      }
      float floatPt() const { return float(hwPt) / CaloCluster::PT_SCALE; }
      float floatPtErr() const { return float(hwPtErr) / CaloCluster::PT_SCALE; }
      float floatEta() const { return float(hwEta) / CaloCluster::ETAPHI_SCALE; }
      float floatPhi() const { return float(hwPhi) / CaloCluster::ETAPHI_SCALE; }
  };

  // https://twiki.cern.ch/twiki/bin/view/CMS/L1TriggerPhase2InterfaceSpecifications
  struct InputTrack {
      uint16_t hwInvpt;
      int32_t  hwVtxEta;
      int32_t  hwVtxPhi;
      bool     hwCharge;
      int16_t  hwZ0;
      uint16_t hwFlags;
      static constexpr float INVPT_SCALE   = 2E4;    // 1%/pt @ 100 GeV is 2 bits 
      static constexpr float VTX_PHI_SCALE = 1/2.5E-6; // 5 micro rad is 2 bits
      static constexpr float VTX_ETA_SCALE = 1/1E-5;   // no idea, but assume it's somewhat worse than phi
      static constexpr float Z0_SCALE      = 20;     // 1mm is 2 bits
      // filling from floating point
      void fillInput(float pt, float eta, float phi, int charge, float dz, unsigned int flags) {
          hwInvpt  = round(1/pt  * InputTrack::INVPT_SCALE);
          hwVtxEta = round(eta * InputTrack::VTX_ETA_SCALE);
          hwVtxPhi = round(phi * InputTrack::VTX_PHI_SCALE);
          hwCharge = (charge > 0);
          hwZ0     = round(dz  * InputTrack::Z0_SCALE);
          hwFlags = flags;
      }
      float floatVtxPt() const { return 1/(float(hwInvpt) / InputTrack::INVPT_SCALE); }
      float floatVtxEta() const { return float(hwVtxEta) / InputTrack::VTX_ETA_SCALE; }
      float floatVtxPhi() const { return float(hwVtxPhi) / InputTrack::VTX_PHI_SCALE; }
      float floatDZ()     const { return float(hwZ0) / InputTrack::Z0_SCALE; }
      int intCharge()     const { return hwCharge ? +1 : -1; }
  };

  struct PropagatedTrack : public InputTrack {
      int16_t  hwPt;
      int16_t  hwPtErr;
      int16_t  hwEta; // at calo
      int16_t  hwPhi; // at calo
      bool     muonLink;
      // sorting
      bool operator<(const PropagatedTrack &other) const { return hwPt > other.hwPt; }
      void fillPropagated(float pt, float ptErr, float eta, float phi, unsigned int flags) {
          hwPt  = round(pt  * CaloCluster::PT_SCALE);
          hwPtErr = round(ptErr  * CaloCluster::PT_SCALE);
          hwEta = round(eta * CaloCluster::ETAPHI_SCALE);
          hwPhi = int16_t(round(phi * CaloCluster::ETAPHI_SCALE)) % CaloCluster::PHI_WRAP;
          muonLink = false;
      }
      float floatPt() const { return float(hwPt) / CaloCluster::PT_SCALE; }
      float floatPtErr() const { return float(hwPtErr) / CaloCluster::PT_SCALE; }
      float floatEta() const { return float(hwEta) / CaloCluster::ETAPHI_SCALE; }
      float floatPhi() const { return float(hwPhi) / CaloCluster::ETAPHI_SCALE; }
  };

  struct Muon {
      int16_t hwPt;   
      int16_t  hwEta;   // at calo
      int16_t  hwPhi;   // at calo
      uint16_t hwFlags;
      bool     hwCharge;
      // sorting
      bool operator<(const Muon &other) const { return hwPt > other.hwPt; }
      void fill(float pt, float eta, float phi, int charge, unsigned int flags) {
          // we assume we use the same discrete ieta, iphi grid for all particles 
          hwPt  = round(pt  * CaloCluster::PT_SCALE);
          hwEta = round(eta * CaloCluster::ETAPHI_SCALE);
          hwPhi = int16_t(round(phi * CaloCluster::ETAPHI_SCALE)) % CaloCluster::PHI_WRAP;
          hwCharge = (charge > 0);
          hwFlags = flags;
      }
      float floatPt() const { return float(hwPt) / CaloCluster::PT_SCALE; }
      float floatEta() const { return float(hwEta) / CaloCluster::ETAPHI_SCALE; }
      float floatPhi() const { return float(hwPhi) / CaloCluster::ETAPHI_SCALE; }
      int intCharge()     const { return hwCharge ? +1 : -1; }
  };

  struct PFParticle {
      enum PID { CH=0, EL=1, NH=2, GAMMA=3, MU=4 };
      int16_t         hwPt;   
      int16_t         hwEta;  // at calo face 
      int16_t         hwPhi;   
      uint8_t         hwId; // CH=0, EL=1, NH=2, GAMMA=3, MU=4 
      int16_t         hwVtxEta;  // propagate back to Vtx for charged particles (if useful?)
      int16_t         hwVtxPhi;   
      uint16_t        hwFlags;
      CaloCluster     cluster;
      PropagatedTrack track;
      bool            chargedPV;
      uint16_t        hwPuppiWeight;
      static constexpr float PUPPI_SCALE = 100;
      // sorting
      bool operator<(const PFParticle &other) const { return hwPt > other.hwPt; }
      float floatPt() const { return float(hwPt) / CaloCluster::PT_SCALE; }
      float floatEta() const { return float(hwEta) / CaloCluster::ETAPHI_SCALE; }
      float floatPhi() const { return float(hwPhi) / CaloCluster::ETAPHI_SCALE; }
      float floatDZ() const { return float(track.hwZ0) / InputTrack::Z0_SCALE; }
      int intCharge()     const { return (track.hwPt > 0 ? track.intCharge() : 0); }
      void setPuppiW(float w) {
            hwPuppiWeight = std::round(w * PUPPI_SCALE);
      }
  };

  struct Region {
    std::vector<CaloCluster>      calo;
    std::vector<PropagatedTrack>  track;
    std::vector<Muon>             muon;
    std::vector<PFParticle>       pf;
    std::vector<PFParticle>       puppi;
    unsigned int caloOverflow, trackOverflow, muonOverflow, pfOverflow, puppiOverflow;

    const float etaMin, etaMax, phiMin, phiMax;
    const unsigned int ncaloMax, ntrackMax, nmuonMax, npfMax, npuppiMax;
    Region(float etamin, float etamax, float phimin, float phimax,
           unsigned int ncalomax, unsigned int ntrackmax, unsigned int nmuonmax, unsigned int npfmax, unsigned int npuppimax) :
        etaMin(etamin), etaMax(etamax), phiMin(phimin), phiMax(phimax),
        ncaloMax(ncalomax), ntrackMax(ntrackmax), nmuonMax(nmuonmax), npfMax(npfmax), npuppiMax(npuppimax) {}

    bool contains(float eta, float phi) const { return (etaMin <= eta && eta <= etaMax && phiMin <= phi && phi <= phiMax); }

    void zero() {
        calo.clear(); track.clear(); muon.clear(); pf.clear(); puppi.clear();
        caloOverflow = 0; trackOverflow = 0; muonOverflow = 0; pfOverflow = 0; puppiOverflow = 0;
    }
    void inputSort() {
        std::sort(calo.begin(),  calo.end());
        std::sort(track.begin(), track.end());
        std::sort(muon.begin(),  muon.end());
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

        void clear() { for (Region & r : regions_) r.zero(); }
        std::vector<Region> & regions() { return regions_; }

        std::vector<l1tpf::Particle> fetch(bool puppi=true, float ptMin=0.01) const ;
        std::vector<l1tpf::Particle> fetchCalo(float ptMin=0.01) const ;
        std::vector<l1tpf::Particle> fetchTracks(float ptMin=0.01) const ;
    protected:
        std::vector<Region> regions_;
  };

  class PFAlgo {
    public:
        PFAlgo( const edm::ParameterSet& ) ;
        void runPF(Region &r) const ;
        void runPuppi(Region &r, float z0, float npu, float alphaCMed, float alphaCRms, float alphaFMed, float alphaFRms) const ;
    private:
        bool skipMuons_;
        float etaCharged_, puppiPtCutC_, puppiPtCutF_, vtxCut_;
        int debug_;
        PFParticle & addTrackToPF(Region &r, const PropagatedTrack &tk) const ;
        PFParticle & addCaloToPF(Region &r, const CaloCluster &calo) const ;
        void mergeTkCalo(Region &r, const PropagatedTrack &tk, CaloCluster & calo) const ;
        void makeChargedPV(Region &r, float z0) const ;
        void computePuppiWeights(Region &r, float alphaCMed, float alphaCRms, float alphaFMed, float alphaFRms) const ;
        void fillPuppi(Region &r) const ;
  };

} // end namespace

#endif
