#ifndef FASTPUPPI_NTUPLERPRODUCER_FIRMWARE2DISCRETEPF_H
#define FASTPUPPI_NTUPLERPRODUCER_FIRMWARE2DISCRETEPF_H

/// NOTE: this include is not standalone, since the path to DiscretePFInputs is different in CMSSW & Vivado_HLS

#include "firmware/data.h"
#include <vector>
#include <cassert>

namespace fw2dpf {

    // convert inputs from discrete to firmware
    void convert(const PFChargedObj & src, const l1tpf_int::PropagatedTrack & track, std::vector<l1tpf_int::PFParticle> &out) {
        l1tpf_int::PFParticle pf;
        pf.hwPt = src.hwPt;
        pf.hwEta = src.hwEta;
        pf.hwPhi = src.hwPhi;
        pf.hwVtxEta = src.hwEta; // FIXME: get from the track
        pf.hwVtxPhi = src.hwPhi; // before propagation
        pf.track = track; // FIXME: ok only as long as there is a 1-1 mapping
        pf.cluster.hwPt = 0;
        switch(src.hwId) {
            case PID_Electron: pf.hwId =  l1tpf::Particle::EL; break;
            case PID_Muon: pf.hwId =  l1tpf::Particle::MU; break;
            default: pf.hwId =  l1tpf::Particle::CH; break;
        };
        pf.hwStatus = 0;
        out.push_back(pf);
    }
    void convert(const PFNeutralObj & src, std::vector<l1tpf_int::PFParticle> &out) {
        l1tpf_int::PFParticle pf;
        pf.hwPt = src.hwPt;
        pf.hwEta = src.hwEta;
        pf.hwPhi = src.hwPhi;
        pf.hwVtxEta = src.hwEta;
        pf.hwVtxPhi = src.hwPhi;
        pf.track.hwPt = 0;
        pf.cluster.hwPt = src.hwPt;
        switch(src.hwId) {
            case PID_Photon: pf.hwId = l1tpf::Particle::GAMMA; break;
            default: pf.hwId = l1tpf::Particle::NH; break;
        }
        pf.hwStatus = 0;
        out.push_back(pf);
    }

    template<unsigned int NMAX, typename In>
    void convert(const In in[NMAX], std::vector<l1tpf_int::PFParticle> &out) {
        for (unsigned int i = 0; i < NMAX; ++i) {
            if (in[i].hwPt > 0) convert(in[i], out);
        }
    } 
    template<unsigned int NMAX>
    void convert(const PFChargedObj in[NMAX], std::vector<l1tpf_int::PropagatedTrack> srctracks, std::vector<l1tpf_int::PFParticle> &out) {
        for (unsigned int i = 0; i < NMAX; ++i) {
            if (in[i].hwPt > 0) {
                assert(i < srctracks.size());
                convert(in[i], srctracks[i], out);
            }
        }
    }

} // namespace

#endif
