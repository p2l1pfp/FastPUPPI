#include "../interface/BitwisePF.h"
#include "firmware/simple_fullpfalgo.h"
#include "DiscretePF2Firmware.h"
#include "Firmware2DiscretePF.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"

int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
    etaphi_t deta = (eta1-eta2);
    etaphi_t dphi = (phi1-phi2);
    return deta*deta + dphi*dphi;
}

using namespace l1tpf_int;

BitwisePF::BitwisePF( const edm::ParameterSet & iConfig ) :
    PFAlgo(iConfig)
{
    debug_ = iConfig.getUntrackedParameter<int>("bitwiseDebug", debug_);
    pfalgo3_full_ref_set_debug(debug_);
}

void BitwisePF::runPF(Region &r) const {
    initRegion(r);

    if (debug_) {
        printf("BW\nBW  region Eta [ %+5.2f , %+5.2f ],  Phi [ %+5.2f , %+5.2f ] \n", r.etaMin, r.etaMax, r.phiCenter-r.phiHalfWidth, r.phiCenter+r.phiHalfWidth );
        printf("BW  \t N(track) %3lu   N(em) %3lu   N(calo) %3lu   N(mu) %3lu\n", r.track.size(), r.emcalo.size(), r.calo.size(), r.muon.size());
        for (int itk = 0, ntk = r.track.size(); itk < ntk; ++itk) {
            const auto & tk = r.track[itk]; 
            printf("BW  \t track %3d: pt %7.2f +- %5.2f  vtx eta %+5.2f  vtx phi %+5.2f  calo eta %+5.2f  calo phi %+5.2f calo ptErr %7.2f stubs %2d chi2 %7.1f\n", 
                                itk, tk.floatPt(), tk.floatPtErr(), tk.floatVtxEta(), tk.floatVtxPhi(), tk.floatEta(), tk.floatPhi(), tk.floatCaloPtErr(), int(tk.hwStubs), tk.hwChi2*0.1f);
        }
        for (int iem = 0, nem = r.emcalo.size(); iem < nem; ++iem) {
            const auto & em = r.emcalo[iem];
            printf("BW  \t EM    %3d: pt %7.2f +- %5.2f  vtx eta %+5.2f  vtx phi %+5.2f  calo eta %+5.2f  calo phi %+5.2f calo ptErr %7.2f\n", 
                                iem, em.floatPt(), em.floatPtErr(), em.floatEta(), em.floatPhi(), em.floatEta(), em.floatPhi(), em.floatPtErr());
        } 
        for (int ic = 0, nc = r.calo.size(); ic < nc; ++ic) {
            auto & calo = r.calo[ic]; 
            printf("BW  \t calo  %3d: pt %7.2f +- %5.2f  vtx eta %+5.2f  vtx phi %+5.2f  calo eta %+5.2f  calo phi %+5.2f calo ptErr %7.2f em pt %7.2f \n", 
                                ic, calo.floatPt(), calo.floatPtErr(), calo.floatEta(), calo.floatPhi(), calo.floatEta(), calo.floatPhi(), calo.floatPtErr(), calo.floatEmPt());
        }
        for (int im = 0, nm = r.muon.size(); im < nm; ++im) {
            auto & mu = r.muon[im]; 
            printf("BW  \t muon  %3d: pt %7.2f           vtx eta %+5.2f  vtx phi %+5.2f  calo eta %+5.2f  calo phi %+5.2f \n", 
                                im, mu.floatPt(), mu.floatEta(), mu.floatPhi(), mu.floatEta(), mu.floatPhi());
        }
    }

    EmCaloObj emcalo[NEMCALO];
    HadCaloObj hadcalo[NCALO];
    TkObj track[NTRACK];
    MuObj mu[NMU];
    PFChargedObj outch[NTRACK];
    PFNeutralObj outpho[NPHOTON];
    PFNeutralObj outne[NSELCALO];
    PFChargedObj outmu[NMU];

    dpf2fw::convert<NTRACK>(r.track, track);
    dpf2fw::convert<NCALO>(r.calo, hadcalo);
    dpf2fw::convert<NEMCALO>(r.emcalo, emcalo);
    dpf2fw::convert<NMU>(r.muon, mu);

    pfalgo3_full_ref(emcalo, hadcalo, track, mu, outch, outpho, outne, outmu);

    r.pf.clear();
    fw2dpf::convert<NTRACK>(outch, r.track, r.pf); // FIXME works only with a 1-1 mapping
    fw2dpf::convert<NSELCALO>(outne, r.pf);
    fw2dpf::convert<NPHOTON>(outpho, r.pf);
    // muons and electrons are already included in outch
}

