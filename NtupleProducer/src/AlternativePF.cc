#include "../interface/AlternativePF.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Math/ProbFunc.h"

namespace {
    template <typename T1, typename T2>
    float floatDR(const T1 &t1, const T2 &t2) { return deltaR(t1.floatEta(), t1.floatPhi(), t2.floatEta(), t2.floatPhi()); }
}

using namespace l1tpf_int;

PFAlgo3::PFAlgo3( const edm::ParameterSet & iConfig ) :
    PFAlgo(iConfig),
    drMatchEm_(0.2), drMatchEmHad_(0.2)
{
    debug_ = iConfig.getUntrackedParameter<int>("altDebug", debug_);
    edm::ParameterSet linkcfg = iConfig.getParameter<edm::ParameterSet>("linking");
    drMatchEm_ = linkcfg.getParameter<double>("trackEmDR");
    drMatchEmHad_ = linkcfg.getParameter<double>("emCaloDR");
    caloReLinkStep_ = linkcfg.getParameter<bool>("caloReLink");
    caloReLinkDr_ = linkcfg.getParameter<double>("caloReLinkDR");
    caloReLinkThreshold_ = linkcfg.getParameter<double>("caloReLinkThreshold");
    sumTkCaloErr2_ = linkcfg.getParameter<bool>("sumTkCaloErr2");
    ecalPriority_ = linkcfg.getParameter<bool>("ecalPriority");
}

void PFAlgo3::runPF(Region &r) const {
    initRegion(r);

    /// ------------- first step (can all go in parallel) ----------------

    muonTrackLink(r);

    if (debug_) {
        printf("ALT \t region Eta [ %+5.2f , %+5.2f ],  Phi [ %+5.2f , %+5.2f ] \n", r.etaMin, r.etaMax, r.phiCenter-r.phiHalfWidth, r.phiCenter+r.phiHalfWidth );
        printf("ALT \t N(track) %3lu   N(em) %3lu   N(calo) %3lu\n", r.track.size(), r.emcalo.size(), r.calo.size());
        for (int itk = 0, ntk = r.track.size(); itk < ntk; ++itk) {
            const auto & tk = r.track[itk]; 
            printf("ALT \t track %3d: pt %7.2f +- %5.2f  vtx eta %+5.2f  vtx phi %+5.2f  calo eta %+5.2f  calo phi %+5.2f calo ptErr %7.2f\n", 
                                itk, tk.floatPt(), tk.floatPtErr(), tk.floatVtxEta(), tk.floatVtxPhi(), tk.floatEta(), tk.floatPhi(), tk.floatCaloPtErr());
        }
        for (int iem = 0, nem = r.emcalo.size(); iem < nem; ++iem) {
            const auto & em = r.emcalo[iem];
            printf("ALT \t EM    %3d: pt %7.2f +- %5.2f  vtx eta %+5.2f  vtx phi %+5.2f  calo eta %+5.2f  calo phi %+5.2f calo ptErr %7.2f\n", 
                                iem, em.floatPt(), em.floatPtErr(), em.floatEta(), em.floatPhi(), em.floatEta(), em.floatPhi(), em.floatPtErr());
        } 
        for (int ic = 0, nc = r.calo.size(); ic < nc; ++ic) {
            auto & calo = r.calo[ic]; 
            printf("ALT \t calo  %3d: pt %7.2f +- %5.2f  vtx eta %+5.2f  vtx phi %+5.2f  calo eta %+5.2f  calo phi %+5.2f calo ptErr %7.2f em pt %7.2f \n", 
                                ic, calo.floatPt(), calo.floatPtErr(), calo.floatEta(), calo.floatPhi(), calo.floatEta(), calo.floatPhi(), calo.floatPtErr(), calo.floatEmPt());
        }
    }
    // match all tracks to the closest EM cluster
    std::vector<int> tk2em(r.track.size(), -1);
    for (int itk = 0, ntk = r.track.size(); itk < ntk; ++itk) {
        const auto & tk = r.track[itk]; 
        //if (tk.muonLink) continue; // not necessary I think
        float drbest = drMatchEm_;
        for (int iem = 0, nem = r.emcalo.size(); iem < nem; ++iem) {
            const auto & em = r.emcalo[iem];
            float dr = floatDR(tk, em);
            if (dr < drbest) { tk2em[itk] = iem; drbest = dr; }
        }
        if (debug_ && tk2em[itk] != -1) printf("ALT \t track %3d (pt %7.2f) matches to EM   %3d (pt %7.2f) with dr %.3f\n", itk, tk.floatPt(), tk2em[itk], tk2em[itk] == -1 ? 0.0 : r.emcalo[tk2em[itk]].floatPt(), drbest );
    }

    // match all em to the closest had (can happen in parallel to the above)
    std::vector<int> em2calo(r.emcalo.size(), -1);
    for (int iem = 0, nem = r.emcalo.size(); iem < nem; ++iem) {
        const auto & em = r.emcalo[iem];
        float drbest = drMatchEmHad_;
        for (int ic = 0, nc = r.calo.size(); ic < nc; ++ic) {
            const auto & calo = r.calo[ic]; 
            float dr = floatDR(calo, em);
            if (dr < drbest) { em2calo[iem] = ic; drbest = dr; }
        }
        if (debug_ && em2calo[iem] != -1) printf("ALT \t EM    %3d (pt %7.2f) matches to calo %3d (pt %7.2f) with dr %.3f\n", iem, em.floatPt(), em2calo[iem], em2calo[iem] == -1 ? 0.0 : r.calo[em2calo[iem]].floatPt(), drbest );
    }

    /// ------------- next step (needs the previous) ----------------

    // for each EM cluster, count and add up the pt of all the corresponding tracks (skipping muons)
    std::vector<int> em2ntk(r.emcalo.size(), 0);
    std::vector<float> em2sumtkpt(r.emcalo.size(), 0);
    std::vector<float> em2sumtkpterr(r.emcalo.size(), 0);
    for (int iem = 0, nem = r.emcalo.size(); iem < nem; ++iem) {
        const auto & em = r.emcalo[iem];
        if (std::abs(em.floatEta()) > 2.5) continue; 
        for (int itk = 0, ntk = r.track.size(); itk < ntk; ++itk) {
            if (tk2em[itk] == iem) {
                const auto & tk = r.track[itk]; 
                if (tk.muonLink) continue; 
                em2ntk[iem]++;
                em2sumtkpt[iem] += tk.floatPt();
                em2sumtkpterr[iem] += tk.floatPtErr();
            }
        }
    }

    /// ------------- next step (needs the previous) ----------------

    // process ecal clusters after linking
    for (int iem = 0, nem = r.emcalo.size(); iem < nem; ++iem) {
        auto & em = r.emcalo[iem];
        em.isEM = false; em.used = false; em.hwFlags = 0; 
        if (std::abs(em.floatEta()) > 2.5) continue; 
        if (debug_) printf("ALT \t EM    %3d (pt %7.2f) has %2d tracks (sumpt %7.2f, sumpterr %7.2f), ptdif %7.2f +- %7.2f\n", iem, em.floatPt(), em2ntk[iem], em2sumtkpt[iem], em2sumtkpterr[iem], em.floatPt() - em2sumtkpt[iem], std::max<float>(em2sumtkpterr[iem],em.floatPtErr()));
        if (em2ntk[iem] == 0) { // Photon
            em.isEM = true;
            addCaloToPF(r, em);
            em.used = true;
            if (debug_) printf("ALT \t EM    %3d (pt %7.2f)    ---> promoted to photon\n", iem, em.floatPt());
            continue;
        } 
        float ptdiff = em.floatPt() - em2sumtkpt[iem];
        float pterr  = std::max<float>(em2sumtkpterr[iem],em.floatPtErr());
        if (ptdiff > -ptMatchLow_*pterr) {
            em.isEM = true;
            em.used = true;
            // convert leftover to a photon if significant
            if (ptdiff > +ptMatchHigh_*pterr) {
                auto & p = addCaloToPF(r, em);
                p.setFloatPt(ptdiff);
                if (debug_) printf("ALT \t EM    %3d (pt %7.2f)    ---> promoted to electron(s) + photon (pt %7.2f)\n", iem, em.floatPt(), ptdiff);
            } else {
                em.hwFlags = 1; // may use calo momentum
                if (debug_) printf("ALT \t EM    %3d (pt %7.2f)    ---> promoted to electron(s)\n", iem, em.floatPt());
            }
        } else {
            em.isEM = false;
            em.used = false;
            em.hwFlags = 0; 
            discardCalo(r, em, 2);
        } 
    }

    /// ------------- next step (needs the previous) ----------------

    // promote all flagged tracks to electrons
    for (int itk = 0, ntk = r.track.size(); itk < ntk; ++itk) {
        auto & tk = r.track[itk]; 
        if (tk2em[itk] == -1 || tk.muonLink) continue;
        const auto & em = r.emcalo[tk2em[itk]];
        if (em.isEM) {
            auto & p = addTrackToPF(r, tk);
            // FIXME to check if this is useful
            if (em2ntk[tk2em[itk]] == 1 && em.floatPt() > 20 && em.hwFlags == 1) { 
                p.setFloatPt(em.floatPt()); // calo has better resolution above 20 GeV
            }
            if (debug_) printf("ALT \t track %3d (pt %7.2f) matched to EM   %3d (pt %7.2f) promoted to electron with pt %7.2f\n", itk, tk.floatPt(), tk2em[itk], em.floatPt(), p.floatPt() );
            p.hwId = l1tpf::Particle::EL;
            tk.used = true;
        }
    }

    // subtract EM component from Calo clusters for all photons and electrons (within tracker coverage)
    // kill clusters that end up below their own uncertainty, or that loose 90% of the energy
    for (int ic = 0, nc = r.calo.size(); ic < nc; ++ic) {
        auto & calo = r.calo[ic]; 
        float pt0 = calo.floatPt(), ept0 = calo.floatEmPt(), pt = pt0, ept = ept0;
        for (int iem = 0, nem = r.emcalo.size(); iem < nem; ++iem) {
            if (em2calo[iem] == ic) {
                const auto & em = r.emcalo[iem];
                if (em.isEM) {
                    if (debug_) printf("ALT \t EM    %3d (pt %7.2f) subtracted from calo %3d (pt %7.2f)\n", iem, em.floatPt(), ic, calo.floatPt());
                    pt  -= em.floatPt();
                    ept -= em.floatPt();
                }
            }
        }
        if (pt < pt0) {
            if (debug_) printf("ALT \t calo  %3d (pt %7.2f +- %7.2f) has a subtracted pt of %7.2f, empt %7.2f -> %7.2f\n", ic, calo.floatPt(), calo.floatPtErr(), pt, ept0, ept);
            calo.setFloatPt(pt);
            calo.setFloatEmPt(ept);
            if (pt < calo.floatPtErr() || pt < 0.1*pt0 || (calo.isEM && ept < 0.1*ept0)) {
                if (debug_) printf("ALT \t calo  %3d (pt %7.2f)    ----> discarded\n", ic, calo.floatPt());
                calo.used = true;
                calo.setFloatPt(pt0); discardCalo(r, calo, 1);  // log this as discarded, for debugging
            }
        }
    }

    /// ------------- next step (needs the previous) ----------------
   
    // track to calo matching (first iteration, with a lower bound on the calo pt; there may be another one later)
    std::vector<int> tk2calo(r.track.size(), -1);
    for (int itk = 0, ntk = r.track.size(); itk < ntk; ++itk) {
        const auto & tk = r.track[itk]; 
        if (tk.muonLink || tk.used) continue; // not necessary but just a waste of CPU otherwise
        float drbest = drMatch_;
        float minCaloPt = tk.floatPt() - ptMatchLow_*tk.floatCaloPtErr();
        if (debug_) printf("ALT \t track %3d (pt %7.2f) to be matched to calo, min pT %7.2f\n", itk, tk.floatPt(), minCaloPt );
        for (int ic = 0, nc = r.calo.size(); ic < nc; ++ic) {
            auto & calo = r.calo[ic]; 
            if (calo.used || calo.floatPt() < minCaloPt) continue;
            float dr = floatDR(tk, calo);
            if (dr < drbest) { tk2calo[itk] = ic; drbest = dr; }
        }
        if (debug_ && tk2calo[itk] != -1) printf("ALT \t track %3d (pt %7.2f) matches to calo %3d (pt %7.2f) with dr %.3f\n", itk, tk.floatPt(), tk2calo[itk], tk2calo[itk] == -1 ? 0.0 : r.calo[tk2calo[itk]].floatPt(), drbest );
    }

    /// ------------- next step (needs the previous) ----------------

    // for each calo, compute the sum of the track pt
    std::vector<int> calo2ntk(r.calo.size(), 0);
    std::vector<float> calo2sumtkpt(r.calo.size(), 0);
    std::vector<float> calo2sumtkpterr(r.calo.size(), 0);
    for (int ic = 0, nc = r.calo.size(); ic < nc; ++ic) {
        const auto & calo = r.calo[ic];
        if (std::abs(calo.floatEta()) > 2.5) continue; 
        for (int itk = 0, ntk = r.track.size(); itk < ntk; ++itk) {
            if (tk2calo[itk] == ic) {
                const auto & tk = r.track[itk]; 
                if (tk.muonLink || tk.used) continue; 
                calo2ntk[ic]++;
                calo2sumtkpt[ic] += tk.floatPt();
                calo2sumtkpterr[ic] += std::pow(tk.floatCaloPtErr(), sumTkCaloErr2_ ? 2 : 1);
            }
        }
        if (sumTkCaloErr2_ && calo2sumtkpterr[ic] > 0)  calo2sumtkpterr[ic] = std::sqrt(calo2sumtkpterr[ic]);
    }


    // in the meantime, promote unlinked low pt tracks to hadrons
    for (int itk = 0, ntk = r.track.size(); itk < ntk; ++itk) {
        auto & tk = r.track[itk]; 
        if (tk2calo[itk] != -1 || tk.muonLink || tk.used) continue;
        if (tk.floatPt() < maxInvisiblePt_) {
            if (debug_) printf("ALT \t track %3d (pt %7.2f) not matched to calo, kept as charged hadron\n", itk, tk.floatPt());
            addTrackToPF(r, tk);
            tk.used = true;
        } else {
            if (debug_) printf("ALT \t track %3d (pt %7.2f) not matched to calo, dropped\n", itk, tk.floatPt());
            discardTrack(r, tk, 1); // log this as discarded, for debugging
        }
    }

    /// ------------- next step (needs the previous) ----------------

    /// OPTIONAL STEP: try to recover split hadron showers (v1.0): 
    //     take hadrons that are not track matched, close by a hadron which has an excess of track pt vs calo pt 
    //     add this pt to the calo pt of the other cluster 
    if (caloReLinkStep_) {
        std::vector<float> addtopt(r.calo.size(), 0);
        for (int ic = 0, nc = r.calo.size(); ic < nc; ++ic) {
            auto & calo = r.calo[ic];
            if (calo2ntk[ic] != 0 || calo.used || std::abs(calo.floatEta()) > 2.5) continue;
            int i2best = -1; float drbest = caloReLinkDr_;
            for (int ic2 = 0; ic2 < nc; ++ic2) {
                const auto & calo2 = r.calo[ic2];
                if (calo2ntk[ic2] == 0 || calo2.used || std::abs(calo2.floatEta()) > 2.5) continue;
                float dr = floatDR(calo,calo2);
                //// uncomment below for more verbose debugging
                //if (debug_ && dr < 0.5) printf("ALT \t calo  %3d (pt %7.2f) with no tracks is at dr %.3f from calo %3d with pt %7.2f (sum tk pt %7.2f), track excess %7.2f +- %7.2f\n", ic, calo.floatPt(), dr, ic2, calo2.floatPt(), calo2sumtkpt[ic2], calo2sumtkpt[ic2] - calo2.floatPt(), useTrackCaloSigma_ ? calo2sumtkpterr[ic2] : calo2.floatPtErr());
                if (dr < drbest) {
                    float ptdiff = calo2sumtkpt[ic2] - calo2.floatPt() + (useTrackCaloSigma_ ? calo2sumtkpterr[ic2] : calo2.floatPtErr());
                    if (ptdiff >= caloReLinkThreshold_*calo.floatPt()) {
                        i2best = ic2; drbest = dr;
                    }
                }
            }
            if (i2best != -1) {
                const auto & calo2 = r.calo[i2best];
                if (debug_) printf("ALT \t calo  %3d (pt %7.2f) with no tracks matched within dr %.3f with calo %3d with pt %7.2f (sum tk pt %7.2f), track excess %7.2f +- %7.2f\n", ic, calo.floatPt(), drbest, i2best, calo2.floatPt(), calo2sumtkpt[i2best], calo2sumtkpt[i2best] - calo2.floatPt(), useTrackCaloSigma_ ? calo2sumtkpterr[i2best] : calo2.floatPtErr());
                calo.used = true;
                addtopt[i2best] += calo.floatPt();
            }
        }
        // we do this at the end, so that the above loop is parallelizable
        for (int ic = 0, nc = r.calo.size(); ic < nc; ++ic) {
            if (addtopt[ic]) {
                auto & calo = r.calo[ic];
                if (debug_) printf("ALT \t calo  %3d (pt %7.2f, sum tk pt %7.2f) is increased to pt %7.2f after merging\n", ic, calo.floatPt(), calo2sumtkpt[ic], calo.floatPt() + addtopt[ic]);
                calo.setFloatPt(calo.floatPt() + addtopt[ic]);
            }
        }
    }
 


    /// ------------- next step (needs the previous) ----------------

    // process matched calo clusters, compare energy to sum track pt
    std::vector<float> calo2alpha(r.calo.size(), 1);
    for (int ic = 0, nc = r.calo.size(); ic < nc; ++ic) {
        auto & calo = r.calo[ic];
        if (calo2ntk[ic] == 0 || calo.used) continue;
        float ptdiff = calo.floatPt() - calo2sumtkpt[ic];
        float pterr  = useTrackCaloSigma_ ? calo2sumtkpterr[ic] : calo.floatPtErr();
        if (debug_) printf("ALT \t calo  %3d (pt %7.2f +- %7.2f, empt %7.2f) has %2d tracks (sumpt %7.2f, sumpterr %7.2f), ptdif %7.2f +- %7.2f\n", ic, calo.floatPt(), calo.floatPtErr(), calo.floatEmPt(), calo2ntk[ic], calo2sumtkpt[ic], calo2sumtkpterr[ic], ptdiff, pterr);
        if (ptdiff > +ptMatchHigh_*pterr) {
            if (ecalPriority_) {
                if (calo.floatEmPt() > 1) {
                    float emptdiff = std::min(ptdiff, calo.floatEmPt());
                    if (debug_) printf("ALT \t calo  %3d (pt %7.2f, empt %7.2f)    ---> make photon with pt %7.2f, reduce ptdiff to %7.2f +- %7.2f\n", ic, calo.floatPt(), calo.floatEmPt(), emptdiff, ptdiff-emptdiff, pterr);
                    auto & p = addCaloToPF(r, calo);
                    p.setFloatPt(emptdiff);
                    p.hwId = l1tpf::Particle::GAMMA;
                    ptdiff -= emptdiff;
                }
                if (ptdiff > 2) {
                    if (debug_) printf("ALT \t calo  %3d (pt %7.2f, empt %7.2f)    ---> make also neutral hadron with pt %7.2f\n", ic, calo.floatPt(), calo.floatEmPt(), ptdiff);
                    auto & p = addCaloToPF(r, calo);
                    p.setFloatPt(ptdiff);
                    p.hwId = l1tpf::Particle::NH;
                }
            } else {
                if (debug_) printf("ALT \t calo  %3d (pt %7.2f)    ---> promoted to neutral with pt %7.2f\n", ic, calo.floatPt(), ptdiff);
                auto & p = addCaloToPF(r, calo);
                p.setFloatPt(ptdiff);
                calo.hwFlags = 0;
            }
        } else if (ptdiff > -ptMatchLow_*pterr) {
            // nothing to do (weighted average happens when we process the tracks)
            calo.hwFlags = 1;
            if (debug_) printf("ALT \t calo  %3d (pt %7.2f)    ---> to be deleted, will use tracks instead\n", ic, calo.floatPt());
            discardCalo(r, calo, 0); // log this as discarded, for debugging
        } else {
            // tracks overshoot, rescale to tracks to calo
            calo2alpha[ic] = calo.floatPt() / calo2sumtkpt[ic];
            calo.hwFlags = 2;
            if (debug_) printf("ALT \t calo  %3d (pt %7.2f)    ---> tracks overshoot and will be scaled down by %.4f\n", ic, calo.floatPt(), calo2alpha[ic]);
        }
        calo.used = true;
    }

    /// ------------- next step (needs the previous) ----------------

    // process matched tracks, if necessary rescale or average
    for (int itk = 0, ntk = r.track.size(); itk < ntk; ++itk) {
        auto & tk = r.track[itk]; 
        if (tk2calo[itk] == -1 || tk.muonLink || tk.used) continue;
        auto & p = addTrackToPF(r, tk);
        tk.used = true;
        const auto & calo = r.calo[tk2calo[itk]];
        if (calo.hwFlags == 1) {
            // can do weighted average if there's just one track
            if (calo2ntk[tk2calo[itk]] == 1) { 
                float ptavg = tk.floatPt();
                if (tk.floatPtErr() > 0) {
                    float wcalo = 1.0/std::pow(tk.floatCaloPtErr(), 2);
                    float wtk   = 1.0/std::pow(tk.floatPtErr(), 2);
                    ptavg = ( calo.floatPt() * wcalo + tk.floatPt() * wtk ) / (wcalo + wtk );
                }
                p.setFloatPt(ptavg);
                if (debug_) printf("ALT \t track %3d (pt %7.2f +- %7.2f) combined with calo %3d (pt %7.2f +- %7.2f (from tk) yielding candidate of pt %7.2f\n", 
                                    itk, tk.floatPt(), tk.floatPtErr(), tk2calo[itk], calo.floatPt(), tk.floatCaloPtErr(), ptavg );
            } else {
                if (debug_) printf("ALT \t track %3d (pt %7.2f) linked to calo %3d promoted to charged hadron\n", itk, tk.floatPt(), tk2calo[itk]);
            }
        } else if (calo.hwFlags == 2) {
            // must rescale
            p.setFloatPt(tk.floatPt() * calo2alpha[tk2calo[itk]]); 
            if (debug_) printf("ALT \t track %3d (pt %7.2f) linked to calo %3d promoted to charged hadron with pt %7.2f after rescaling\n", itk, tk.floatPt(), tk2calo[itk], p.floatPt());
        }
    }

    // process unmatched calo clusters 
    for (int ic = 0, nc = r.calo.size(); ic < nc; ++ic) {
        if (!r.calo[ic].used) {
            addCaloToPF(r, r.calo[ic]);
            if (debug_) printf("ALT \t calo  %3d (pt %7.2f) not linked, promoted to neutral\n", ic, r.calo[ic].floatPt());
        }
    }
    
    // finally do muons
    if(!skipMuons_) { 
        for (int itk = 0, ntk = r.track.size(); itk < ntk; ++itk) {
            if (r.track[itk].muonLink) addTrackToPF(r, r.track[itk]);
        }
    }
}

