#include <cmath>
#include <cstdio>
#include <cassert>
#include <vector>
#include <TH1F.h>
#include <TEfficiency.h>
#include <TGraph.h>
#include <Math/PtEtaPhiM4D.h>
#include <Math/LorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include "L1Trigger/Phase2L1ParticleFlow/src/corrector.h"

unsigned int fillTH1FastGenCut(TH1F *ret, const std::vector<float> & corrArr, const std::vector<float> & genArr, float genThr) {
    assert(genArr.size() == corrArr.size());
    unsigned int nsel = 0;
    for (unsigned int i = 0, n = corrArr.size(); i < n; ++i) {
        if (genArr[i] > genThr) { ret->Fill(corrArr[i]); nsel++; }
    }
    unsigned int nbins = ret->GetNbinsX();
    ret->SetBinContent(nbins, ret->GetBinContent(nbins) + ret->GetBinContent(nbins+1));
    ret->GetBinContent(nbins+1, 0);
    return nsel;
}

unsigned int fillTH1Fast(TH1F *ret, const std::vector<float> & corrArr) {
    for (unsigned int i = 0, n = corrArr.size(); i < n; ++i) {
        ret->Fill(corrArr[i]);
    }
    unsigned int nbins = ret->GetNbinsX();
    ret->SetBinContent(nbins, ret->GetBinContent(nbins) + ret->GetBinContent(nbins+1));
    ret->GetBinContent(nbins+1, 0);
    return corrArr.size();
}

void fillTEffFast(TEfficiency *eff, const std::vector<float> & refArr, const std::vector<float> & corrArr, float corrThr) {
    assert(refArr.size() == corrArr.size());
    for (unsigned int i = 0, n = refArr.size(); i < n; ++i) {
        eff->Fill(corrArr[i] > corrThr, refArr[i]);
    }
}

TGraph *makeROCFast(TH1 *effsig, TH1 *effbkg) {
    TGraph *graph = new TGraph(effsig->GetNbinsX());
    for (unsigned int i = 1, n = graph->GetN()+1; i < n; ++i) {
        graph->SetPoint(i-1, effsig->GetBinContent(i), effbkg->GetBinContent(i));
    }
    return graph;
}

class JetCalcBase {
    public:
        typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> Jet;
        JetCalcBase() {}
        virtual ~JetCalcBase() {}
        virtual float operator()(const std::vector<Jet> & jets) const = 0;
};
class CalcHT : public JetCalcBase {
    public:
        CalcHT() : JetCalcBase()  {}
        virtual float operator()(const std::vector<Jet> & jets) const {
            float ret = 0;
            for (const auto & j : jets) ret += j.pt();
            return ret;
        }
};
class CalcMHT : public JetCalcBase {
    public:
        CalcMHT() : JetCalcBase()  {}
        virtual float operator()(const std::vector<Jet> & jets) const {
            float retx = 0, rety = 0;
            for (const auto & j : jets) {
                retx += j.pt() * std::cos(j.phi());
                rety += j.pt() * std::sin(j.phi());
            }
            return std::hypot(retx,rety);
        }
};
class CalcJ : public JetCalcBase {
    public:
        CalcJ(unsigned int n) : JetCalcBase(), n_(n)  {}
        virtual float operator()(const std::vector<Jet> & jets) const {
            if (jets.size() < n_) return 0.;
            pts.clear();
            for (const auto & j : jets) pts.push_back(-j.pt());
            std::sort(pts.begin(), pts.end());
            return -pts[n_-1];
        }
    protected:
        unsigned int n_;
        mutable std::vector<float> pts;
};
class CalcMJJ : public JetCalcBase {
    public:
        CalcMJJ() : JetCalcBase() {}
        virtual float operator()(const std::vector<Jet> & jets) const {
            if (jets.size() < 2) return 0.;
            float ret = 0;
            for (unsigned int i1 = 0, n = jets.size(); i1 < n; ++i1) {
                for (unsigned int i2 = i1+1; i2 < n; ++i2) {
                    ret = std::max(ret, (jets[i1]+jets[i2]).M());
                }
            }
            return ret;
        }
};
class CalcJ2_MJJcut : public JetCalcBase {
    public:
        CalcJ2_MJJcut(float mjj) : JetCalcBase(), mjj_(mjj)  {}
        virtual float operator()(const std::vector<Jet> & jets) const {
            if (jets.size() < 2) return 0.;
            float ret = 0;
            for (unsigned int i1 = 0, n = jets.size(); i1 < n; ++i1) {
                for (unsigned int i2 = i1+1; i2 < n; ++i2) {
                    if ((jets[i1]+jets[i2]).M() > mjj_) {
                        ret = std::max(ret, std::min(jets[i1].pt(), jets[i2].pt()));
                    }
                }
            }
            return ret;
        }
    protected:
        float mjj_;
};




std::vector<float> makeJetArray(TTree *tree, const std::string & obj, float ptCut, float etaCut, const JetCalcBase &calc, const l1tpf::corrector *jetcorr = nullptr) {
    std::vector<float> ret; 
    ret.reserve(tree->GetEntries());

    TTreeReader reader(tree);
    TTreeReaderArray<float> jet_pt(reader, (obj+"Jets_pt").c_str()), jet_eta(reader, (obj+"Jets_eta").c_str()), jet_phi(reader, (obj+"Jets_phi").c_str());

    //printf("Reading from tree %s (%llu entries) the %sJets, ptCut %g, etaCut %g\n", tree->GetDirectory()->GetFile()->GetName(), tree->GetEntries(), obj.c_str(), ptCut, etaCut);
    std::vector<JetCalcBase::Jet> jets;
    //unsigned int iev = 0;
    while (reader.Next()) {
        jets.clear();
        unsigned int njets = jet_pt.GetSize();
        for (unsigned int i = 0; i < njets; ++i) {
            //if (iev < 10) printf("        ev %3u   rawjet %2u: %10.3f %+5.2f %+5.2f\n", iev, i, jet_pt[i], jet_eta[i], jet_phi[i]);
            float pt = jet_pt[i], eta = jet_eta[i];
            if (std::abs(eta) > etaCut) continue;
            if (jetcorr) pt = jetcorr->correctedPt(pt, eta);
            if (pt > ptCut) {
                jets.emplace_back(pt, eta, jet_phi[i], 0.); // jet_mass[i]);
                //if (iev < 10) printf("        ev %3u      jet %2u: %10.3f %+5.2f %+5.2f\n", iev, i, pt, eta, jet_phi[i]);
            }
        }
        ret.push_back(calc(jets));
        //if (iev < 10) printf("        ev %3u calc = %10.3f\n", iev, ret.back());
        //iev++;
    }
    return ret;
}

std::vector<float> makeMetArray(TTree *tree, const std::string & obj) {
    std::vector<float> ret; 
    ret.reserve(tree->GetEntries());

    TTreeReader reader(tree);
    TTreeReaderValue<float> met_pt(reader, (obj+"_pt").c_str());

    while (reader.Next()) {
        ret.push_back(*met_pt);
    }
    return ret;
}

std::vector<float> makeMinimum(const std::vector<float> & v1, const std::vector<float> & v2) {
    std::vector<float> ret(v1.size()); 
    for (unsigned int i = 0, n = v1.size(); i < n; ++i) {
        ret[i] = std::min(v1[i], v2[i]);
    }
    return ret;
}
