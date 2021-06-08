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
std::pair<unsigned,unsigned> inclusiveEffFastGenCut(const std::vector<float> & refArr, float refThr, const std::vector<float> & corrArr, float corrThr) {
    assert(refArr.size() == corrArr.size());
    std::pair<unsigned,unsigned> ret = std::pair<unsigned,unsigned>(0,0);
    for (unsigned int i = 0, n = refArr.size(); i < n; ++i) {
        if (refArr[i] > refThr) {
            ret.second++;
            if (corrArr[i] > corrThr) ret.first++;
        }
    }
    return ret;
}
std::pair<unsigned,unsigned> inclusiveEffFast(const std::vector<float> & corrArr, float corrThr) {
    std::pair<unsigned,unsigned> ret = std::pair<unsigned,unsigned>(0,0);
    for (unsigned int i = 0, n = corrArr.size(); i < n; ++i) {
        ret.second++;
        if (corrArr[i] > corrThr) ret.first++;
    }
    return ret;
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


class LepCalcBase {
    public:
        class Lep : public ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> {
            public:
                Lep(float pt, float eta, float phi, int charge, float iso, float id) :
                     ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>(pt,eta,phi,0.),
                     charge_(charge), iso_(iso), id_(id) {}
                int charge() const { return charge_; }
                float iso() const { return iso_; }
                float id() const { return id_; }
            private:
                int charge_;
                float id_, iso_;
        };
        LepCalcBase() {}
        virtual ~LepCalcBase() {}
        virtual float operator()(const std::vector<Lep> & leps) const = 0;
};
class CalcL : public LepCalcBase {
    public:
        enum RetType { RetPt=0, RetIso=1, RetId=2 };
        CalcL(unsigned int n=1, RetType ret = RetPt) : 
            LepCalcBase(), n_(n), ret_(ret) {} 
        virtual float operator()(const std::vector<Lep> & leps) const {
            if (leps.size() < n_) return 0.;
            pts.clear();
            switch(ret_) {
                case RetPt: 
                    for (const auto & l : leps) pts.push_back(-l.pt());
                    break;
                case RetIso: 
                    for (const auto & l : leps) pts.push_back(-l.iso());
                    break;
                case RetId: 
                    for (const auto & l : leps) pts.push_back(l.id());
                    break;
            }
            std::sort(pts.begin(), pts.end());
            return (ret_ == RetId ? +1 : -1)*pts[n_-1];
        }
    protected:
        unsigned int n_; RetType ret_;
        mutable std::vector<float> pts;
};


std::vector<float> makeLepArray(TTree *tree, const std::string & obj, float ptCut, float etaCut, const LepCalcBase &calc, 
                                const std::string & idName, float idCut, bool idIsFloat, const std::string & isoName, float isoCut, bool makeIsoRel) {
    std::vector<float> ret; 
    ret.reserve(tree->GetEntries());

    bool hasId = !idName.empty(), hasIso = !isoName.empty();
    std::string idf = (hasId && idIsFloat ? idName : "pt"), idi = (hasId && !idIsFloat ? idName : "charge");
    std::string isof = hasIso ? isoName : "pt";
    TTreeReader reader(tree);
    TTreeReaderArray<float> lep_pt(reader, (obj+"_pt").c_str()), lep_eta(reader, (obj+"_eta").c_str()), lep_phi(reader, (obj+"_phi").c_str());
    TTreeReaderArray<int> lep_charge(reader, (obj+"_charge").c_str()), lep_idi(reader, (obj+"_"+idi).c_str());
    TTreeReaderArray<float> lep_iso(reader, (obj+"_"+isof).c_str()), lep_idf(reader, (obj+"_"+idf).c_str());

    std::vector<LepCalcBase::Lep> leps;
    //unsigned int iev = 0;
    while (reader.Next()) {
        leps.clear();
        unsigned int nleps = lep_pt.GetSize();
        for (unsigned int i = 0; i < nleps; ++i) {
            float pt = lep_pt[i], eta = lep_eta[i];
            if (pt < ptCut || std::abs(eta) > etaCut) continue;
            float id = hasId ? (idIsFloat ? lep_idf[i] : lep_idi[i]) : 1, iso = hasIso ? lep_iso[i] : 0;
            if (makeIsoRel) iso /= pt;
            if (hasId && id < idCut) continue;
            if (hasIso && iso > isoCut) continue;
            leps.emplace_back(pt, eta, lep_phi[i], lep_charge[i], iso, id); 
        }
        ret.push_back(calc(leps));
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
