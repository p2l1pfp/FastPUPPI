// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include <algorithm>

class L1PFJetTableProducer : public edm::global::EDProducer<>  {
    public:
        explicit L1PFJetTableProducer(const edm::ParameterSet&);
        ~L1PFJetTableProducer();

    private:
        virtual void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

        edm::EDGetTokenT<reco::CandidateView> gen_;
        StringCutObjectSelector<reco::Candidate> sel_;

        struct ExtraVar {
            std::string name, expr;
            StringObjectFunction<reco::Candidate> func;
            ExtraVar(const std::string & n, const std::string & expr) : name(n), expr(expr), func(expr, true) {}
        };
        std::vector<ExtraVar> extraVars_;

        struct JetRecord {
            public:
                std::string coll;
                edm::EDGetTokenT<reco::CandidateView> src;
                StringCutObjectSelector<reco::Candidate> sel;
                
                JetRecord(const std::string & name, const edm::EDGetTokenT<reco::CandidateView> & tag, const edm::ParameterSet & pset) :
                    coll(name), src(tag), 
                    sel(pset.existsAs<std::string>(name+"_sel") ? pset.getParameter<std::string>(name+"_sel") : "", true) {}
        };
        std::vector<JetRecord> jets_;

        float minPt_, dr2Max_, minPtRatio_;
};

L1PFJetTableProducer::L1PFJetTableProducer(const edm::ParameterSet& iConfig) :
    gen_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("gen"))),
    sel_(iConfig.getParameter<std::string>("commonSel"), true),
    dr2Max_(std::pow(iConfig.getParameter<double>("drMax"), 2)),
    minPtRatio_(iConfig.getParameter<double>("minRecoPtOverGenPt"))
{
    edm::ParameterSet jets = iConfig.getParameter<edm::ParameterSet>("jets");
    auto jetnames = jets.getParameterNamesForType<edm::InputTag>();
    for (const std::string & name : jetnames) {
        jets_.emplace_back(name, consumes<reco::CandidateView>(jets.getParameter<edm::InputTag>(name)), jets);
        produces<nanoaod::FlatTable>(name+"Jets");
    }

    if (iConfig.existsAs<edm::ParameterSet>("moreVariables")) {
        edm::ParameterSet vars = iConfig.getParameter<edm::ParameterSet>("moreVariables");
        auto morenames = vars.getParameterNamesForType<std::string>();
        for (const std::string & name : morenames) {
            extraVars_.emplace_back(name, vars.getParameter<std::string>(name));
        }
    }
 }

L1PFJetTableProducer::~L1PFJetTableProducer() { }

// ------------ method called for each event  ------------
    void
L1PFJetTableProducer::produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{
    edm::Handle<reco::CandidateView> gens;
    iEvent.getByToken(gen_, gens);

    edm::Handle<reco::CandidateView> src;
    std::vector<const reco::Candidate *> selected, matched;
    std::vector<float> vals_pt, vals_eta, vals_phi, vals_mass;
    for (auto & jets : jets_) {
        // get and select
        iEvent.getByToken(jets.src, src);
        for (const auto & j : *src) {
            if (sel_(j) && jets.sel(j)) {
                selected.push_back(&j);
            }
        }
        // mc-match
        matched.resize(selected.size());
        std::fill(matched.begin(), matched.end(), nullptr);
        std::vector<std::pair<float,std::pair<const reco::Candidate*, unsigned int>>> matches;
        for (const auto & genj : *gens) {
            for (unsigned int ireco = 0, nreco = selected.size(); ireco < nreco; ++ireco) {
                const reco::Candidate & reco = *(selected[ireco]);
                if (reco.pt() <= minPtRatio_*genj.pt()) continue;
                float dr2 = deltaR2(reco, genj);
                if (dr2 < dr2Max_) matches.emplace_back(dr2, std::make_pair(&genj, ireco));
            }
        }
        std::sort(matches.begin(), matches.end());
        for (unsigned int im = 0, nm = matches.size(); im < nm; ++im) {
            const auto & match = matches[im];
            if (match.second.first == nullptr) continue;
            // set the match
            matched[match.second.second] = match.second.first;
            // remove any other candidate pair that overlaps with this one
            for (unsigned int im2 = im+1; im2 < nm; ++im2) {
                auto & match2 = matches[im2];
                if (match2.second.first == nullptr) continue;
                if (match2.second.first == match.second.first ||
                    match2.second.second == match.second.second) {
                    match2.second.first = nullptr;
                }
            }
        }
        
        // create the table
        unsigned int njets = selected.size();
        auto out = std::make_unique<nanoaod::FlatTable>(njets, jets.coll+"Jets", false);

        // fill basic info
        vals_pt.resize(njets); 
        vals_eta.resize(njets); 
        vals_phi.resize(njets); 
        vals_mass.resize(njets); 
        for (unsigned int i = 0; i < njets; ++i) {
            vals_pt[i] = selected[i]->pt();
            vals_eta[i] = selected[i]->eta();
            vals_phi[i] = selected[i]->phi();
            vals_mass[i] = selected[i]->mass();
        }
        out->addColumn<float>("pt", vals_pt, "pt of jet", nanoaod::FlatTable::FloatColumn);
        out->addColumn<float>("eta", vals_eta, "eta of jet", nanoaod::FlatTable::FloatColumn);
        out->addColumn<float>("phi", vals_phi, "phi of jet", nanoaod::FlatTable::FloatColumn);
        out->addColumn<float>("mass", vals_mass, "mass of jet", nanoaod::FlatTable::FloatColumn);

        // fill gen match (if this is not the gen selection)
        if (src.id() != gens.id()) {
            for (unsigned int i = 0; i < njets; ++i) {
                if (matched[i] != nullptr) {
                    vals_pt[i] = matched[i]->pt();
                    vals_eta[i] = deltaR(*selected[i], *matched[i]);;
                } else {
                    vals_pt[i] = 0;
                    vals_eta[i] = 99.9;
                }
            }
            out->addColumn<float>("genpt", vals_pt, "pt of matched gen jet", nanoaod::FlatTable::FloatColumn);
            out->addColumn<float>("gendr", vals_eta, "dr of matched gen jet", nanoaod::FlatTable::FloatColumn);
        }

        // fill extra vars
        for (const auto & evar : extraVars_) {
            for (unsigned int i = 0; i < njets; ++i) {
                vals_pt[i] = evar.func(*selected[i]);
            }
            out->addColumn<float>(evar.name, vals_pt, evar.expr, nanoaod::FlatTable::FloatColumn);
        }
        
        // save to the event branches
        iEvent.put(std::move(out), jets.coll+"Jets");

        // clear
        selected.clear();
        matched.clear();
    }
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PFJetTableProducer);
