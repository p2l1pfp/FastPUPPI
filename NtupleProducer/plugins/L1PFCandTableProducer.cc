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

class L1PFCandTableProducer : public edm::global::EDProducer<>  {
    public:
        explicit L1PFCandTableProducer(const edm::ParameterSet&);
        ~L1PFCandTableProducer();

    private:
        virtual void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

        StringCutObjectSelector<reco::Candidate> sel_;

        struct ExtraVar {
            std::string name, expr;
            StringObjectFunction<reco::Candidate> func;
            ExtraVar(const std::string & n, const std::string & expr) : name(n), expr(expr), func(expr, true) {}
        };
        std::vector<ExtraVar> extraVars_;

        struct CandRecord {
            public:
                std::string coll;
                edm::EDGetTokenT<reco::CandidateView> src;
                StringCutObjectSelector<reco::Candidate> sel;
                
                CandRecord(const std::string & name, const edm::EDGetTokenT<reco::CandidateView> & tag, const edm::ParameterSet & pset) :
                    coll(name), src(tag), 
                    sel(pset.existsAs<std::string>(name+"_sel") ? pset.getParameter<std::string>(name+"_sel") : "", true) {}
        };
        std::vector<CandRecord> cands_;
};

L1PFCandTableProducer::L1PFCandTableProducer(const edm::ParameterSet& iConfig) :
    sel_(iConfig.getParameter<std::string>("commonSel"), true)
{
    edm::ParameterSet cands = iConfig.getParameter<edm::ParameterSet>("cands");
    auto candnames = cands.getParameterNamesForType<edm::InputTag>();
    for (const std::string & name : candnames) {
        cands_.emplace_back(name, consumes<reco::CandidateView>(cands.getParameter<edm::InputTag>(name)), cands);
        produces<nanoaod::FlatTable>(name+"Cands");
    }

    if (iConfig.existsAs<edm::ParameterSet>("moreVariables")) {
        edm::ParameterSet vars = iConfig.getParameter<edm::ParameterSet>("moreVariables");
        auto morenames = vars.getParameterNamesForType<std::string>();
        for (const std::string & name : morenames) {
            extraVars_.emplace_back(name, vars.getParameter<std::string>(name));
        }
    }
 }

L1PFCandTableProducer::~L1PFCandTableProducer() { }

// ------------ method called for each event  ------------
    void
L1PFCandTableProducer::produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{
    edm::Handle<reco::CandidateView> src;
    std::vector<const reco::Candidate *> selected;
    std::vector<float> vals_pt, vals_eta, vals_phi, vals_mass;
    for (auto & cands : cands_) {
        // get and select
        iEvent.getByToken(cands.src, src);
        for (const auto & j : *src) {
            if (sel_(j) && cands.sel(j)) {
                selected.push_back(&j);
            }
        }
        
        // create the table
        unsigned int ncands = selected.size();
        auto out = std::make_unique<nanoaod::FlatTable>(ncands, cands.coll+"Cands", false);

        // fill basic info
        vals_pt.resize(ncands); 
        vals_eta.resize(ncands); 
        vals_phi.resize(ncands); 
        vals_mass.resize(ncands); 
        for (unsigned int i = 0; i < ncands; ++i) {
            vals_pt[i] = selected[i]->pt();
            vals_eta[i] = selected[i]->eta();
            vals_phi[i] = selected[i]->phi();
            vals_mass[i] = selected[i]->mass();
        }
        out->addColumn<float>("pt", vals_pt, "pt of cand", nanoaod::FlatTable::FloatColumn);
        out->addColumn<float>("eta", vals_eta, "eta of cand", nanoaod::FlatTable::FloatColumn);
        out->addColumn<float>("phi", vals_phi, "phi of cand", nanoaod::FlatTable::FloatColumn);
        out->addColumn<float>("mass", vals_mass, "mass of cand", nanoaod::FlatTable::FloatColumn);

        // fill extra vars
        for (const auto & evar : extraVars_) {
            for (unsigned int i = 0; i < ncands; ++i) {
                vals_pt[i] = evar.func(*selected[i]);
            }
            out->addColumn<float>(evar.name, vals_pt, evar.expr, nanoaod::FlatTable::FloatColumn);
        }
        
        // save to the event branches
        iEvent.put(std::move(out), cands.coll+"Cands");

        // clear
        selected.clear();
    }
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PFCandTableProducer);
