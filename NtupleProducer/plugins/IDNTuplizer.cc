// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TPFUtils.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include <cstdint>
#include <TTree.h>

class IDNTuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns>  {
    public:
        explicit IDNTuplizer(const edm::ParameterSet&);
        ~IDNTuplizer();

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

        virtual void beginRun(edm::Run const&, edm::EventSetup const& iSetup) override {
            //edm::ESHandle<MagneticField> magneticField;
            //iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
            //bZ_ = magneticField->inTesla(GlobalPoint(0,0,0)).z();
            bZ_ = 3.8112; // avoid loading the event setup
        }
        virtual void endRun(edm::Run const&, edm::EventSetup const& iSetup) override { } // framework wants this to be implemented

        edm::EDGetTokenT<reco::CandidateView> src_;
        StringCutObjectSelector<reco::Candidate> sel_;
        
        edm::EDGetTokenT<std::vector<reco::GenParticle>> genparticles_;
        bool prop_;
        float dr2Max_, minPtRatio_;
        bool onlyMatched_;
        TTree *tree_;
        uint32_t run_, lumi_; uint64_t event_;

        struct McVars {
            float pt, eta, phi, dr;
            int   id;
            void makeBranches(TTree *tree) {
                tree->Branch("genid", &id, "genid/I");
                tree->Branch("gendr", &dr, "gendr/F");
                tree->Branch("genpt", &pt, "genpt/F");
                tree->Branch("geneta", &eta, "geneta/F");
                tree->Branch("genphi", &phi, "genphi/F");
            }
            void clear() {
                pt = 0; eta = 0; phi = 0; dr = 0; id = 0;
            }
            void fill(const reco::GenParticle &c, bool propagate, float bz) {
                id = c.pdgId();
                pt = c.pt(); 
                if (propagate && c.charge() != 0) {
                    math::XYZTLorentzVector vertex(c.vx(),c.vy(),c.vz(),0.);
                    auto caloetaphi = l1tpf::propagateToCalo(c.p4(),vertex,c.charge(),bz);
                    eta = caloetaphi.first; phi = caloetaphi.second;
                } else {
                    eta = c.eta(); phi = c.phi();
                }
            }
        } mc_;

        class RecoVar {
            public:
                RecoVar(const std::string & name, const std::string & expr) : name_(name), expr_(expr,true) {}
                void makeBranch(TTree *tree) {
                    tree->Branch(name_.c_str(), &val_, (name_+"/F").c_str());
                }
                void fill(const reco::Candidate & c) {
                    val_ = expr_(c);
                }
            private:
                std::string name_;
                StringObjectFunction<reco::Candidate> expr_;
                float val_;
        };
        std::vector<RecoVar> reco_;

        struct ExtVar { 
                ExtVar(const std::string & name, edm::EDGetTokenT<edm::ValueMap<float>> token) : name_(name), token_(token) {}
                void makeBranch(TTree *tree) {
                    tree->Branch(name_.c_str(), &val_, (name_+"/F").c_str());
                }
                void init(const edm::Event & iEvent) {
                    iEvent.getByToken(token_, handle_);
                }
                void fill(const reco::CandidatePtr & c) {
                    val_ = (*handle_)[c];
                }
            private:
                std::string name_;
                edm::EDGetTokenT<edm::ValueMap<float>> token_;
                edm::Handle<edm::ValueMap<float>> handle_;
                float val_;
        };
        std::vector<ExtVar> ext_;

        float bZ_;
 
};

IDNTuplizer::IDNTuplizer(const edm::ParameterSet& iConfig) :
    src_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("src"))),
    sel_(iConfig.getParameter<std::string>("cut"), true),
    genparticles_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    prop_(iConfig.getParameter<bool>("propagateToCalo")),
    dr2Max_(std::pow(iConfig.getParameter<double>("drMax"), 2)),
    minPtRatio_(iConfig.getParameter<double>("minRecoPtOverGenPt")),
    onlyMatched_(iConfig.getParameter<bool>("onlyMatched"))
{
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("tree","tree");
    tree_->Branch("run",  &run_,  "run/i");
    tree_->Branch("lumi", &lumi_, "lumi/i");
    tree_->Branch("event", &event_, "event/l");

    edm::ParameterSet vars = iConfig.getParameter<edm::ParameterSet>("variables");
    auto reconames = vars.getParameterNamesForType<std::string>();
    for (const std::string & name : reconames) {
        reco_.emplace_back(name, vars.getParameter<std::string>(name));
    }

    if (iConfig.existsAs<edm::ParameterSet>("extVariables")) {
        edm::ParameterSet evars = iConfig.getParameter<edm::ParameterSet>("extVariables");
        auto enames = evars.getParameterNamesForType<edm::InputTag>();
        for (const std::string & name : enames) {
            ext_.emplace_back(name, consumes<edm::ValueMap<float>>(evars.getParameter<edm::InputTag>(name)));
        }

    }
}

IDNTuplizer::~IDNTuplizer() { }

// ------------ method called once each job just before starting event loop  ------------
    void 
IDNTuplizer::beginJob()
{
    mc_.makeBranches(tree_);
    for (auto & v : reco_) v.makeBranch(tree_);
    for (auto & v : ext_) v.makeBranch(tree_);
}


// ------------ method called for each event  ------------
    void
IDNTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    run_  = iEvent.id().run();
    lumi_ = iEvent.id().luminosityBlock();
    event_ = iEvent.id().event();

    edm::Handle<std::vector<reco::GenParticle>> genparticles;
    iEvent.getByToken(genparticles_, genparticles);

    edm::Handle<reco::CandidateView> src;
    iEvent.getByToken(src_, src);

    for (auto & v : ext_) v.init(iEvent);

    std::vector<reco::CandidatePtr> selected;
    for (unsigned int i = 0, n = src->size(); i < n; ++i) {
        const auto & c = (*src)[i];
        if (sel_(c)) selected.push_back(src->ptrAt(i));
    }

    std::vector<bool> matched(selected.size(), 0);

    for (const reco::GenParticle &gen : *genparticles) {
        mc_.fill(gen, prop_, bZ_);
        float dr2best = dr2Max_; int ibest = -1;
        for (unsigned int ireco = 0, nreco = selected.size(); ireco < nreco; ++ireco) {
            const reco::Candidate & reco = *(selected[ireco]);
            if (reco.pt() <= minPtRatio_*gen.pt()) continue;
            float dr2 = deltaR2(reco.eta(), reco.phi(), mc_.eta, mc_.phi);
            if (dr2 < dr2best) { dr2best = dr2; ibest = ireco; }
        }
        if (ibest != -1) {
            matched[ibest] = true;
            mc_.dr = std::sqrt(dr2best);
            const reco::CandidatePtr & reco = (selected[ibest]);
            for (auto & v : reco_) v.fill(*reco);
            for (auto & v : ext_) v.fill(reco);
            tree_->Fill();
        }
    }
    if (!onlyMatched_) {
        mc_.clear();
        for (unsigned int ireco = 0, nreco = selected.size(); ireco < nreco; ++ireco) {
            if (matched[ireco]) continue;
            const reco::CandidatePtr & reco = (selected[ireco]);
            for (auto & v : reco_) v.fill(*reco);
            for (auto & v : ext_) v.fill(reco);
            tree_->Fill();
        }
    }
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(IDNTuplizer);
