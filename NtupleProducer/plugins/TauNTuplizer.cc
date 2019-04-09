// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Phase2L1ParticleFlow/interface/PFTau.h"
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
//#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include <cstdint>
#include <TTree.h>
#include <TLorentzVector.h>

class TauNTuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns>  {
    public:
        explicit TauNTuplizer(const edm::ParameterSet&);
        ~TauNTuplizer();

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void beginRun(edm::Run const&, edm::EventSetup const& iSetup) override {}
        virtual void endRun(edm::Run const&, edm::EventSetup const& iSetup) override { } // framework wants this to be implemented

        edm::EDGetTokenT<std::vector<l1t::PFTau>> L1PFTaus_;
        edm::EDGetTokenT<std::vector<reco::GenParticle>> genparticles_;
        float dr2Max_, minPtRatio_;
        TTree *tree_;
        uint32_t run_, lumi_; uint64_t event_; uint64_t eventcount_;

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
                pt = 0; eta = 0; phi = 0; dr = -999; id = 0;
            }
	  void fill(const reco::GenParticle &c) { 
            TLorentzVector lVec; 
	    lVec.SetPtEtaPhiM(pt,eta,phi,0);
	    TLorentzVector lNewVec;
	    lNewVec.SetPtEtaPhiM(c.pt(),c.eta(),c.phi(),0);
	    lVec += lNewVec;
	    id  = c.pdgId();
	    pt  = lVec.Pt(); 
	    eta = lVec.Eta(); 
	    phi = lVec.Phi();
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

 
};

TauNTuplizer::TauNTuplizer(const edm::ParameterSet& iConfig) :
  L1PFTaus_    (consumes<std::vector<l1t::PFTau>>(iConfig.getParameter<edm::InputTag>("src"))),
  genparticles_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  dr2Max_(std::pow(iConfig.getParameter<double>("drMax"), 2)),
  minPtRatio_(float(iConfig.getParameter<double>("minRecoPtOverGenPt"))) {
 
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("tree","tree");
    tree_->Branch("run",  &run_,  "run/i");
    tree_->Branch("lumi", &lumi_, "lumi/i");
    tree_->Branch("event", &event_, "event/l");
    tree_->Branch("event2", &eventcount_, "eventcount/l");

    edm::ParameterSet vars = iConfig.getParameter<edm::ParameterSet>("variables");
    auto reconames = vars.getParameterNamesForType<std::string>();
    for (const std::string & name : reconames) {
        reco_.emplace_back(name, vars.getParameter<std::string>(name));
    }
}

TauNTuplizer::~TauNTuplizer() { }
void  TauNTuplizer::beginJob() {
    mc_.makeBranches(tree_);
    for (auto & v : reco_) v.makeBranch(tree_);
    eventcount_=0;
}
void TauNTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    run_  = iEvent.id().run();
    lumi_ = iEvent.id().luminosityBlock();
    event_ = iEvent.id().event();
    eventcount_++;
    
    edm::Handle<std::vector<reco::GenParticle>> genparticles;
    iEvent.getByToken(genparticles_, genparticles);

    edm::Handle<  l1t::PFTauCollection > l1PFTaus;
    try { 
      iEvent.getByToken( L1PFTaus_, l1PFTaus);
    } catch(...) { 
      return;
    }
    l1t::PFTau dummy;
    std::vector<int> matchedGen;
    for (unsigned int i = 0, n = l1PFTaus->size(); i < n; ++i) {
        const auto & c = (*l1PFTaus)[i];
	float dr2best = dr2Max_; 
	int index = 0;
	mc_.clear();
	for (const reco::GenParticle &gen : *genparticles) {
	  index++;
	  if(!gen.isDirectPromptTauDecayProductFinalState()) continue;
	  if(fabs(gen.pdgId()) > 10 && fabs(gen.pdgId()) < 17) continue;
	  if (c.pt() <= minPtRatio_*gen.pt()) continue;
	  float dr2 = deltaR2(c.eta(), c.phi(), gen.eta(), gen.phi());
	  if (dr2 < dr2best) {
	    //dr2best = dr2; 
	    mc_.fill(gen);
	    mc_.dr = std::sqrt(dr2best);
	    matchedGen.push_back(index);
	  }
	}
	for (auto & v : reco_) v.fill(c);
	tree_->Fill();
    }
    int pCount=0;
    for (const reco::GenParticle &gen : *genparticles) {
      pCount++;
      if(!gen.isDirectPromptTauDecayProductFinalState()) continue;
      if(fabs(gen.pdgId()) > 10 && fabs(gen.pdgId()) < 17) continue;
      bool pMatch = false;
      for(unsigned int i0 = 0; i0 < matchedGen.size(); i0++) if(matchedGen[i0] == pCount) pMatch = true;
      if(pMatch) continue;
      mc_.fill(gen);
      mc_.dr=-99;
      for (auto & v : reco_) v.fill(dummy);
      tree_->Fill();
    }
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauNTuplizer);
