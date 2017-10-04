// -*- C++ -*-
//
// Package:    Giovanni/NTuplizer
// Class:      NTuplizer
// 
/**\class NTuplizer NTuplizer.cc Giovanni/NTuplizer/plugins/NTuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giovanni Petrucciani
//         Created:  Thu, 01 Sep 2016 11:30:38 GMT
//
//

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "FastPUPPI/NtupleProducer/interface/SimpleCalibrations.h"

#include <cstdint>
#include <TTree.h>
#include <TRandom3.h>

class JetMetNTuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit JetMetNTuplizer(const edm::ParameterSet&);
      ~JetMetNTuplizer();

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

      struct JetInput {
         std::string name;
         edm::EDGetTokenT<edm::View<reco::Jet>> token;
         l1tpf::SimpleCorrEm jec;
         JetInput() {}
         JetInput(const std::string &n, const edm::EDGetTokenT<edm::View<reco::Jet>> &t, const l1tpf::SimpleCorrEm & c) : name(n), token(t), jec(c) {}
      };
      std::vector<JetInput> jets_;
      std::vector<std::pair<std::string,StringCutObjectSelector<reco::Jet>>> jetSels_;

      struct MetInput {
         std::string name;
         edm::EDGetTokenT<edm::View<reco::MET>> token;
         MetInput() {}
         MetInput(const std::string &n, const edm::EDGetTokenT<edm::View<reco::MET>> &t) : name(n), token(t) {}
      };
      std::vector<MetInput> mets_;

      struct JetBranch {
        std::string name;
        int   count_raw, count_corr;
        float ht_raw, ht_corr;
        float mht_raw, mht_corr;
        JetBranch(const std::string &n) : name(n) {}
        void makeBranches(TTree *tree) {
            tree->Branch((name+"_ht_raw").c_str(),   &ht_raw,   (name+"_ht_raw/F").c_str());
            tree->Branch((name+"_mht_raw").c_str(),   &mht_raw,   (name+"_mht_raw/F").c_str());
            tree->Branch((name+"_num_raw").c_str(), &count_raw, (name+"_num_raw/I").c_str());
            tree->Branch((name+"_ht_corr").c_str(),   &ht_corr,   (name+"_ht_corr/F").c_str());
            tree->Branch((name+"_mht_corr").c_str(),   &mht_corr,   (name+"_mht_corr/F").c_str());
            tree->Branch((name+"_num_corr").c_str(), &count_corr, (name+"_num_corr/I").c_str());
        }
        void fill(const edm::View<reco::Jet> & jets, const  l1tpf::SimpleCorrEm & corr, const StringCutObjectSelector<reco::Jet> & sel) {  
            count_raw = 0; ht_raw = 0;
            count_corr = 0; ht_corr = 0;
            reco::Particle::LorentzVector sum_raw, sum_corr;
            for (reco::Jet jet : jets) {
                if (sel(jet)) { count_raw++; ht_raw += std::min(jet.pt(),500.); sum_raw += jet.p4(); }
                jet.setP4(reco::Particle::PolarLorentzVector(corr(jet.pt(), std::abs(jet.eta())), jet.eta(), jet.phi(), jet.mass()));
                if (sel(jet)) { count_corr++; ht_corr += std::min(jet.pt(),500.); sum_corr += jet.p4(); }
            }
        }
      };
      std::vector<JetBranch> jetbranches_;

      struct MetBranch {
       std::string name;
       float met, met_phi;
       MetBranch(const std::string &n) : name(n) {}
       void makeBranches(TTree *tree) {
         tree->Branch((name).c_str(),   &met,   (name+"/F").c_str());
         tree->Branch((name+"_phi").c_str(), &met_phi, (name+"_phi/F").c_str());
       }
       void fill(const edm::View<reco::MET> & recoMets) {  
         met = 0; met_phi = 0;
         for (reco::MET recoMet : recoMets) { met = recoMet.et(); met_phi = recoMet.phi(); }
       }
      };
      std::vector<MetBranch> metbranches_;      
       
      TTree *tree_;
      uint32_t run_, lumi_; uint64_t event_;      
};

JetMetNTuplizer::JetMetNTuplizer(const edm::ParameterSet& iConfig) 
{
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("tree","tree");
    tree_->Branch("run",  &run_,  "run/i");
    tree_->Branch("lumi", &lumi_, "lumi/i");
    tree_->Branch("event", &event_, "event/l");

    // jet branches
    edm::ParameterSet jets = iConfig.getParameter<edm::ParameterSet>("jets");
    edm::ParameterSet jecs = iConfig.getParameter<edm::ParameterSet>("jecs");
    for (const std::string & name : jets.getParameterNamesForType<edm::InputTag>()) {
        jets_.emplace_back(name, consumes<edm::View<reco::Jet>>(jets.getParameter<edm::InputTag>(name)), l1tpf::SimpleCorrEm(jecs, name));
    }
    edm::ParameterSet sels = iConfig.getParameter<edm::ParameterSet>("sels");
    for (const std::string & name : sels.getParameterNamesForType<std::string>()) {
        jetSels_.emplace_back(name, StringCutObjectSelector<reco::Jet>(sels.getParameter<std::string>(name)));
    }
    jetbranches_.reserve(jets_.size() * jetSels_.size());
    for (const auto & j : jets_) {
        for (const auto & s : jetSels_) {
            jetbranches_.emplace_back(j.name + s.first);
        }
    }
    for (auto & b : jetbranches_) b.makeBranches(tree_);

    //met branches
    edm::ParameterSet mets = iConfig.getParameter<edm::ParameterSet>("mets");
    for (const std::string & name : mets.getParameterNamesForType<edm::InputTag>()) {
        mets_.emplace_back(name, consumes<edm::View<reco::MET> >(mets.getParameter<edm::InputTag>(name)));
    }
    metbranches_.reserve(mets_.size());
    for (const auto & m : mets_) metbranches_.emplace_back(m.name);
    for (auto & b : metbranches_) b.makeBranches(tree_);    
}

JetMetNTuplizer::~JetMetNTuplizer() { }

// ------------ method called for each event  ------------
void
JetMetNTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    run_  = iEvent.id().run();
    lumi_ = iEvent.id().luminosityBlock();
    event_ = iEvent.id().event();
    
    int ibranch = 0;
    for (const auto & j : jets_) {
        edm::Handle<edm::View<reco::Jet>> jets;
        iEvent.getByToken(j.token, jets);
        for (const auto & s : jetSels_) {
            jetbranches_[ibranch++].fill(*jets, j.jec, s.second);
        } 
    }

    ibranch = 0;
    for (const auto & m : mets_) {
        edm::Handle<edm::View<reco::MET>> mets;
        iEvent.getByToken(m.token, mets);
        metbranches_[ibranch++].fill(*mets);
    }
    
    tree_->Fill();
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetMetNTuplizer);
