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

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include <cstdint>
#include <TTree.h>
#include <TRandom3.h>

class JetNTuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit JetNTuplizer(const edm::ParameterSet&);
      ~JetNTuplizer();

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

      std::vector<std::pair<std::string,edm::EDGetTokenT<edm::View<reco::Jet>>>> jets_;
      std::vector<std::pair<std::string,StringCutObjectSelector<reco::Jet>>> jetSels_;

      TTree *tree_;
      uint32_t run_, lumi_; uint64_t event_;

      struct JetBranch {
        std::string name;
        int   count;
        float ht;
        JetBranch(const std::string &n) : name(n), count(), ht() {}
        void makeBranches(TTree *tree) {
            tree->Branch((name+"_ht").c_str(),   &ht,   (name+"_ht/F").c_str());
            tree->Branch((name+"_num").c_str(), &count, (name+"_num/I").c_str());
        }
        void fill(const edm::View<reco::Jet> & jets, const StringCutObjectSelector<reco::Jet> & sel) {  
            count = 0; ht = 0;
            for (const reco::Jet & jet : jets) {
                if (sel(jet)) { count++; ht += std::min(jet.pt(),500.); }
            }
        }
      };
      std::vector<JetBranch> branches_;
};

JetNTuplizer::JetNTuplizer(const edm::ParameterSet& iConfig) 
{
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("tree","tree");
    tree_->Branch("run",  &run_,  "run/i");
    tree_->Branch("lumi", &lumi_, "lumi/i");
    tree_->Branch("event", &event_, "event/l");

    edm::ParameterSet jets = iConfig.getParameter<edm::ParameterSet>("jets");
    for (const std::string & name : jets.getParameterNamesForType<edm::InputTag>()) {
        jets_.emplace_back(name, consumes<edm::View<reco::Jet>>(jets.getParameter<edm::InputTag>(name)));
    }
    edm::ParameterSet sels = iConfig.getParameter<edm::ParameterSet>("sels");
    for (const std::string & name : sels.getParameterNamesForType<std::string>()) {
        jetSels_.emplace_back(name, StringCutObjectSelector<reco::Jet>(sels.getParameter<std::string>(name)));
    }
    branches_.reserve(jets_.size() * jetSels_.size());
    for (const auto & j : jets_) {
        for (const auto & s : jetSels_) {
            branches_.emplace_back(j.first + s.first);
        }
    }
    for (auto & b : branches_) b.makeBranches(tree_);
    
}

JetNTuplizer::~JetNTuplizer() { }

// ------------ method called for each event  ------------
void
JetNTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    run_  = iEvent.id().run();
    lumi_ = iEvent.id().luminosityBlock();
    event_ = iEvent.id().event();
    
    int ibranch = 0;
    for (const auto & j : jets_) {
        edm::Handle<edm::View<reco::Jet>> jets;
        iEvent.getByToken(j.second, jets);
        for (const auto & s : jetSels_) {
            branches_[ibranch++].fill(*jets, s.second);
        } 
    }
    tree_->Fill();
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetNTuplizer);
