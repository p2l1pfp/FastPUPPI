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

class MetNTuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MetNTuplizer(const edm::ParameterSet&);
      ~MetNTuplizer();

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;    

      struct MetInput {
         std::string name;
         edm::EDGetTokenT<edm::View<reco::MET>> token;
         MetInput() {}
         MetInput(const std::string &n, const edm::EDGetTokenT<edm::View<reco::MET>> &t) : name(n), token(t) {}
      };
      std::vector<MetInput> mets_;
      
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
      std::vector<MetBranch> branches_; 
        
      TTree *tree_;
      uint32_t run_, lumi_; uint64_t event_;      

};

MetNTuplizer::MetNTuplizer(const edm::ParameterSet& iConfig) 
{
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("tree","tree");
    tree_->Branch("run",  &run_,  "run/i");
    tree_->Branch("lumi", &lumi_, "lumi/i");
    tree_->Branch("event", &event_, "event/l");
    
    edm::ParameterSet mets = iConfig.getParameter<edm::ParameterSet>("mets");
    for (const std::string & name : mets.getParameterNamesForType<edm::InputTag>()) {
        mets_.emplace_back(name, consumes<edm::View<reco::MET> >(mets.getParameter<edm::InputTag>(name)));
    }
    branches_.reserve(mets_.size());
    for (const auto & m : mets_) branches_.emplace_back(m.name);
    for (auto & b : branches_) b.makeBranches(tree_);
        
}

MetNTuplizer::~MetNTuplizer() { }

// ------------ method called for each event  ------------
void
MetNTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    run_  = iEvent.id().run();
    lumi_ = iEvent.id().luminosityBlock();
    event_ = iEvent.id().event();

    int ibranch = 0;
    for (const auto & m : mets_) {
        edm::Handle<edm::View<reco::MET>> mets;
        iEvent.getByToken(m.token, mets);
        branches_[ibranch++].fill(*mets);
    }
    
    tree_->Fill();    
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MetNTuplizer);    
