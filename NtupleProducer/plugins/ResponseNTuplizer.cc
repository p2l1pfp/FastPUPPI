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

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
//#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include <cstdint>
#include <TTree.h>
#include <TRandom3.h>

namespace {
    struct SimpleObject {
        float pt, eta, phi;
        SimpleObject(float apt, float aneta, float aphi) : pt(apt), eta(aneta), phi(aphi) {}
        bool operator<(const SimpleObject &other) const { return eta < other.eta; }
        bool operator<(const float &other) const { return eta < other; }
    };
    class MultiCollection {
        public:
            MultiCollection(const edm::ParameterSet &iConfig, const std::string &name, edm::ConsumesCollector && coll) :
                name_(name)
            {
                const std::vector<edm::InputTag> & tags = iConfig.getParameter< std::vector<edm::InputTag>>(name);
                for (const auto & tag : tags) tokens_.push_back(coll.consumes<reco::CandidateView>(tag));
            }
            const std::string & name() const { return name_; }
            void get(const edm::Event &iEvent) {
                edm::Handle<reco::CandidateView> handle;
                for (const auto & token : tokens_) {
                    iEvent.getByToken(token, handle);
                    for (const reco::Candidate & c : *handle) {
                        objects_.emplace_back(c.pt(), c.eta(), c.phi());
                    }
                }
                std::sort(objects_.begin(), objects_.end());
            }        
            const std::vector<SimpleObject> & objects() const { return objects_; }
            void clear() { objects_.clear(); }
        private:
            std::string name_;
            std::vector<edm::EDGetTokenT<reco::CandidateView>> tokens_;
            std::vector<SimpleObject> objects_;
    };
    class InCone {
        public:
            InCone(const std::vector<SimpleObject> & objects, float eta, float phi, float dr) {
                auto first = std::lower_bound(objects.begin(), objects.end(), eta-dr-0.01f); // small offset to avoid dealing with ==
                auto end   = std::lower_bound(objects.begin(), objects.end(), eta+dr+0.01f);
                float dr2 = dr*dr;
                for (auto it = first; it < end; ++it) {
                    float mydr2 = ::deltaR2(eta,phi, it->eta,it->phi);
                    if (mydr2 < dr2) ptdr2.emplace_back(it->pt, mydr2);
                }
            }
            float sum() const { 
                float mysum = 0;
                for (const auto & p : ptdr2) mysum += p.first;
                return mysum;
            }
            float sum(float dr) const { 
                float dr2 = dr*dr;
                float mysum = 0;
                for (const auto & p : ptdr2) {
                    if (p.second < dr2) mysum += p.first;
                }
                return mysum;
            }
            float nearest() const {
                std::pair<float,float> best(0,9999);
                for (const auto & p : ptdr2) {
                    if (p.second < best.second) best = p;
                }
                return best.first;
            }
            float max() const {
                float best = 0;
                for (const auto & p : ptdr2) {
                    if (p.first > best) best = p.first;
                }
                return best;
            }

        private:
            std::vector<std::pair<float,float>> ptdr2;
    };

}
class ResponseNTuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ResponseNTuplizer(const edm::ParameterSet&);
      ~ResponseNTuplizer();

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

      edm::EDGetTokenT<std::vector<reco::GenJet>> genjets_;
      edm::EDGetTokenT<std::vector<reco::GenParticle>> genparticles_;
      bool isParticleGun_;
      TRandom3 * random_;
      TTree *tree_;
      uint32_t run_, lumi_; uint64_t event_;
      struct McVars {
         float pt, pt02, eta, phi, iso02, iso04;
         int   id;
         void makeBranches(TTree *tree) {
            tree->Branch("mc_pt", &pt, "mc_pt/F");
            tree->Branch("mc_pt02", &pt02, "mc_pt02/F");
            tree->Branch("mc_eta", &eta, "mc_eta/F");
            tree->Branch("mc_phi", &phi, "mc_phi/F");
            tree->Branch("mc_iso02", &iso02, "mc_iso02/F");
            tree->Branch("mc_iso04", &iso04, "mc_iso04/F");
            tree->Branch("mc_id", &id, "mc_id/I");
         }
         void fillP4(const reco::Candidate &c) {
             pt = c.pt(); eta = c.eta(); phi = c.phi();
         }
      } mc_;
      struct RecoVars {
         float pt, pt02, ptbest, pthighest;
         void makeBranches(const std::string &prefix, TTree *tree) {
             tree->Branch((prefix+"_pt").c_str(),   &pt,   (prefix+"_pt/F").c_str());
             tree->Branch((prefix+"_pt02").c_str(), &pt02, (prefix+"_pt02/F").c_str());
             tree->Branch((prefix+"_ptbest").c_str(), &ptbest, (prefix+"_ptbest/F").c_str());
             tree->Branch((prefix+"_pthighest").c_str(), &pthighest, (prefix+"_pthighest/F").c_str());
         }
         void fill(const std::vector<::SimpleObject> & objects, float eta, float phi) {
             ::InCone incone(objects, eta, phi, 0.4);
             pt = incone.sum();
             pt02 = incone.sum(0.2);
             ptbest = incone.nearest();
             pthighest = incone.max();
         }
      };
      std::vector<std::pair<::MultiCollection,RecoVars>> reco_;
};

ResponseNTuplizer::ResponseNTuplizer(const edm::ParameterSet& iConfig) :
    genjets_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
    genparticles_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    isParticleGun_(iConfig.getParameter<bool>("isParticleGun")),
    random_(new TRandom3())
{
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("tree","tree");
    tree_->Branch("run",  &run_,  "run/i");
    tree_->Branch("lumi", &lumi_, "lumi/i");
    tree_->Branch("event", &event_, "event/l");

    auto reconames = iConfig.getParameterNamesForType<std::vector<edm::InputTag>>();
    for (const std::string & name : reconames) {
        reco_.emplace_back(::MultiCollection(iConfig,name,consumesCollector()),RecoVars());
    }
}

ResponseNTuplizer::~ResponseNTuplizer() { }

// ------------ method called once each job just before starting event loop  ------------
void 
ResponseNTuplizer::beginJob()
{
    mc_.makeBranches(tree_);
    for (auto & pair : reco_) {
        pair.second.makeBranches(pair.first.name(), tree_);
    }
}


// ------------ method called for each event  ------------
void
ResponseNTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    run_  = iEvent.id().run();
    lumi_ = iEvent.id().luminosityBlock();
    event_ = iEvent.id().event();

    edm::Handle<std::vector<reco::GenJet>> genjets;
    edm::Handle<std::vector<reco::GenParticle>> genparticles;
    iEvent.getByToken(genjets_, genjets);
    iEvent.getByToken(genparticles_, genparticles);

    std::vector<const reco::GenParticle *> prompts, taus;
    for (const reco::GenParticle &gen : *genparticles) {
        if (isParticleGun_) {
            if (gen.statusFlags().isPrompt() == 1) prompts.push_back(&gen);
            continue;
        }
        if ((gen.isPromptFinalState() || gen.isDirectPromptTauDecayProductFinalState()) && (std::abs(gen.pdgId()) == 11 || std::abs(gen.pdgId()) == 13) && gen.pt() > 5) {
            prompts.push_back(&gen);
        } else if (gen.isPromptFinalState() && std::abs(gen.pdgId()) == 22 && gen.pt() > 10) {
            prompts.push_back(&gen);
        } else if (abs(gen.pdgId()) == 15 && gen.isPromptDecayed()) {
            taus.push_back(&gen);
        }
    }

    for (auto & recopair : reco_) {
        recopair.first.get(iEvent);
    }
    for (const reco::GenJet & j : *genjets) {
        bool ok = true;
        const reco::Candidate * match = nullptr;
        for (const reco::GenParticle * p : prompts) {
            if (::deltaR2(*p, j) < 0.16f) {
                if (match != nullptr) { ok = false; break; }
                else { match = p; }
            }
        }
        if (!ok) continue;
        if (!match) {
            // look for a tau
            for (const reco::GenParticle * p : taus) {
                if (::deltaR2(*p, j) < 0.16f) {
                    if (match != nullptr) { ok = false; break; }
                    else { match = p; }
                }
            }
            if (!ok) continue;
            if (match != nullptr && match->numberOfDaughters() == 2 && std::abs(match->daughter(0)->pdgId()) + std::abs(match->daughter(1)->pdgId()) == 211+16) {
                // one-prong tau, consider it a pion
                match = (std::abs(match->daughter(0)->pdgId()) == 211 ? match->daughter(0) : match->daughter(1));
            }
        }
        if (match != nullptr) {
            if (std::abs(match->pdgId()) == 15) {
                reco::Particle::LorentzVector pvis;
                for (unsigned int i = 0, n = match->numberOfDaughters(); i < n; ++i) {
                    const reco::Candidate *dau = match->daughter(i);
                    if (std::abs(dau->pdgId()) == 12 || std::abs(dau->pdgId()) == 14 || std::abs(dau->pdgId()) == 16) {
                        continue;
                    } 
                    pvis += dau->p4();
                }
                mc_.pt  = pvis.Pt();
                mc_.eta = pvis.Eta();
                mc_.phi = pvis.Phi();
            } else {
                mc_.fillP4(*match);
            }
            mc_.id = std::abs(match->pdgId());
            mc_.iso04 = j.pt()/mc_.pt - 1;
            mc_.iso02 = 0;
            for (const auto &dptr : j.daughterPtrVector()) {
                if (::deltaR2(*dptr, *match) < 0.04f) {
                    mc_.iso02 += dptr->pt();
                }
            }
            mc_.iso02 = mc_.iso02/mc_.pt - 1;
        } else {
            if (j.pt() < 20) continue;
            mc_.fillP4(j);
            mc_.id = 0;
            mc_.iso02 = 0;
            mc_.iso04 = 0;
        }
        for (auto & recopair : reco_) {
            recopair.second.fill(recopair.first.objects(), mc_.eta, mc_.phi);
        }
        tree_->Fill();
    }
    for (auto & recopair : reco_) {
        recopair.first.clear();
    }

}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ResponseNTuplizer);
