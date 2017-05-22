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

#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <FastPUPPI/NtupleProducer/interface/L1TPFUtils.h>
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
                name_(name),prop_(false)
            {
                if      (name.find("Ecal") != std::string::npos) prop_ = true;
                else if (name.find("Hcal") != std::string::npos) prop_ = true;
                else if (name.find("Calo") != std::string::npos) prop_ = true;
                const std::vector<edm::InputTag> & tags = iConfig.getParameter< std::vector<edm::InputTag>>(name);
                for (const auto & tag : tags) tokens_.push_back(coll.consumes<reco::CandidateView>(tag));
            }
            const std::string & name() const { return name_; }
            bool  prop() const { return prop_; }
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
            bool prop_;
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
            int number(float dr, float threshold) const { 
                float dr2 = dr*dr, absthreshold = sum()*threshold;
                int mysum = 0;
                for (const auto & p : ptdr2) {
                    if (p.second < dr2 && p.first > absthreshold) mysum++;
                }
                return mysum;
            }
            float mindr(float threshold) const {
                float best = 9999, absthreshold = sum()*threshold;
                for (const auto & p : ptdr2) {
                    if (p.second < best && p.first > absthreshold) best = p.second;
                }
                return std::sqrt(best);
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
class ResponseNTuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns>  {
   public:
      explicit ResponseNTuplizer(const edm::ParameterSet&);
      ~ResponseNTuplizer();

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


      edm::EDGetTokenT<std::vector<reco::GenJet>> genjets_;
      edm::EDGetTokenT<std::vector<reco::GenParticle>> genparticles_;
      bool isParticleGun_;
      std::unique_ptr<TRandom> random_;
      bool doRandom_;
      TTree *tree_;
      uint32_t run_, lumi_; uint64_t event_;
      struct McVars {
         float pt, pt02, eta, phi, iso02, iso04;
         int charge; float caloeta, calophi;
         int   id;
         void makeBranches(TTree *tree) {
            tree->Branch("mc_pt", &pt, "mc_pt/F");
            tree->Branch("mc_pt02", &pt02, "mc_pt02/F");
            tree->Branch("mc_eta", &eta, "mc_eta/F");
            tree->Branch("mc_phi", &phi, "mc_phi/F");
            tree->Branch("mc_iso02", &iso02, "mc_iso02/F");
            tree->Branch("mc_iso04", &iso04, "mc_iso04/F");
            tree->Branch("mc_id", &id, "mc_id/I");
            tree->Branch("mc_q", &charge, "mc_q/I");
            tree->Branch("mc_caloeta", &caloeta, "mc_caloeta/F");
            tree->Branch("mc_calophi", &calophi, "mc_calophi/F");
         }
         void fillP4(const reco::Candidate &c) {
             pt = c.pt(); eta = c.eta(); phi = c.phi();
             caloeta = eta; calophi = phi; charge = 0;
         }
         void fillPropagated(const reco::Candidate &c, float bz) {
             if (c.charge() != 0) {
                charge = c.charge();
                math::XYZTLorentzVector momentum(c.px(),c.py(),c.pz(),c.energy());
                math::XYZTLorentzVector vertex(c.vx(),c.vy(),c.vz(),0.);
                std::vector<double> lVars;
                l1tpf::propagate(1,lVars,momentum,vertex,c.charge(),bz);
                caloeta = lVars[4]; calophi = lVars[5];
            }
         }

      } mc_;
      struct RecoVars {
         float pt, pt02, ptbest, pthighest, mindr025; int n025, n010;
         void makeBranches(const std::string &prefix, TTree *tree) {
             tree->Branch((prefix+"_pt").c_str(),   &pt,   (prefix+"_pt/F").c_str());
             tree->Branch((prefix+"_pt02").c_str(), &pt02, (prefix+"_pt02/F").c_str());
             tree->Branch((prefix+"_ptbest").c_str(), &ptbest, (prefix+"_ptbest/F").c_str());
             tree->Branch((prefix+"_pthighest").c_str(), &pthighest, (prefix+"_pthighest/F").c_str());
             tree->Branch((prefix+"_mindr025").c_str(), &mindr025, (prefix+"_mindr025/F").c_str());
             tree->Branch((prefix+"_n025").c_str(), &n025, (prefix+"_n025/I").c_str());
             tree->Branch((prefix+"_n010").c_str(), &n010, (prefix+"_n010/I").c_str());
         }
         void fill(const std::vector<::SimpleObject> & objects, float eta, float phi) {
             ::InCone incone(objects, eta, phi, 0.4);
             pt = incone.sum();
             pt02 = incone.sum(0.2);
             ptbest = incone.nearest();
             pthighest = incone.max();
             mindr025 =  incone.mindr(0.25);
             n025 = incone.number(0.2, 0.25);
             n010 = incone.number(0.2, 0.10);
         }
      };
      std::vector<std::pair<::MultiCollection,RecoVars>> reco_;
      template<typename T> struct 
      CopyT {
            std::string name;
            edm::EDGetTokenT<T> token;
            T buffer;
            CopyT(const edm::InputTag &tag, edm::EDGetTokenT<T> tok) :
                name(tag.label()+tag.instance()),
                token(tok),
                buffer() {}
            void get(const edm::Event &iEv) {
                edm::Handle<T> handle;
                iEv.getByToken(token, handle);
                buffer = *handle;
            }
            void makeBranches(TTree *tree) {
                tree->Branch((name).c_str(), &buffer, (name+"/"+typecode()).c_str());
            }
            static std::string typecode() { assert(false); }
      };
      std::vector<CopyT<unsigned>> copyUInts_;
      float bZ_;
};
template<> std::string ResponseNTuplizer::CopyT<unsigned>::typecode() { return "i"; }
template<> std::string ResponseNTuplizer::CopyT<float>::typecode() { return "F"; }

ResponseNTuplizer::ResponseNTuplizer(const edm::ParameterSet& iConfig) :
    genjets_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
    genparticles_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    isParticleGun_(iConfig.getParameter<bool>("isParticleGun")),
    random_(new TRandom3()), doRandom_(iConfig.getParameter<bool>("doRandom"))
{
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("tree","tree");
    tree_->Branch("run",  &run_,  "run/i");
    tree_->Branch("lumi", &lumi_, "lumi/i");
    tree_->Branch("event", &event_, "event/l");

    edm::ParameterSet objs = iConfig.getParameter<edm::ParameterSet>("objects");
    auto reconames = objs.getParameterNamesForType<std::vector<edm::InputTag>>();
    for (const std::string & name : reconames) {
        reco_.emplace_back(::MultiCollection(objs,name,consumesCollector()),RecoVars());
    }
    for (const edm::InputTag &tag : iConfig.getParameter<std::vector<edm::InputTag>>("copyUInts")) {
        copyUInts_.emplace_back(tag.instance(), consumes<unsigned>(tag));
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
    for (auto & c : copyUInts_) c.makeBranches(tree_);
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
    for (auto & c : copyUInts_) c.get(iEvent);
    mc_.id = 998; tree_->Fill(); // so that we write only one per event
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
                mc_.fillPropagated(*match, bZ_);
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
            recopair.second.fill(recopair.first.objects(),
                recopair.first.prop() ? mc_.caloeta : mc_.eta,
                recopair.first.prop() ? mc_.calophi : mc_.phi);
        }
        tree_->Fill();
    }
    // now let's throw a few "random" cones
    if (doRandom_) {
        for (int icone = 0; icone < 7; ++icone) {
            mc_.pt  = 10;
            mc_.eta = random_->Rndm() * 10 - 5;
            mc_.phi = M_PI*(random_->Rndm() * 2 - 1);
            mc_.id = 999;
            mc_.iso02 = 0;
            mc_.iso04 = 0;
            bool badcone = false;
            for (const reco::GenParticle &gen : *genparticles) {
                if ((gen.isPromptFinalState() || gen.isDirectPromptTauDecayProductFinalState()) && gen.pt() > 0.5) {
                    if (::deltaR2(mc_.eta,mc_.phi,gen.eta(),gen.phi()) < 0.36f) {
                        badcone = true; break;
                    }
                }
            }
            if (badcone) continue;
            for (const reco::GenJet & j : *genjets) {
                if (j.pt() > 5 && ::deltaR2(mc_.eta,mc_.phi,j.eta(),j.phi()) < 0.64f) {
                    badcone = true; break;
                }
            }
            if (badcone) continue;
            for (auto & recopair : reco_) {
                recopair.second.fill(recopair.first.objects(), mc_.eta, mc_.phi);
            }
            tree_->Fill();

        }
    }
    for (auto & recopair : reco_) {
        recopair.first.clear();
    }

}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ResponseNTuplizer);
