// -*- C++ -*-
//
// Package:    FastPUPPI/L1PFCompare
// Class:      L1PFCompare
// 
/**\class L1PFCompare L1PFCompare.cc FastPUPPI/NtuplerProducer/plugins/L1PFCompare.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Dylan Sheldon Rankin
//         Created:  Fri, 22 Jun 2018 19:18:53 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include  "DataFormats/L1Trigger/interface/Vertex.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/BJetId.h"


#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class L1PFCompare : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit L1PFCompare(const edm::ParameterSet&);
      ~L1PFCompare();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      static bool sortff(std::pair<float,float>, std::pair<float,float>);
      static bool sortPt(std::pair<TLorentzVector,int>, std::pair<TLorentzVector,int>);
      static bool sortTupPt(std::tuple<TLorentzVector,int,int>, std::tuple<TLorentzVector,int,int>);
      static bool sortTupfffPt(std::tuple<TLorentzVector,int,float,float,float>, std::tuple<TLorentzVector,int,float,float,float>);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void searchDaughters(reco::GenParticle, bool&, int);

      // ----------member data ---------------------------

      TTree* objects;

      std::vector<std::pair<TLorentzVector,int>> emcalo_;
      std::vector<std::pair<TLorentzVector,int>> egcalo_;
      std::vector<std::pair<TLorentzVector,int>> calo_;
      std::vector<std::pair<TLorentzVector,int>> pf_;
      std::vector<std::pair<TLorentzVector,int>> pup_;
      std::vector<std::pair<TLorentzVector,int>> gen_;
      std::vector<std::pair<TLorentzVector,int>> l1jet_;
      std::vector<std::pair<TLorentzVector,int>> recojet_;

      std::vector<float> vx_pf_;
      std::vector<float> vy_pf_;
      std::vector<float> vz_pf_;
      std::vector<float> vx_pup_;
      std::vector<float> vy_pup_;
      std::vector<float> vz_pup_;

      std::vector<float> vz_;

      //edm::EDGetTokenT<std::vector<l1t::PFCandidate>> emcaloToken_;
      edm::EDGetTokenT<l1t::HGCalMulticlusterBxCollection> emcaloToken_;
      edm::EDGetTokenT<BXVector<l1t::EGamma>> egcaloToken_;
      edm::EDGetTokenT<std::vector<l1t::PFCandidate>> caloToken_;
      edm::EDGetTokenT<std::vector<l1t::PFCandidate>> pfToken_;
      edm::EDGetTokenT<std::vector<l1t::PFCandidate>> pupToken_;
      edm::EDGetTokenT<std::vector<l1t::Vertex>> vtxToken_;
      edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
      edm::EDGetTokenT<edm::View<l1t::PFJet>> l1jetToken_;
      edm::EDGetTokenT<edm::ValueMap<float>> l1bName_;
      edm::EDGetTokenT<pat::JetCollection> recojetToken_;
      std::string recobName_;
      double minPt_;
      double maxEta_;
      unsigned int maxN_;

      std::vector<int> genIDs_;
      std::vector<int> addGenIDs_;
      std::vector<int> genStatuses_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1PFCompare::L1PFCompare(const edm::ParameterSet& iConfig):
	emcaloToken_(consumes<l1t::HGCalMulticlusterBxCollection>(iConfig.getParameter<edm::InputTag>("emcalo"))),
	egcaloToken_(consumes<BXVector<l1t::EGamma>>(iConfig.getParameter<edm::InputTag>("egcalo"))),
	caloToken_(consumes<std::vector<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("calo"))),
	pfToken_(consumes<std::vector<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pf"))),
	pupToken_(consumes<std::vector<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pup"))),
	vtxToken_(consumes<std::vector<l1t::Vertex>>(iConfig.getParameter<edm::InputTag>("vtx"))), // "L1PrimaryVertexTag"?
	genToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("generator"))),
	l1jetToken_(consumes<edm::View<l1t::PFJet>>(iConfig.getParameter<edm::InputTag>("l1jet"))),
        l1bName_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("l1bName"))),
	recojetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("recojet"))),
        recobName_(iConfig.getParameter<std::string>("recobName")),
        minPt_(iConfig.getParameter<double>("minPt")),
        maxEta_(iConfig.getParameter<double>("maxEta")),
        maxN_(iConfig.getParameter<unsigned int>("maxN")),
        genIDs_(iConfig.getParameter<std::vector<int>>("genIDs")),
        addGenIDs_(iConfig.getParameter<std::vector<int>>("addGenIDs")),
        genStatuses_(iConfig.getParameter<std::vector<int>>("genStatuses"))
{
   //now do what ever initialization is needed
   usesResource("TFileService");

   edm::Service<TFileService> fs;
   objects = fs->make<TTree>("objects","objects");

   objects->Branch("emcalo", "std::vector<std::pair<TLorentzVector,int>>", &emcalo_, 32000, 0);
   objects->Branch("egcalo", "std::vector<std::pair<TLorentzVector,int>>", &egcalo_, 32000, 0);
   objects->Branch("calo", "std::vector<std::pair<TLorentzVector,int>>", &calo_, 32000, 0);
   objects->Branch("pf", "std::vector<std::pair<TLorentzVector,int>>", &pf_, 32000, 0);
   objects->Branch("pup", "std::vector<std::pair<TLorentzVector,int>>", &pup_, 32000, 0);
   objects->Branch("gen", "std::vector<std::pair<TLorentzVector,int>>", &gen_, 32000, 0);
   objects->Branch("l1jet", "std::vector<std::pair<TLorentzVector,int>>", &l1jet_, 32000, 0);
   objects->Branch("recojet", "std::vector<std::pair<TLorentzVector,int>>", &recojet_, 32000, 0);

   objects->Branch("pf_vx", &vx_pf_);
   objects->Branch("pf_vy", &vy_pf_);
   objects->Branch("pf_vz", &vz_pf_);

   objects->Branch("pup_vx", &vx_pup_);
   objects->Branch("pup_vy", &vy_pup_);
   objects->Branch("pup_vz", &vz_pup_);

   objects->Branch("vz", &vz_);

}




L1PFCompare::~L1PFCompare()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
void L1PFCompare::searchDaughters(reco::GenParticle p, bool& flagged, int iter) {
    if (not flagged and iter<=2) {
        iter++;
        const reco::GenParticleRefVector& daughterRefs = p.daughterRefVector();
        for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) {
            bool goodID = false;
            for (unsigned int d1 = 0; d1 < addGenIDs_.size(); d1++) {
                if (abs((*idr)->pdgId())==addGenIDs_[d1]) goodID = true;
            }
            if (goodID) {flagged = true; break;}
            else searchDaughters((*(*idr)),flagged,iter);
        }
    }
    return;
}

// ------------ method called for each event  ------------
void
L1PFCompare::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   edm::Handle<l1t::HGCalMulticlusterBxCollection> emcaloHGClusters;
   bool doEmCalo = true;

   if (iEvent.getByToken(emcaloToken_, emcaloHGClusters)) {}
   else doEmCalo = false;
 
   edm::Handle<BXVector<l1t::EGamma>> egClusters;
   bool doEg = true;

   if (iEvent.getByToken(egcaloToken_, egClusters)) {}
   else doEg = false;
 
   edm::Handle<std::vector<l1t::PFCandidate>> caloParticles;
   bool doCalo = true;

   if (iEvent.getByToken(caloToken_, caloParticles)) {}
   else doCalo = false;

   edm::Handle<std::vector<l1t::PFCandidate>> pfParticles;
   iEvent.getByToken(pfToken_, pfParticles);

   edm::Handle<std::vector<l1t::PFCandidate>> pupParticles;
   iEvent.getByToken(pupToken_, pupParticles);

   edm::Handle<std::vector<l1t::Vertex>> vtxHandle;
   iEvent.getByToken(vtxToken_, vtxHandle);
   
   edm::Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByToken(genToken_, genParticles);
   
   edm::Handle<edm::View<l1t::PFJet>> l1jets;
   edm::Handle<edm::ValueMap<float> > l1ValueMap_;
   bool doL1Jets = true;

   if(iEvent.getByToken(l1jetToken_, l1jets)) {}
   else doL1Jets = false;
   if(iEvent.getByToken(l1bName_, l1ValueMap_)) {}
   else doL1Jets = false;

   edm::Handle<pat::JetCollection> recojets;
   bool doRecoJets = true;

   if (iEvent.getByToken(recojetToken_, recojets)) {}
   else doRecoJets = false;

   TLorentzVector tmp;
   std::pair<TLorentzVector, int> dummy;
   std::pair<float, float> dummyff;
   std::tuple<TLorentzVector, int, int> dumtup;
   std::tuple<TLorentzVector, int, float, float, float> dumtupfff;
   
   emcalo_.clear();
   egcalo_.clear();
   calo_.clear();
   pf_.clear();
   pup_.clear();
   gen_.clear();
   l1jet_.clear();
   recojet_.clear();

   vx_pf_.clear();
   vy_pf_.clear();
   vz_pf_.clear();
   vx_pup_.clear();
   vy_pup_.clear();
   vz_pup_.clear();

   vz_.clear();

   std::vector<std::pair<TLorentzVector,int>> emcaloTmp;
   std::vector<std::pair<TLorentzVector,int>> egcaloTmp;
   std::vector<std::pair<TLorentzVector,int>> caloTmp;
   std::vector<std::tuple<TLorentzVector,int,float,float,float>> pTmpTup;
   std::vector<std::pair<TLorentzVector,int>> pfTmp;
   std::vector<std::pair<TLorentzVector,int>> pupTmp;
   std::vector<std::pair<TLorentzVector,int>> l1jetTmp;
   std::vector<std::pair<TLorentzVector,int>> recojetTmp;
   std::vector<std::pair<float,float>> z0pTtmp;
   
   for (const l1t::Vertex/*TkPrimaryVertex*/ & vtx : *vtxHandle) {
       dummyff.first = vtx.pt(); //vtx.sum();
       dummyff.second = vtx.z0();//vtx.zvertex();
       z0pTtmp.push_back(dummyff);
   }
   std::sort(z0pTtmp.begin(),z0pTtmp.end(),L1PFCompare::sortff);
   for (size_t i = 0; i < std::min(z0pTtmp.size(),(size_t)(5)); i++) {vz_.push_back(z0pTtmp[i].second);}

   if (doEmCalo) {
       for (auto it = emcaloHGClusters->begin(0), ed = emcaloHGClusters->end(0); it != ed; ++it){
           if (it->pt()<minPt_ or fabs(it->eta())>maxEta_) continue;
           tmp.SetPtEtaPhiE(it->pt(),it->eta(),it->phi(),it->energy());
           dummy.first = tmp;
           dummy.second = int(it->hOverE()*1000.);
           emcaloTmp.push_back(dummy);
       }
       std::sort(emcaloTmp.begin(),emcaloTmp.end(),L1PFCompare::sortPt);
       for (size_t i = 0; i < std::min(emcaloTmp.size(),(size_t)maxN_); i++) {emcalo_.push_back(emcaloTmp[i]);}
   }
   if (doEg) {
       for (auto it = egClusters->begin(0), ed = egClusters->end(0); it != ed; ++it){
           if (it->pt()<minPt_ or fabs(it->eta())>maxEta_) continue;
           tmp.SetPtEtaPhiE(it->pt(),it->eta(),it->phi(),it->energy());
           dummy.first = tmp;
           dummy.second = it->hwQual();
           egcaloTmp.push_back(dummy);
       }
       std::sort(egcaloTmp.begin(),egcaloTmp.end(),L1PFCompare::sortPt);
       for (size_t i = 0; i < std::min(egcaloTmp.size(),(size_t)maxN_); i++) {egcalo_.push_back(egcaloTmp[i]);}
   }
   if (doCalo) {
       for (size_t i = 0; i < caloParticles->size(); i++){
           if ((*caloParticles)[i].pt()<minPt_ or fabs((*caloParticles)[i].eta())>maxEta_) continue;
           tmp.SetPtEtaPhiE((*caloParticles)[i].pt(),(*caloParticles)[i].eta(),(*caloParticles)[i].phi(),(*caloParticles)[i].energy());
           dummy.first = tmp;
           dummy.second = (*caloParticles)[i].pdgId();
           caloTmp.push_back(dummy);
       }
       std::sort(caloTmp.begin(),caloTmp.end(),L1PFCompare::sortPt);
       for (size_t i = 0; i < std::min(caloTmp.size(),(size_t)maxN_); i++) {calo_.push_back(caloTmp[i]);}
   }
   for (size_t i = 0; i < pfParticles->size(); i++){
       if ((*pfParticles)[i].pt()<minPt_ or fabs((*pfParticles)[i].eta())>maxEta_) continue;
       tmp.SetPtEtaPhiE((*pfParticles)[i].pt(),(*pfParticles)[i].eta(),(*pfParticles)[i].phi(),(*pfParticles)[i].energy());
       std::get<0>(dumtupfff) = tmp;
       std::get<1>(dumtupfff) = (*pfParticles)[i].pdgId();
       float x0_ = 0.;
       float y0_ = 0.;
       float z0_ = 0.;
       if ((*pfParticles)[i].pfTrack().isNonnull()) {
           x0_ = (*pfParticles)[i].pfTrack()->vx();
           y0_ = (*pfParticles)[i].pfTrack()->vy();
           z0_ = (*pfParticles)[i].pfTrack()->vz();
       }
       std::get<2>(dumtupfff) = x0_;
       std::get<3>(dumtupfff) = y0_;
       std::get<4>(dumtupfff) = z0_;
       pTmpTup.push_back(dumtupfff);
   }
   std::sort(pTmpTup.begin(),pTmpTup.end(),L1PFCompare::sortTupfffPt);
   for (size_t i = 0; i < std::min(pTmpTup.size(),(size_t)maxN_); i++) {
       dummy.first = std::get<0>(pTmpTup[i]);
       dummy.second = std::get<1>(pTmpTup[i]);
       pf_.push_back(dummy);
       vx_pf_.push_back(std::get<2>(pTmpTup[i]));
       vy_pf_.push_back(std::get<3>(pTmpTup[i]));
       vz_pf_.push_back(std::get<4>(pTmpTup[i]));
   }
   pTmpTup.clear();
   for (size_t i = 0; i < pupParticles->size(); i++){
       if ((*pupParticles)[i].pt()<minPt_ or fabs((*pupParticles)[i].eta())>maxEta_) continue;
       tmp.SetPtEtaPhiE((*pupParticles)[i].pt(),(*pupParticles)[i].eta(),(*pupParticles)[i].phi(),(*pupParticles)[i].energy());
       std::get<0>(dumtupfff) = tmp;
       std::get<1>(dumtupfff) = (*pupParticles)[i].pdgId();
       float x0_ = 0.;
       float y0_ = 0.;
       float z0_ = 0.;
       if ((*pupParticles)[i].pfTrack().isNonnull()) {
           x0_ = (*pupParticles)[i].pfTrack()->vx();
           y0_ = (*pupParticles)[i].pfTrack()->vy();
           z0_ = (*pupParticles)[i].pfTrack()->vz();
       }
       std::get<2>(dumtupfff) = x0_;
       std::get<3>(dumtupfff) = y0_;
       std::get<4>(dumtupfff) = z0_;
       pTmpTup.push_back(dumtupfff);
   }
   std::sort(pTmpTup.begin(),pTmpTup.end(),L1PFCompare::sortTupfffPt);
   for (size_t i = 0; i < std::min(pTmpTup.size(),(size_t)maxN_); i++) {
       dummy.first = std::get<0>(pTmpTup[i]);
       dummy.second = std::get<1>(pTmpTup[i]);
       pup_.push_back(dummy);
       vx_pup_.push_back(std::get<2>(pTmpTup[i]));
       vy_pup_.push_back(std::get<3>(pTmpTup[i]));
       vz_pup_.push_back(std::get<4>(pTmpTup[i]));
   }

   std::vector<std::tuple<TLorentzVector,int,int>> genFull;
   std::vector<edm::Ptr<reco::GenParticle>> lTauMothers;
   for (size_t i = 0; i < genParticles->size(); i++){
       const reco::GenParticle & p = (*genParticles).at(i);
       if (p.pt()<minPt_ or fabs(p.eta())>maxEta_) continue;
       bool badID = true;
       for (unsigned int d = 0; d < genIDs_.size(); d++) {
           if (p.pdgId()==genIDs_[d] and p.status()==genStatuses_[d]) {
               badID = false;
               break;
           }
       }
       if (badID) continue;
       bool flagged = false;
       int iter = 0;
       searchDaughters(p, flagged, iter);
       if (flagged) continue;
       int tauIndex = -1;
       if( (p.vx()*p.vx() + p.vy()*p.vy()) > 1. or p.vz() > 30.) continue;
       tmp.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy());
       bool already = false;
       for (unsigned int j = 0; j < genFull.size(); j++) {if (tmp.DeltaR(std::get<0>(genFull[j]))<0.3 && fabs((tmp.Pt()/std::get<0>(genFull[j]).Pt())-1.)<0.1 && std::get<1>(genFull[j])==p.pdgId()) {already = true; break;}}
       if (already) continue;
       std::get<0>(dumtup) = tmp;
       std::get<1>(dumtup) = p.pdgId();
       std::get<2>(dumtup) = tauIndex;
       genFull.push_back(dumtup);
   }
   for (size_t i = 0; i < genParticles->size(); i++){
       if (genFull.size()>=maxN_) break;
       const reco::GenParticle & p = (*genParticles).at(i);
       if (p.pt()<minPt_ or fabs(p.eta())>maxEta_) continue;
       if ( !(abs(p.pdgId())==5) ) continue;
       if( !p.isHardProcess()) continue;
       int tauIndex = -1;
       if( (p.vx()*p.vx() + p.vy()*p.vy()) > 1. or p.vz() > 30.) continue;
       tmp.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy());
       bool already = false;
       for (unsigned int j = 0; j < genFull.size(); j++) {if (tmp.DeltaR(std::get<0>(genFull[j]))<0.05 && fabs((tmp.Pt()/std::get<0>(genFull[j]).Pt())-1.)<0.1 && std::get<1>(genFull[j])==p.pdgId()) {already = true; break;}}
       if (already) continue;
       std::get<0>(dumtup) = tmp;
       std::get<1>(dumtup) = p.pdgId();
       std::get<2>(dumtup) = tauIndex;
       genFull.push_back(dumtup);
   }
   std::sort(genFull.begin(),genFull.end(),L1PFCompare::sortTupPt);
   for (size_t i = 0; i < std::min(genFull.size(),(size_t)maxN_); i++) {
       dummy.first = std::get<0>(genFull[i]);
       dummy.second = std::get<1>(genFull[i]);
       gen_.push_back(dummy);
   }

   if (doL1Jets) {
       size_t il1 = 0;
       for (const auto& srcjet : *l1jets) {
           if (srcjet.pt()<20. || fabs(srcjet.eta())>maxEta_) {
               il1++;
               continue;
           }
           tmp.SetPtEtaPhiE(srcjet.pt(), srcjet.eta(), srcjet.phi(), srcjet.energy());
           dummy.first = tmp;
           dummy.second = int(((*l1ValueMap_)[l1jets->refAt(il1)])*1000.);
           l1jetTmp.push_back(dummy);
           il1++;
       }
       std::sort(l1jetTmp.begin(),l1jetTmp.end(),L1PFCompare::sortPt);
       for (size_t i = 0; i < std::min(l1jetTmp.size(),(size_t)maxN_); i++) {l1jet_.push_back(l1jetTmp[i]);}
   }
//`pat::Jet::bDiscriminator(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&) const'
   if (doRecoJets) {
       for (const pat::Jet &j : *recojets) {
           if (j.pt() < 20. || fabs(j.eta())>maxEta_) {
               continue;
           }
           tmp.SetPtEtaPhiM(j.pt(), j.eta(), j.phi(), j.mass());
           dummy.first = tmp;
           dummy.second = int(std::max(0.f,j.bDiscriminator(recobName_))*1000.);
           recojetTmp.push_back(dummy);
       }
       std::sort(recojetTmp.begin(),recojetTmp.end(),L1PFCompare::sortPt);
       for (size_t i = 0; i < std::min(recojetTmp.size(),(size_t)maxN_); i++) {recojet_.push_back(recojetTmp[i]);}
   }

   objects->Fill();

}

bool L1PFCompare::sortff(std::pair<float,float> i, std::pair<float,float> j) {return i.first > j.first;}
bool L1PFCompare::sortPt(std::pair<TLorentzVector,int> i, std::pair<TLorentzVector,int> j) {return i.first.Pt() > j.first.Pt();}
bool L1PFCompare::sortTupPt(std::tuple<TLorentzVector,int,int> i, std::tuple<TLorentzVector,int,int> j) {return std::get<0>(i).Pt() > std::get<0>(j).Pt();}
bool L1PFCompare::sortTupfffPt(std::tuple<TLorentzVector,int,float,float,float> i, std::tuple<TLorentzVector,int,float,float,float> j) {return std::get<0>(i).Pt() > std::get<0>(j).Pt();}

// ------------ method called once each job just before starting event loop  ------------
void 
L1PFCompare::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1PFCompare::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1PFCompare::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1PFCompare);
