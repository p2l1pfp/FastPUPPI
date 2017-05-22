// -*- C++ -*-
//
// Package:    CaloNtupleProducer
// Class:      CaloNtupleProducer
// 
/**\class CaloNtupleProducer CaloNtupleProducer.cc Ntuplizer/CaloNtupleProducer/plugins/CaloNtupleProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Stanislava Sevova
//         Created:  Fri, 24 Jun 2016 14:57:31 GMT
// $Id$
//
//


// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FastPUPPI/NtupleProducer/interface/corrector.hh"
#include "FastPUPPI/NtupleProducer/interface/CaloClusterer.h"
#include "FastPUPPI/NtupleProducer/interface/SimpleCalibrations.h"
#include "FastPUPPI/NtupleProducer/interface/L1TPFUtils.h"

#include "DataFormats/Math/interface/deltaR.h"

// ROOT classes
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObject.h>

typedef std::vector<reco::PFCandidate> PFOutputCollection;
typedef std::vector<l1tpf::Particle>   L1PFCollection;

//--------------------------------------------------------------------------------------------------
class CaloNtupleProducer : public edm::EDProducer {
public:
  explicit CaloNtupleProducer(const edm::ParameterSet&);
  ~CaloNtupleProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  const bool zeroSuppress_, ecalOnly_;
  const std::vector<edm::InputTag> EcalTPTags_;
  const std::vector<edm::InputTag> HcalTPTags_;
  const edm::InputTag GenParTag_;
  const std::string CorrectorTag_;
  const std::string ECorrectorTag_;

private:
  inline std::string getFilePath(const edm::ParameterSet & pset, const std::string & name) const {
    std::string ret = pset.getParameter<std::string>(name);
    if (ret[0] != '/') ret = edm::FileInPath(ret).fullPath();
    return ret;
  }

  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void genMatch(std::vector<double> &iGenVars, int iType, double iEta, double iPhi, double iPt, const reco::GenParticleCollection &iGenParticles);
  void addPF(const L1PFCollection &iCandidates, const std::string &iLabel, edm::Event& iEvent);

  edm::EDGetTokenT<reco::GenParticleCollection>   TokGenPar_;
  std::vector<edm::EDGetTokenT<L1PFCollection>>   TokEcalTPTags_;
  std::vector<edm::EDGetTokenT<L1PFCollection>>   TokHcalTPTags_;
  std::unique_ptr<corrector> corrector_;
  std::unique_ptr<corrector> ecorrector_;
  // simplified corrections
  l1tpf::SimpleCorrEm simpleCorrEm_;
  l1tpf::SimpleCorrHad simpleCorrHad_;

  // new calo clusterer (float)
  l1pf_calo::SingleCaloClusterer ecalClusterer_, hcalClusterer_;
  l1pf_calo::SimpleCaloLinker caloLinker_;
     
  // declare variables for output file
  std::string fOutputName;
  TFile *fOutputFile;
  TH1D  *fTotalEvents;

  TTree *fEcalInfoTree;
  TTree *fHcalInfoTree;
  float runNum, lumiSec, evtNum;
  float genPt, genEta, genPhi, genId;

  float ecal_subdet, ecal_ieta, ecal_iphi, ecal_curTwrEta, ecal_curTwrPhi, ecal_et, ecal_num,ecal_dr;
  float hcal_subdet, hcal_ieta, hcal_iphi, hcal_TwrR, hcal_et, hcal_num, hcal_ecal_et,hcal_ecal_etcorr,hcal_ecal_eta,hcal_ecal_phi,hcal_dr;
  float ecal_clust_et,ecal_clust_eta,ecal_clust_phi,ecal_corr_et;
  float hcal_clust_et,hcal_clust_eta,hcal_clust_phi,hcal_clust_emf,hcal_corr_et,hcal_corr_emf;
  float ecal_genPt, ecal_genEta, ecal_genPhi, ecal_genId;
  float hcal_genPt, hcal_genEta, hcal_genPhi, hcal_genId;
};

//
// constructors and destructor
//
CaloNtupleProducer::CaloNtupleProducer(const edm::ParameterSet& iConfig):
  zeroSuppress_         (iConfig.getParameter<bool>("zeroSuppress")),
  ecalOnly_             (iConfig.existsAs<bool>("ecalOnly") ? iConfig.getParameter<bool>("ecalOnly") : false),
  EcalTPTags_           (iConfig.getParameter<std::vector<edm::InputTag>>("EcalTPTags")),
  HcalTPTags_           (!ecalOnly_ ? iConfig.getParameter<std::vector<edm::InputTag>>("HcalTPTags") : std::vector<edm::InputTag>()),
  GenParTag_            (iConfig.getParameter<edm::InputTag>("genParTag")),
  CorrectorTag_         (getFilePath(iConfig,"corrector")),
  ECorrectorTag_        (getFilePath(iConfig,"ecorrector")),
  corrector_            (new corrector(CorrectorTag_,11,iConfig.getUntrackedParameter<int>("debug",0))),
  ecorrector_           (new corrector(ECorrectorTag_,1,iConfig.getUntrackedParameter<int>("debug",0))),
  simpleCorrEm_         (iConfig, "simpleCorrEm"),
  simpleCorrHad_        (iConfig, "simpleCorrHad"),
  ecalClusterer_        (iConfig.getParameter<edm::ParameterSet>("caloClusterer").getParameter<edm::ParameterSet>("ecal")),
  hcalClusterer_        (iConfig.getParameter<edm::ParameterSet>("caloClusterer").getParameter<edm::ParameterSet>("hcal")),
  caloLinker_           (iConfig.getParameter<edm::ParameterSet>("caloClusterer").getParameter<edm::ParameterSet>("linker"), ecalClusterer_, hcalClusterer_),
  fOutputName           (iConfig.getUntrackedParameter<std::string>("outputName", "ntuple.root")),
  fOutputFile           (0),
  fTotalEvents          (0),
  fEcalInfoTree         (0),
  fHcalInfoTree         (0)
{
  TokGenPar_ = consumes<reco::GenParticleCollection>(GenParTag_);
  for (const edm::InputTag &tag : EcalTPTags_) {
    TokEcalTPTags_.push_back(consumes<L1PFCollection>(tag));
  }
  for (const edm::InputTag &tag : HcalTPTags_) {
    TokHcalTPTags_.push_back(consumes<L1PFCollection>(tag));
  }
  produces<L1PFCollection>("emCalibrated");
  produces<L1PFCollection>("emUncalibrated");

  if (ecalOnly_) return;

  produces<L1PFCollection>("uncalibrated");
  produces<L1PFCollection>("calibrated");
  produces<PFOutputCollection>("RawCalo");
  produces<PFOutputCollection>("Calo");
}

CaloNtupleProducer::~CaloNtupleProducer()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
CaloNtupleProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //fOutputFile->cd();
  if (!fOutputName.empty()) fTotalEvents->Fill(1);  
  using namespace edm;

  edm::Handle<reco::GenParticleCollection> hGenParProduct;
  iEvent.getByToken(TokGenPar_,hGenParProduct);
  const reco::GenParticleCollection genParticles = *(hGenParProduct.product());  

  ecalClusterer_.clear();
  hcalClusterer_.clear();
  /// ----------------ECAL INFO-------------------
  edm::Handle<L1PFCollection> ecals;
  for (const auto & token : TokEcalTPTags_) {
      iEvent.getByToken(token, ecals);
      for (const l1tpf::Particle & it : *ecals) {
          ecalClusterer_.add(it);
      }
  }

  ecalClusterer_.run();
  //// ---- Dummy calibration == no calibration
  // ecalClusterer_.correct( [](const l1pf_calo::Cluster &c, int ieta, int iphi) -> double { return c.et; } );
  //    
  //// ---- Calibration from Phil's workflow
  if (simpleCorrEm_.empty()) {
      ecalClusterer_.correct( [&](const l1pf_calo::Cluster &c, int ieta, int iphi) -> double { 
              return ecorrector_->correct(0., c.et, ieta, iphi);
              } );
  } else {
      ecalClusterer_.correct( [&](const l1pf_calo::Cluster &c, int ieta, int iphi) -> double { 
              return simpleCorrEm_(c.et, std::abs(c.eta));
              } );
  }
  // write debug output tree
  if (!fOutputName.empty()) {
      unsigned int ne = 0;
      const auto & ecraw = ecalClusterer_.raw();
      const auto & ecals = ecalClusterer_.clusters();
      for (unsigned int i = 0, ncells = ecals.size(); i < ncells; ++i) {
          if (ecals[i].et == 0) continue; 
          ecal_num = ne++;
          ecal_et = ecraw[i];
          ecal_ieta = ecals.ieta(i);
          ecal_iphi = ecals.iphi(i);
          ecal_curTwrEta = ecals.eta(i);
          ecal_curTwrPhi = ecals.phi(i);
          ecal_clust_et = ecals[i].et;
          ecal_clust_eta = ecals[i].eta;
          ecal_clust_phi = ecals[i].phi;
          ecal_corr_et = ecals[i].et_corr;

          std::vector<double> lGenVars;
          genMatch(lGenVars,1,ecal_clust_eta,ecal_clust_phi,ecal_clust_et,genParticles);
          ecal_genPt=0; ecal_genEta=0; ecal_genPhi=0; ecal_genId=0; ecal_dr = 0;
          if(lGenVars.size() > 3) { 
              ecal_genPt   = float(lGenVars[0]);
              ecal_genEta  = float(lGenVars[1]);
              ecal_genPhi  = float(lGenVars[2]);
              ecal_genId   = float(lGenVars[3]);
              ecal_dr = reco::deltaR( ecal_genEta, ecal_genPhi, ecal_clust_eta,ecal_clust_phi );
          }
          if (zeroSuppress_) {
              if(ecal_genPt > 1.) fEcalInfoTree->Fill();      
          } else {
              if(ecal_et > 1.) fEcalInfoTree->Fill();      
          }
      }
  }

  // Get also ecal-only from the clusterer
  std::unique_ptr<L1PFCollection> RawEcal  = ecalClusterer_.fetch(false);
  std::unique_ptr<L1PFCollection> CorrEcal = ecalClusterer_.fetch(true);

  iEvent.put(std::move(RawEcal),  "emUncalibrated");
  iEvent.put(std::move(CorrEcal), "emCalibrated");
  if (ecalOnly_) return;


  // / ----------------HCAL INFO-------------------

  edm::Handle<L1PFCollection> hcals;
  for (const auto & token : TokHcalTPTags_) {
      iEvent.getByToken(token, hcals);
      for (const l1tpf::Particle & it : *hcals) {
          hcalClusterer_.add(it);
      }
  }
  hcalClusterer_.run();

  // Calorimeter linking
  caloLinker_.run();
  //// ---- Dummy calibration (no calibration at all)
  // caloLinker_.correct( [](const l1pf_calo::CombinedCluster &c, int ieta, int iphi) -> double { return c.et; } );
  //
  //// ---- Calibration from Phil's workflow
  if (simpleCorrHad_.empty()) {
      caloLinker_.correct( [&](const l1pf_calo::CombinedCluster &c, int ieta, int iphi) -> double { 
              return corrector_->correct(c.et, c.ecal_et, ieta, iphi); 
              } );
  } else {
      caloLinker_.correct( [&](const l1pf_calo::CombinedCluster &c, int ieta, int iphi) -> double { 
              return simpleCorrHad_(c.et, std::abs(c.eta), c.ecal_et/c.et); 
              } );
  }
  // write debug output tree
  if (!fOutputName.empty()) {
      const auto & clusters = caloLinker_.clusters();
      unsigned int nh = 0;
      for (unsigned int i = 0, ncells = clusters.size(); i < ncells; ++i) {
          if (clusters[i].et == 0) continue; 
          hcal_num = nh++;
          hcal_et = clusters[i].hcal_et;
          hcal_ecal_et = clusters[i].ecal_et;
          hcal_ecal_eta = -999; // FIXME missing
          hcal_ecal_phi = -999; // FIXME missing
          hcal_ieta = clusters.ieta(i);
          hcal_iphi = clusters.iphi(i);
          hcal_clust_et  = clusters[i].et;  // all these are of the combined
          hcal_clust_eta = clusters[i].eta; // cluster (ecal+hcal)
          hcal_clust_phi = clusters[i].phi;
          hcal_clust_emf = clusters[i].ecal_et / clusters[i].et; // note: this is raw EMF
          hcal_corr_et  = clusters[i].et_corr;
          hcal_corr_emf = -999; // FIXME
          std::vector<double> lGenVars;
          genMatch(lGenVars,1,hcal_clust_eta,hcal_clust_phi,hcal_clust_et,genParticles);
          hcal_genPt=0; hcal_genEta=0; hcal_genPhi=0; hcal_genId=0; hcal_dr = 0;
          if(lGenVars.size() > 3) { 
              hcal_genPt   = float(lGenVars[0]);
              hcal_genEta  = float(lGenVars[1]);
              hcal_genPhi  = float(lGenVars[2]);
              hcal_genId   = float(lGenVars[3]);
              hcal_dr = reco::deltaR( hcal_genEta, hcal_genPhi, hcal_clust_eta,hcal_clust_phi );
          }
          if (zeroSuppress_) {
              if(hcal_genPt > 1.) fHcalInfoTree->Fill();      
          } else {
              if(hcal_et > 1.) fHcalInfoTree->Fill();      
          }
      }
  }

  // Get particles from the clusterer
  std::unique_ptr<L1PFCollection> RawCalo  = caloLinker_.fetch(false);
  std::unique_ptr<L1PFCollection> CorrCalo = caloLinker_.fetch(true);

  addPF(*RawCalo,  "RawCalo",  iEvent);
  addPF(*CorrCalo, "Calo",     iEvent);

  iEvent.put(std::move(RawCalo),  "uncalibrated");
  iEvent.put(std::move(CorrCalo), "calibrated");


}

void CaloNtupleProducer::addPF(const L1PFCollection &iCandidates, const std::string &iLabel, edm::Event& iEvent) { 
  std::unique_ptr<PFOutputCollection > corrCandidates( new PFOutputCollection );
  corrCandidates->reserve(iCandidates.size());
  for(unsigned int i0 = 0; i0 < iCandidates.size(); i0++) { 
    reco::PFCandidate::ParticleType id = reco::PFCandidate::ParticleType::X; 
    int pCharge=0; 
    if(iCandidates[i0].pdgId() == l1tpf::Particle::CH) id = reco::PFCandidate::ParticleType::h;
    if(iCandidates[i0].pdgId() == l1tpf::Particle::EL) id = reco::PFCandidate::ParticleType::e;
    if(iCandidates[i0].pdgId() == l1tpf::Particle::NH) id = reco::PFCandidate::ParticleType::h0;
    if(iCandidates[i0].pdgId() == l1tpf::Particle::GAMMA) id = reco::PFCandidate::ParticleType::gamma;
    if(iCandidates[i0].pdgId() == l1tpf::Particle::MU) id = reco::PFCandidate::ParticleType::mu;
    if(iCandidates[i0].pdgId() == l1tpf::Particle::CH || iCandidates[i0].pdgId() == l1tpf::Particle::EL)  pCharge = 1;
    if(iCandidates[i0].pdgId() == l1tpf::Particle::MU)  pCharge = iCandidates[i0].charge();
    reco::PFCandidate pCand(pCharge,iCandidates[i0].p4(),id);
    corrCandidates->push_back(pCand);
  }
  //Fill!
  iEvent.put(std::move(corrCandidates),iLabel);
}

//////////////////////////// ------------------------------------------------------
//////////////////////////// ------------------------------------------------------

//--- Gen Matching
void CaloNtupleProducer::genMatch(std::vector<double> &iGenVars, int iType, double iEta, double iPhi, double iPt, const reco::GenParticleCollection &iGenParticles) { 
  int lId = -999;
  double lPt = -1;
  double lDeltaRMin = 100;
  TLorentzVector lVec; lVec.SetPtEtaPhiM(0,0,0,0);
  for (reco::GenParticleCollection::const_iterator itGenP = iGenParticles.begin(); itGenP!=iGenParticles.end(); ++itGenP) {
    if(iType == 0 && itGenP->charge() == 0) continue;
    double deltaEta = itGenP->eta()-iEta;
    double deltaPhi = fabs(itGenP->phi()-iPhi); if(deltaPhi > 2.*TMath::Pi()-deltaPhi) deltaPhi = 2.*TMath::Pi()-deltaPhi;
    double deltaR   = sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);
    if(deltaR > 0.2) continue;
    if(iType == 0 && deltaR > lDeltaRMin) continue;
    if(iType == 0) lId = itGenP->pdgId();
    if(iType == 1 && itGenP->pt() > lPt) { 
      lPt = itGenP->pt();
      lId = itGenP->pdgId();
    }
    lDeltaRMin = deltaR;
    TLorentzVector pVec; pVec.SetPtEtaPhiM(itGenP->pt(),itGenP->eta(),itGenP->phi(),itGenP->mass());
    if(iType == 0) lVec  = pVec;
    if(iType == 1) lVec += pVec;
  }
  if(lVec.Pt() > 0) iGenVars.push_back(lVec.Pt());
  if(lVec.Pt() > 0) iGenVars.push_back(lVec.Eta());
  if(lVec.Pt() > 0) iGenVars.push_back(lVec.Phi());
  if(lVec.Pt() > 0) iGenVars.push_back(lId);
}

// ------------ method called once each job just before starting event loop  ------------
void 
CaloNtupleProducer::beginJob()
{
  if (fOutputName.empty()) return;
  //
  // Create output file, trees, and histograms
  //
  fOutputFile = new TFile(fOutputName.c_str(), "RECREATE");
  fTotalEvents = new TH1D("TotalEvents","TotalEvents",1,-10,10);
  fEcalInfoTree     = new TTree("EcalInfo",   "EcalInfo");
  fHcalInfoTree     = new TTree("HcalInfo",   "HcalInfo");

  fEcalInfoTree->Branch("ecal_ieta", &ecal_ieta, "ecal_ieta/F");
  fEcalInfoTree->Branch("ecal_iphi", &ecal_iphi, "ecal_iphi/F");
  fEcalInfoTree->Branch("ecal_curTwrEta", &ecal_curTwrEta, "ecal_curTwrEta/F");
  fEcalInfoTree->Branch("ecal_curTwrPhi", &ecal_curTwrPhi, "ecal_curTwrPhi/F");
  fEcalInfoTree->Branch("ecal_et", &ecal_et, "ecal_et/F");
  fEcalInfoTree->Branch("ecal_num",&ecal_num, "ecal_num/F");
  fEcalInfoTree->Branch("ecal_clust_et",  &ecal_clust_et , "ecal_clust_et/F");
  fEcalInfoTree->Branch("ecal_clust_eta", &ecal_clust_eta, "ecal_clust_eta/F");
  fEcalInfoTree->Branch("ecal_clust_phi", &ecal_clust_phi, "ecal_clust_phi/F");
  fEcalInfoTree->Branch("ecal_corr_et"  , &ecal_corr_et  , "ecal_corr_et/F");
  fEcalInfoTree->Branch("genPt",   &ecal_genPt,"ecal_genPt/F");
  fEcalInfoTree->Branch("genEta",  &ecal_genEta,"ecal_genEta/F");
  fEcalInfoTree->Branch("genPhi",  &ecal_genPhi,"ecal_genPhi/F");
  fEcalInfoTree->Branch("genid",   &ecal_genId, "ecal_genid/F");
  fEcalInfoTree->Branch("gendr",   &ecal_dr, "ecal_dr/F");


  //fHcalInfoTree->Branch("hcal_subdet", &hcal_subdet, "hcal_subdet/F");
  fHcalInfoTree->Branch("hcal_ieta", &hcal_ieta, "hcal_ieta/F");
  fHcalInfoTree->Branch("hcal_iphi", &hcal_iphi, "hcal_iphi/F");
  //fHcalInfoTree->Branch("hcal_TwrR", &hcal_TwrR, "hcal_TwrR/F");
  fHcalInfoTree->Branch("hcal_num", &hcal_num, "hcal_num/F");
  fHcalInfoTree->Branch("hcal_et", &hcal_et, "hcal_et/F");
  fHcalInfoTree->Branch("hcal_clust_et",  &hcal_clust_et , "hcal_clust_et/F");
  fHcalInfoTree->Branch("hcal_clust_eta", &hcal_clust_eta, "hcal_clust_eta/F");
  fHcalInfoTree->Branch("hcal_clust_phi", &hcal_clust_phi, "hcal_clust_phi/F");
  fHcalInfoTree->Branch("hcal_clust_emf" , &hcal_clust_emf , "hcal_clust_emf/F");
  fHcalInfoTree->Branch("hcal_corr_et"  , &hcal_corr_et  , "hcal_corr_et/F");
  fHcalInfoTree->Branch("hcal_corr_emf" , &hcal_corr_emf , "hcal_corr_emf/F");
  fHcalInfoTree->Branch("ecal_et",      &hcal_ecal_et,     "hcal_ecal_et/F");
  fHcalInfoTree->Branch("ecal_etcorr",  &hcal_ecal_etcorr, "hcal_ecal_etcorr/F");
  //fHcalInfoTree->Branch("ecal_eta", &hcal_ecal_eta, "hcal_ecal_eta/F"); // FIXME
  //fHcalInfoTree->Branch("ecal_phi", &hcal_ecal_phi, "hcal_ecal_phi/F"); // FIXME
  fHcalInfoTree->Branch("genPt",   &hcal_genPt ,"hcal_genPt/F");
  fHcalInfoTree->Branch("genEta",  &hcal_genEta,"hcal_genEta/F");
  fHcalInfoTree->Branch("genPhi",  &hcal_genPhi,"hcal_genPhi/F");
  fHcalInfoTree->Branch("genid",   &hcal_genId, "hcal_genid/F");
  fHcalInfoTree->Branch("gendr",   &hcal_dr ,"hcal_dr/F");  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CaloNtupleProducer::endJob() {
  if (fOutputName.empty()) return;
  //
  // Save to ROOT file
  //
  fOutputFile->cd();
  fTotalEvents->Write();
  fOutputFile->Write();
  fOutputFile->Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CaloNtupleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CaloNtupleProducer);
