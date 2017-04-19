// -*- C++ -*-
//
// Package:    NtupleProducer
// Class:      NtupleProducer
// 
/**\class NtupleProducer NtupleProducer.cc Ntuplizer/NtupleProducer/plugins/NtupleProducer.cc

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
#include "FastPUPPI/NtupleProducer/interface/combiner.hh"
#include "FastPUPPI/NtupleProducer/interface/metanalyzer.hh"
#include "FastPUPPI/NtupleProducer/interface/jetanalyzer.hh"
#include "FastPUPPI/NtupleProducer/interface/L1TPFUtils.h"
#include "FastPUPPI/NtupleProducer/interface/DiscretePF.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"
#include "DataFormats/Math/interface/deltaR.h"


const double PI = 3.1415926535897;

// ROOT classes
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObject.h>

typedef math::XYZTLorentzVector                        LorentzVector;
typedef std::vector< reco::PFCandidate >               PFOutputCollection;
typedef std::vector<l1tpf::Particle>                   L1PFCollection;

//--------------------------------------------------------------------------------------------------
class NtupleProducer : public edm::EDProducer {
public:
  explicit NtupleProducer(const edm::ParameterSet&);
  ~NtupleProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  const bool          zeroSuppress_;
  const edm::InputTag L1TrackTag_;
  const edm::InputTag EcalTPTag_;
  const edm::InputTag HGEcalTPTag_;
  const edm::InputTag HcalTPTag_;
  const edm::InputTag HGHcalTPTag_;
  const edm::InputTag BHHcalTPTag_;
  const edm::InputTag HFTPTag_;
  const edm::InputTag MuonTPTag_;
  const edm::InputTag GenParTag_;
  const std::string CorrectorTag_;
  const std::string Corrector2Tag_;
  const std::string ECorrectorTag_;
  const std::string TrackResTag_;
  const std::string EleResTag_;
  const std::string PionResTag_;
  std::unique_ptr<PFOutputCollection > corrCandidates_;
 
private:
  inline std::string getFilePath(const edm::ParameterSet & pset, const std::string & name) const {
    std::string ret = pset.getParameter<std::string>(name);
    if (ret[0] != '/') ret = edm::FileInPath(ret).fullPath();
    return ret;
  }

  struct MyEcalCluster { 
    int ieta, iphi; float et, corr_et, eta, phi; 
    MyEcalCluster(int iIeta, int iIphi, float iEt, float iCorr_et, float iEta, float iPhi) :
        ieta(iIeta), iphi(iIphi), et(iEt), corr_et(iCorr_et), eta(iEta), phi(iPhi) {}
  };

  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void genMatch(std::vector<double> &iGenVars,int iType,double iEta,double iPhi,double iPt,const reco::GenParticleCollection &iGenParticles);
  TLorentzVector getVector(double iPt[][73],int iEta,int iPhi,int iEta0,int iPhi0,double iNSigma=2,double iENoise=0.2);
  void simpleCluster(std::vector<TLorentzVector> &iClusters,double  iEta,double iPhi,double iPt[][73],double iNSigma=2,double iENoise=0.2);
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  void addPF(std::vector<combiner::Particle> &iCandidates,std::string iLabel,edm::Event& iEvent);

  edm::EDGetTokenT<reco::GenParticleCollection>   TokGenPar_;
  edm::EDGetTokenT<L1PFCollection>                TokL1TrackTPTag_;
  edm::EDGetTokenT<L1PFCollection>                TokEcalTPTag_;
  edm::EDGetTokenT<L1PFCollection>                TokHGEcalTPTag_;
  edm::EDGetTokenT<L1PFCollection>                TokHcalTPTag_;
  edm::EDGetTokenT<L1PFCollection>                TokHGHcalTPTag_;
  edm::EDGetTokenT<L1PFCollection>                TokBHHcalTPTag_;
  edm::EDGetTokenT<L1PFCollection>                TokHFTPTag_;
  edm::EDGetTokenT<l1t::MuonBxCollection>         TokMuonTPTag_;
  double trkPt_;
  bool   metRate_;
  double etaCharged_;
  double puppiPtCut_;
  double vtxRes_;
  corrector* corrector_;
  corrector* corrector2_;
  corrector* ecorrector_;
  combiner * connector_;
  combiner * rawconnector_;
  metanalyzer* metanalyzer_;
  jetanalyzer* jetanalyzer_;
  // discretized version
  l1tpf_int::RegionMapper l1regions_;
  l1tpf_int::PFAlgo       l1pfalgo_;
  // debug flag
  int fDebug;
     
  // declare variables for output file
  std::string           fOutputName;
  TFile                 *fOutputFile;
  TH1D                  *fTotalEvents;

  TTree                 *fTrkInfoTree;
  TTree                 *fEcalInfoTree;
  TTree                 *fHcalInfoTree;
  float runNum, lumiSec, evtNum;
  float trkNum;
  float trkPx, trkPz, trkPy, trkPt, trkEta, trkPhi, trkz0, trkd0;    
  float trkEcalEta, trkEcalPhi, trkEcalR;
  float genPt, genEta, genPhi, genId;

  float ecal_subdet, ecal_ieta, ecal_iphi, ecal_curTwrEta, ecal_curTwrPhi, ecal_et, ecal_num,ecal_dr;
  float hcal_subdet, hcal_ieta, hcal_iphi, hcal_TwrR, hcal_et, hcal_num, hcal_ecal_et,hcal_ecal_etcorr,hcal_ecal_eta,hcal_ecal_phi,hcal_dr;
  float ecal_clust_et,ecal_clust_eta,ecal_clust_phi,ecal_corr_et;
  float hcal_clust_et,hcal_clust_eta,hcal_clust_phi,hcal_corr_et,hcal_corr_emf;
  float ecal_genPt, ecal_genEta, ecal_genPhi, ecal_genId;
  float hcal_genPt, hcal_genEta, hcal_genPhi, hcal_genId;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;  
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
NtupleProducer::NtupleProducer(const edm::ParameterSet& iConfig):
  zeroSuppress_         (iConfig.getParameter<bool>("zeroSuppress")),
  L1TrackTag_           (iConfig.getParameter<edm::InputTag>("L1TrackTag")),
  EcalTPTag_            (iConfig.getParameter<edm::InputTag>("EcalTPTag")),
  HGEcalTPTag_          (iConfig.getParameter<edm::InputTag>("HGEcalTPTag")),
  HcalTPTag_            (iConfig.getParameter<edm::InputTag>("HcalTPTag")),
  HGHcalTPTag_          (iConfig.getParameter<edm::InputTag>("HGHcalTPTag")),
  BHHcalTPTag_          (iConfig.getParameter<edm::InputTag>("BHHcalTPTag")),
  HFTPTag_              (iConfig.getParameter<edm::InputTag>("HFTPTag")),
  MuonTPTag_            (iConfig.getParameter<edm::InputTag>("MuonTPTag")),
  GenParTag_            (iConfig.getParameter<edm::InputTag>("genParTag")),
  CorrectorTag_         (getFilePath(iConfig,"corrector")),
  Corrector2Tag_        (getFilePath(iConfig,"corrector2")),
  ECorrectorTag_        (getFilePath(iConfig,"ecorrector")),
  TrackResTag_          (getFilePath(iConfig,"trackres")),
  EleResTag_            (getFilePath(iConfig,"eleres")),
  PionResTag_           (getFilePath(iConfig,"pionres")),
  trkPt_                (iConfig.getParameter<double>       ("trkPtCut")),
  metRate_              (iConfig.getParameter<bool>         ("metRate")),
  etaCharged_           (iConfig.getParameter<double>       ("etaCharged")),
  puppiPtCut_           (iConfig.getParameter<double>       ("puppiPtCut")),
  vtxRes_               (iConfig.getParameter<double>       ("vtxRes")),
  l1regions_            (iConfig),
  l1pfalgo_             (iConfig),
  fDebug                (iConfig.getUntrackedParameter<int>("debug",0)),
  fOutputName           (iConfig.getUntrackedParameter<std::string>("outputName", "ntuple.root")),
  fOutputFile           (0),
  fTotalEvents          (0),
  fTrkInfoTree          (0),
  fEcalInfoTree         (0),
  fHcalInfoTree         (0)
{
  //now do what ever other initialization is needed
  corrector_  = new corrector(CorrectorTag_,11,fDebug);
  corrector2_ = new corrector(Corrector2Tag_,11,fDebug);
  ecorrector_ = new corrector(ECorrectorTag_,1,fDebug);
  connector_  = new combiner (PionResTag_,EleResTag_,TrackResTag_,"puppi.root",etaCharged_,puppiPtCut_,vtxRes_,fDebug);
  rawconnector_  = new combiner (PionResTag_,EleResTag_,TrackResTag_,"puppiraw.root",etaCharged_,puppiPtCut_,vtxRes_);
  metanalyzer_   = new metanalyzer("MetFile.root");
  jetanalyzer_   = new jetanalyzer("JetFile.root");
  produces<PFOutputCollection>("TK");
  produces<PFOutputCollection>("RawCalo");
  produces<PFOutputCollection>("Calo");
  produces<PFOutputCollection>("PF");
  produces<PFOutputCollection>("Puppi");
  produces<PFOutputCollection>("L1TK");
  produces<PFOutputCollection>("L1Calo");
  produces<PFOutputCollection>("L1PF");
  produces<PFOutputCollection>("L1Puppi");
  TokGenPar_       = consumes<reco::GenParticleCollection>( GenParTag_    );
  TokL1TrackTPTag_ = consumes<L1PFCollection>( L1TrackTag_  );
  TokEcalTPTag_    = consumes<L1PFCollection>( EcalTPTag_   );
  TokHGEcalTPTag_  = consumes<L1PFCollection>( HGEcalTPTag_ );
  TokHcalTPTag_    = consumes<L1PFCollection>( HcalTPTag_   );
  TokHGHcalTPTag_  = consumes<L1PFCollection>( HGHcalTPTag_ );
  TokBHHcalTPTag_  = consumes<L1PFCollection>( BHHcalTPTag_ );
  TokHFTPTag_      = consumes<L1PFCollection>( HFTPTag_     );
  TokMuonTPTag_    = consumes<l1t::MuonBxCollection>( MuonTPTag_  );
}

NtupleProducer::~NtupleProducer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
NtupleProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //fOutputFile->cd();
  fTotalEvents->Fill(1);  
  using namespace edm;

  edm::Handle<reco::GenParticleCollection> hGenParProduct;
  iEvent.getByToken(TokGenPar_,hGenParProduct);
  const reco::GenParticleCollection genParticles = *(hGenParProduct.product());  

  connector_->clear();
  rawconnector_->clear();
  l1regions_.clear(); 
  /// ----------------TRACK INFO-------------------
  edm::Handle<std::vector<l1tpf::Particle>> l1tks;
  iEvent.getByLabel(L1TrackTag_, l1tks);
  trkNum = 0;
  for (l1tpf::Particle tk : *l1tks) { // no const & since we modify it to set the sigma
      tk.setSigma(connector_->getTrkRes(tk.pt(), tk.eta(), tk.phi())); // the combiner knows the sigma, the track producer doesn't
      // adding objects to PF
      if(tk.pt() > trkPt_) l1regions_.addTrack(tk);      
      if(tk.pt() > trkPt_) connector_->addTrack(tk);      
      if(tk.pt() > trkPt_) rawconnector_->addTrack(tk);
      /// filling the tree    
      trkPx  = tk.px();
      trkPz  = tk.py();
      trkPy  = tk.pz();
      trkPt  = tk.pt();
      trkEta = tk.eta();
      trkPhi = tk.phi();
      trkz0  = tk.vertex().Z();
      trkd0  = tk.vertex().Rho();
      trkEcalEta = tk.caloEta();
      trkEcalPhi = tk.caloPhi();
      trkEcalR   = -999.; //FIXME: this wasn't in the Particle class. Do we need it?
      std::vector<double> lGenVars;
      genMatch(lGenVars,0,double(trkEta),double(trkPhi),double(trkPt),genParticles);
      genPt=0; genEta=0; genPhi=0; genId=0;
      if(lGenVars.size() > 3) { 
      	genPt   = float(lGenVars[0]);
      	genEta  = float(lGenVars[1]);
      	genPhi  = float(lGenVars[2]);
      	genId   = float(lGenVars[3]);
      }
      fTrkInfoTree->Fill();   
      trkNum++;
  }
  
  /// ----------------Muon INFO-------------------
  // new muon getting in 91x
  edm::Handle<l1t::MuonBxCollection> muon;    
  iEvent.getByToken(TokMuonTPTag_, muon);
  if (muon.isValid()){ 
    for (int ibx = muon->getFirstBX(); ibx <= muon->getLastBX(); ++ibx) {
      if (ibx != 0) continue; // only the first bunch crossing
      for (auto it=muon->begin(ibx); it!=muon->end(ibx); it++){      
        if (it->et() == 0) continue; // if you don't care about L1T candidates with zero ET.
        l1tpf::Particle mu( it->pt(), it->eta(), ::deltaPhi(it->phi(), 0.), 0.105, combiner::MU ); // the deltaPhi is to get the wrapping correct
        mu.setCharge(it->charge());
        mu.setQuality(it->hwQual());
        mu.setSigma(connector_->getTrkRes(mu.pt(),mu.eta(),mu.phi())); // this needs to be updated with the muon resolutions!      

        connector_->addMuon(mu);
        l1regions_.addMuon(mu);
      }
    }
  } 
  else {
    edm::LogWarning("MissingProduct") << "L1Upgrade muon bx collection not found." << std::endl;
  }

  /// ----------------ECAL INFO-------------------
  std::vector<MyEcalCluster> lEcals;  
  double lEEt[l1tpf::towerNEta()][73];
  double lHEt[l1tpf::towerNEta()][73];
  for(int i0 = 0; i0 < l1tpf::towerNEta(); i0++) { for(int i1 = 0; i1 <l1tpf::towerNPhi(1)+1; i1++) {lEEt[i0][i1]=0; lHEt[i0][i1] = 0;}}
  edm::Handle<L1PFCollection> classicecals;
  iEvent.getByToken(TokEcalTPTag_, classicecals);
  edm::Handle<L1PFCollection> hgecals;
  iEvent.getByToken(TokHGEcalTPTag_, hgecals);
  L1PFCollection * ecals = new L1PFCollection;
  for (const l1tpf::Particle & it : *classicecals) ecals->push_back(it); 
  for (const l1tpf::Particle & it : *hgecals     ) ecals->push_back(it); 
  for (const l1tpf::Particle & it : *ecals) {
     lEEt[it.aEta()][it.aPhi()] = it.pt();
  }
  std::vector<TLorentzVector> pClust;
  int ne = 0;
  
  for (const l1tpf::Particle & it : *ecals) {
      double et = it.pt();
      ecal_subdet = 0;// it->id().subDet(); // FIXME: missing in Particle class
      ecal_ieta = it.iEta();
      ecal_iphi = it.iPhi();
      ecal_et = et;
      ecal_curTwrEta = it.caloEta();
      ecal_curTwrPhi = it.caloPhi();
      ecal_num = ne;
      simpleCluster(pClust,it.aEta(),it.aPhi(),lEEt);
      // FIXME Gio: shouldn't we continue to the next TP here if we don't find a cluster?
      // Phil : Yes, except this ntuple also contains debug info so we can look at energy lost from the cluster (You are right we don't use it)
      ecal_clust_et=0,ecal_clust_eta=0,ecal_clust_phi =0;
      if(pClust.size() > 0) ecal_clust_et =pClust[0].Pt();
      if(pClust.size() > 0) ecal_clust_eta=pClust[0].Eta();
      if(pClust.size() > 0) ecal_clust_phi=pClust[0].Phi();
      pClust.clear();
      ecal_corr_et = ecorrector_->correct(0.,double(ecal_clust_et),ecal_ieta,ecal_iphi);
      //if(ecal_clust_et > 0) std::cout << "=ECalo=> " << ecal_clust_et << "---> " << ecal_corr_et << " ---> " << ecal_et << " -- " << ecal_clust_eta << " -- " << ecal_clust_phi << std::endl;
      if(ecal_clust_et > 0)  { 
	lEcals.push_back(MyEcalCluster(it.iEta(), it.iPhi(), ecal_clust_et, ecal_corr_et, ecal_clust_eta, ecal_clust_phi));
      }
      std::vector<double> lGenVars;
      genMatch(lGenVars,1,double(it.caloEta()),double(it.caloPhi()),double(et),genParticles);
      ecal_genPt=0; ecal_genEta=0; ecal_genPhi=0; ecal_genId=0; ecal_dr = 0;
      if(lGenVars.size() > 3) { 
      	ecal_genPt   = float(lGenVars[0]);
      	ecal_genEta  = float(lGenVars[1]);
      	ecal_genPhi  = float(lGenVars[2]);
      	ecal_genId   = float(lGenVars[3]);
	ecal_dr = reco::deltaR( ecal_genEta, ecal_genPhi, ecal_clust_eta,ecal_clust_phi );
      }
      if(ecal_genPt > 1. && zeroSuppress_) fEcalInfoTree->Fill();      
      //if(!zeroSuppress_) fEcalInfoTree->Fill();
      if(!zeroSuppress_ && ecal_et > 1.) fEcalInfoTree->Fill();
      ne++;
      // if (ne > 20) break;
  }
  // / ----------------HCAL INFO-------------------
  // / Stealing some other code!
  // / https://github.com/cms-sw/cmssw/blob/0397259dd747cee94b68928f17976224c037057a/L1Trigger/L1TNtuples/src/L1AnalysisCaloTP.cc#L40
  // / 72 in phi and 56 in eta = 4032
  // / 144 TP is HF: 4 (in eta) * 18 (in phi) * 2 (sides)
  // / Hcal TPs
  edm::Handle<L1PFCollection> classichcals;
  iEvent.getByToken(TokHcalTPTag_, classichcals);
  edm::Handle<L1PFCollection> hghcals;
  iEvent.getByToken(TokHGHcalTPTag_, hghcals);
  edm::Handle<L1PFCollection> bhhcals;
  iEvent.getByToken(TokBHHcalTPTag_, bhhcals);
  edm::Handle<L1PFCollection> hfhcals;
  iEvent.getByToken(TokHFTPTag_, hfhcals);
  L1PFCollection * hcals = new L1PFCollection;
  for (const l1tpf::Particle & it : *classichcals) hcals->push_back(it); 
  for (const l1tpf::Particle & it : *hfhcals     ) hcals->push_back(it); 
  for (const l1tpf::Particle & it : *hghcals     ) {
    //combine duplicates for towers we use
    bool lFound = false;
    for (l1tpf::Particle & it2 : *hcals     ) {
      if(it2.eta() == it.eta() && it2.phi() == it.phi() && it2.pt() == it.pt()) {std::cout << "copied!!!!!!!!!!!!!!!!!" << std::endl; break;}
      if(it2.iEta() != it.iEta() || it2.iPhi() != it.iPhi()) continue;
      TLorentzVector lVec = it.tp4();
      lVec += it2.tp4();
      it2.setPtEtaPhiM(lVec.Pt(),lVec.Eta(),lVec.Phi(),lVec.M());
      lFound = true; 
      break;
      }
    if(!lFound) hcals->push_back(it);
  } 
  for (const l1tpf::Particle & it : *bhhcals     ) { 
    //Flatten the BH to the HG
    bool lFound = false;
    for (l1tpf::Particle & it2 : *hcals) { 
      if(it2.iPhi() == it.iPhi() && it2.iEta() == it.iEta()) { 
	TLorentzVector lVec = it.tp4();
	lVec += it2.tp4();
	it2.setPtEtaPhiM(lVec.Pt(),lVec.Eta(),lVec.Phi(),lVec.M());
	lFound = true; 
	break;
      }
    }
    if(!lFound) hcals->push_back(it);
  }
  for (const l1tpf::Particle & it : *hcals) lHEt[it.aEta()][it.aPhi()] += it.pt();
  unsigned nh = 0;
  for (const l1tpf::Particle & it : *hcals) {
      hcal_subdet = 0; // it->id().subdet(); // FIXME: missing in Particle class
      hcal_ieta   = it.iEta();
      hcal_iphi   = it.iPhi();
      hcal_TwrR   = 0; // towerR; // FIXME: missing in Particle class
      hcal_num    = nh;
      hcal_et     = it.pt();
      simpleCluster(pClust,it.aEta(),it.aPhi(),lHEt); 
      // FIXME GP: shouldn't there be a iEta < 31 as for ECAL, since lHEt is not fixed outside it?   
      //PH : If I recall by construction it should be limited, but let me follow up. 
      // FIXME Gio: shouldn't we continue to the next TP here if we don't find a cluster?
      // PH : Unfortunately the code is written such that ecal is dealt with here, so a continue will ignore H/E = 0 clusters. I agree its confusing.
      hcal_clust_et=0,hcal_clust_eta=0,hcal_clust_phi=0;
      if(pClust.size() > 0) hcal_clust_et =pClust[0].Pt();
      if(pClust.size() > 0) hcal_clust_eta=pClust[0].Eta();
      if(pClust.size() > 0) hcal_clust_phi=pClust[0].Phi();
      pClust.clear();
      std::vector<double> lGenVars;
      genMatch(lGenVars,1,double(hcal_clust_eta),double(hcal_clust_phi),double(it.pt()),genParticles);
      hcal_genPt=0; hcal_genEta=0; hcal_genPhi=0; hcal_genId=0; hcal_dr=0;
      if(lGenVars.size() > 3) { 
      	hcal_genPt   = float(lGenVars[0]);
      	hcal_genEta  = float(lGenVars[1]);
      	hcal_genPhi  = float(lGenVars[2]);
      	hcal_genId   = float(lGenVars[3]);
	hcal_dr = reco::deltaR( hcal_genEta, hcal_genPhi, hcal_clust_eta,hcal_clust_phi );
      }
      //Linking Ecal and Hcal by Delta 1 in iEta iPhi
      bool isZero = true;  if(it.pt() > 2. * 0.5) isZero = false;
      bool link   = false; if(!isZero && hcal_clust_et > 0) link  = true;
      hcal_ecal_et   = 0; hcal_ecal_etcorr = 0; hcal_ecal_eta = 0; hcal_ecal_phi = 0; 
      for(unsigned int i1 = 0; i1 < lEcals.size(); i1++) { 
	if(!isZero && !link) continue; //Avoid double counting ecal
	if( isZero && (it.iEta() != lEcals[i1].eta || it.iPhi() != lEcals[i1].phi)) continue;
      	if(abs(it.iEta()-lEcals[i1].ieta) > 1) continue;
	int pDPhi = abs(it.iPhi()-lEcals[i1].iphi); if(pDPhi > l1tpf::towerNPhi(it.iEta())-pDPhi) pDPhi =   l1tpf::towerNPhi(it.iEta())-pDPhi;
      	if(abs(pDPhi) > 1) continue;
	if(lEcals[i1].et < 0.2) continue;
      	hcal_ecal_et     = lEcals[i1].et;
      	hcal_ecal_etcorr = lEcals[i1].corr_et;
      	hcal_ecal_eta    = lEcals[i1].eta;
      	hcal_ecal_phi    = lEcals[i1].phi;
      	break;
      }
      hcal_corr_et = 0;
      if(hcal_clust_et > 0) hcal_corr_et   = corrector_ ->correct(double(hcal_clust_et),double(hcal_ecal_et),hcal_ieta,hcal_iphi);
      if(hcal_corr_et  > 0) hcal_corr_et   = corrector2_->correct(double(hcal_corr_et) ,double(hcal_clust_et),double(hcal_ecal_et),hcal_ieta,hcal_iphi);
      if(hcal_corr_et  < 0) hcal_corr_et   = 0;
      if(hcal_corr_et)      hcal_corr_emf  = corrector_->ecalFrac();
      if (hcal_corr_et >= 1 || hcal_ecal_etcorr >= 1) {
          l1tpf::Particle calo = connector_->makeCalo(double(hcal_corr_et) ,double(hcal_ecal_etcorr),double(hcal_clust_eta),double(hcal_clust_phi),double(hcal_ecal_eta),double(hcal_ecal_phi));
          connector_->addCalo(calo);
          l1regions_.addCalo(calo); 
      }
      if (hcal_clust_et >= 1 || hcal_ecal_et >= 1) {
          rawconnector_->addCalo(connector_->makeCalo(double(hcal_clust_et),double(hcal_ecal_et),    double(hcal_clust_eta),double(hcal_clust_phi),double(hcal_ecal_eta),double(hcal_ecal_phi)));
      }
      if(hcal_genPt > 1. && zeroSuppress_) fHcalInfoTree->Fill();      
      if(!zeroSuppress_ && hcal_et > 1.)   fHcalInfoTree->Fill();
      nh++;
      if (nh > 99999) break;
  }
  std::vector<combiner::Particle> lRawCalo      = rawconnector_->candidates();
  std::vector<combiner::Particle> lCorrCalo     = connector_->candidates();
  connector_->link(metRate_);
  std::vector<combiner::Particle> lCands        = connector_->candidates();
  std::vector<combiner::Particle> lTKCands      = connector_->tkcandidates();
  connector_->doVertexing();
  std::vector<combiner::Particle> lTKVtxCands   = connector_->tkvtxcandidates();
  connector_->fetchPuppi();
  connector_->fill();
  std::vector<combiner::Particle> lPupCands     = connector_->puppiFetch();

  // Now we run the discretized version
  // First, get a copy of the discretized inputs (for reference later)
  std::vector<combiner::Particle> ll1CaloCands = l1regions_.fetchCalo();
  std::vector<combiner::Particle> ll1TkCands   = l1regions_.fetchTracks();
  // then get global inputs
  float z0 = connector_->dZ();
  std::pair<float,float> alphaC = connector_->alphaCMedRms(), alphaF = connector_->alphaFMedRms();
  if (fDebug) std::cout << " z0 = " << z0 << "\t alphaC = " << alphaC.first << " +/- " << alphaC.second << "\t alphaF = " << alphaF.first << " +/- " << alphaF.second << std::endl;

  // run PF in each region
  for (auto & l1region : l1regions_.regions()) {
      l1pfalgo_.runPF(l1region);
      l1pfalgo_.runPuppi(l1region, z0, -1., alphaC.first, alphaC.second, alphaF.first, alphaF.second);
  }
  /*
  unsigned int lBase = lCands.size(); 
  bool *lFound = new bool[lBase]; for(unsigned int i0 = 0; i0 < lBase; i0++) lFound[i0] = false;  
  for (reco::GenParticleCollection::const_iterator itGenP = genParticles.begin(); itGenP!=genParticles.end(); ++itGenP) {
    if(itGenP->status() != 1 || itGenP->pt() < 5) continue;
    bool pFound = false;
    for(unsigned int i0   = 0; i0 < lBase; i0++) { 
      double pDPhi = fabs(lCands[i0].phi()-itGenP->phi()); if(pDPhi > 2.*TMath::Pi()-pDPhi) pDPhi = 2.*TMath::Pi()-pDPhi;
      double pDEta = fabs(lCands[i0].eta()-itGenP->eta());
      if(sqrt(pDEta*pDEta+pDPhi*pDPhi) > 0.1) continue;
      lFound[i0] = true;
      pFound = true;
    }
    //if(!pFound) std::cout << "Not Found===> " << itGenP->pt() << " -- " << itGenP->eta() << " -- " << itGenP->phi() << " -- " << itGenP->mass() << " ---> " << itGenP->pdgId() << std::endl;
    if(!pFound) { 
      l1tpf::Particle lPart(itGenP->pt(),itGenP->eta(),itGenP->phi(),itGenP->mass(),itGenP->pdgId(),1,0,itGenP->eta(),itGenP->phi());
      //lCands.push_back(lPart);
    }
  }
  */
  addPF(lRawCalo ,"RawCalo" ,iEvent);
  addPF(lCorrCalo,"Calo"    ,iEvent);
  addPF(lTKCands ,"TK"      ,iEvent);
  addPF(lCands   ,"PF"      ,iEvent);
  addPF(lPupCands,"Puppi"   ,iEvent);

  std::vector<combiner::Particle> ll1PFCands   = l1regions_.fetch(false);
  std::vector<combiner::Particle> ll1PupCands  = l1regions_.fetch(true);
  addPF(ll1CaloCands,"L1Calo" ,iEvent);
  addPF(ll1TkCands,  "L1TK"   ,iEvent);
  addPF(ll1PFCands,  "L1PF"   ,iEvent);
  addPF(ll1PupCands, "L1Puppi",iEvent);
  
  metanalyzer_->clear();
  metanalyzer_->setZ(lTKCands,connector_->dZ());
  metanalyzer_->setMETRecoil(2,lTKCands ,false);
  metanalyzer_->setMETRecoil(0,lCands,false);
  metanalyzer_->setMETRecoil(3,lRawCalo ,true);
  metanalyzer_->setMETRecoil(1,lCorrCalo,true);
  metanalyzer_->setMETRecoil(5,lPupCands,false);
  metanalyzer_->setMETRecoil(4,lTKVtxCands ,false);
  metanalyzer_->setGenMET(genParticles);
  metanalyzer_->fill();

  jetanalyzer_->clear();
  jetanalyzer_->setZ(lTKCands,connector_->dZ());
  jetanalyzer_->setGenJets(genParticles,1);
  jetanalyzer_->setJets(lCands,0);
  jetanalyzer_->setJets(lTKCands,2);
  jetanalyzer_->setJets(lTKVtxCands,3);
  jetanalyzer_->setJets(lRawCalo ,4);
  jetanalyzer_->setJets(lCorrCalo,5);
  jetanalyzer_->setJets(lPupCands,6);
  jetanalyzer_->fill();
}
void NtupleProducer::addPF(std::vector<combiner::Particle> &iCandidates,std::string iLabel,edm::Event& iEvent) { 
  corrCandidates_.reset( new PFOutputCollection );
  for(unsigned int i0 = 0; i0 < iCandidates.size(); i0++) { 
    reco::PFCandidate::ParticleType id = reco::PFCandidate::ParticleType::X; 
    int pCharge=0; 
    if(iCandidates[i0].pdgId() == combiner::CH) id = reco::PFCandidate::ParticleType::h;
    if(iCandidates[i0].pdgId() == combiner::EL) id = reco::PFCandidate::ParticleType::e;
    if(iCandidates[i0].pdgId() == combiner::NH) id = reco::PFCandidate::ParticleType::h0;
    if(iCandidates[i0].pdgId() == combiner::GAMMA) id = reco::PFCandidate::ParticleType::gamma;
    if(iCandidates[i0].pdgId() == combiner::MU) id = reco::PFCandidate::ParticleType::mu;
    if(iCandidates[i0].pdgId() == combiner::CH || iCandidates[i0].pdgId() == combiner::EL)  pCharge = 1;
    if(iCandidates[i0].pdgId() == combiner::MU)  pCharge = iCandidates[i0].charge();
    reco::PFCandidate pCand(pCharge,iCandidates[i0].p4(),id);
    corrCandidates_->push_back(pCand);
  }
  //Fill!
  iEvent.put(std::move(corrCandidates_),iLabel);
}
//////////////////////////// ------------------------------------------------------
//////////////////////////// ------------------------------------------------------

//--- Gen Matching
void NtupleProducer::genMatch(std::vector<double> &iGenVars,int iType,double iEta,double iPhi,double iPt,const reco::GenParticleCollection &iGenParticles) { 
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
//3x3 clusterizer with potential to grow clusters commented out
TLorentzVector NtupleProducer::getVector(double iPt[][73],int iEta,int iPhi,int iEta0,int iPhi0,double iNSigma,double iENoise) { //iEtaC,iPhiC
  double lPtTot = 0;
  std::vector<std::pair<int,int> > lGrow;
  for(int i0=-1; i0 < 2; i0++) {
    for(int i1=-1; i1 < 2; i1++) {
      if(i0 == 0 && i1 == 0) continue;
      if(iEta+i0 < 0 || iEta+i0 > l1tpf::towerNEta()) continue; 
      int pPhi = iPhi+i1;
      if(pPhi < 0)  pPhi=l1tpf::towerNPhi(iEta+i0)+pPhi;
      if(pPhi > l1tpf::towerNPhi(iEta+i0)) pPhi=pPhi-l1tpf::towerNPhi(iEta+i0);
      if(( i1+i0 > -1 && i1 > -1) && iPt[iEta][iPhi] <  iPt[iEta+i0][pPhi]) lPtTot += iPt[iEta+i0][pPhi];
      if(!(i1+i0 > -1 && i1 > -1) && iPt[iEta][iPhi] <= iPt[iEta+i0][pPhi]) lPtTot += iPt[iEta+i0][pPhi];
      /*
      if(( i1+i0 > -1 && i1 > -1) && iPt[iEta][iPhi] <  iPt[iEta+i0][iPhi+i1]) continue;
      if(!(i1+i0 > -1 && i1 > -1) && iPt[iEta][iPhi] <= iPt[iEta+i0][iPhi+i1]) continue;
      double lDEta0 = fabs(iEtaC-iEta);
      double lDPhi0 = fabs(iPhiC-iPhi);
      double lDEta1 = fabs(iEtaC-iEta-i0);
      double lDPhi1 = fabs(iPhiC-iPhi-i1);
      if(sqrt(lDEta0*lEta0+lDPhi0*lDPhi0) > sqrt(lDEta1*lEta1+lDPhi1*lDPhi1)) continue;
      if(sqrt(lDEta1*lEta1+lDPhi1*lDPhi1) > fDistance) continue;
      lGrow.push_back(std::pair<int,int>(iEta+i0,iPhi+i1));
      */
    }
  }      
  double lPt = iPt[iEta][iPhi];
  if(lPtTot >  iNSigma * iENoise && !(iEta0 == iEta && iPhi0 == iPhi)) lPt = iPt[iEta0][iPhi0]/lPtTot * lPt;
  if(lPtTot >  iNSigma * iENoise &&  iEta0 == iEta && iPhi0 == iPhi  ) lPt = 0;
  TLorentzVector lVec;
  lVec.SetPtEtaPhiM(lPt,l1tpf::towerEta(l1tpf::translateAEta(iEta,true)),l1tpf::towerPhi(l1tpf::translateAEta(iEta,true),l1tpf::translateAPhi(iPhi,true)),0);
  //if(lPt > 0) for(unsigned int i0 = 0; i0 < lGrow.size(); i0++) { lVec += getVector(iPt,lGrow[i0].first,lGrow[i0].second,iEta,iPhi,iEta,iPhiC);}
  return lVec;
}

//--- Simple Clustering Start with Local maxima (in 3x3) keep adding neighbors 2sigma above threshold require bottom left to be equal or greater (to avoid points with the same value)
void NtupleProducer::simpleCluster(std::vector<TLorentzVector> &iClusters,double  iEta,double iPhi,double iPt[][73],double iNSigma,double iENoise) { 
  for(int i0 = 0; i0 < l1tpf::towerNEta()-1; i0++) { 
    for(int i1 = 0; i1 < l1tpf::towerNPhi(iEta); i1++) { 
      if(i0 != iEta || i1 != iPhi) continue;
      if (iPt[i0][i1] < iNSigma * iENoise)        continue;
      //Max requirement
      TLorentzVector pVec = getVector(iPt,i0,i1,i0,i1,iNSigma,iENoise);
      if(pVec.Pt() == 0) continue;
      for(int i2=-1; i2 < 2; i2++) {
	for(int i3=-1; i3 < 2; i3++) {
	  if(i2 == 0 && i3 == 0) continue;
	  if(i0+i2 < 0 || i0+i2 > l1tpf::towerNEta()-1) continue; 
	  int pPhi = i1+i3;
	  if(pPhi < 0)  pPhi=l1tpf::towerNPhi(iEta)+pPhi;
	  if(pPhi > l1tpf::towerNPhi(i0+i2)) pPhi=pPhi-l1tpf::towerNPhi(i0+i2);
	  pVec += getVector(iPt,i0+i2,pPhi,i0,i1,iNSigma,iENoise);
	}
      }
      iClusters.push_back(pVec);
    }
  }
}
// ------------ method called once each job just before starting event loop  ------------
void 
NtupleProducer::beginJob()
{
  //
  // Create output file, trees, and histograms
  //
  fOutputFile = new TFile(fOutputName.c_str(), "RECREATE");
  fTotalEvents = new TH1D("TotalEvents","TotalEvents",1,-10,10);
  fTrkInfoTree     = new TTree("TrkInfo",   "TrkInfo");
  fEcalInfoTree     = new TTree("EcalInfo",   "EcalInfo");
  fHcalInfoTree     = new TTree("HcalInfo",   "HcalInfo");

  fTrkInfoTree->Branch("runNum",  &runNum,  "runNum/F");
  fTrkInfoTree->Branch("lumiSec", &lumiSec, "lumiSec/F");
  fTrkInfoTree->Branch("evtNum",  &evtNum,  "evtNum/F");
  fTrkInfoTree->Branch("trkNum",  &trkNum,  "trkNum/F");
  fTrkInfoTree->Branch("trkPx",   &trkPx,   "trkPx/F");
  fTrkInfoTree->Branch("trkPy",   &trkPy,   "trkPy/F");
  fTrkInfoTree->Branch("trkPz",   &trkPz,   "trkPz/F");
  fTrkInfoTree->Branch("trkPt",   &trkPt,   "trkPt/F");
  fTrkInfoTree->Branch("trkEta",  &trkEta,  "trkEta/F");
  fTrkInfoTree->Branch("trkPhi",  &trkPhi,  "trkPhi/F");
  fTrkInfoTree->Branch("trkz0",   &trkz0,   "trkz0/F");
  fTrkInfoTree->Branch("trkd0",   &trkd0,   "trkd0/F");
  fTrkInfoTree->Branch("trkEcalPhi",  &trkEcalPhi, "trkEcalPhi/F");
  fTrkInfoTree->Branch("trkEcalEta",  &trkEcalEta, "trkEcalEta/F");
  fTrkInfoTree->Branch("trkEcalR",    &trkEcalR,   "trkEcalR/F");
  fTrkInfoTree->Branch("genPt",       &genPt,   "genPt/F");
  fTrkInfoTree->Branch("genEta",      &genEta,  "genEta/F");
  fTrkInfoTree->Branch("genPhi",      &genPhi,  "genPhi/F");
  fTrkInfoTree->Branch("genid",       &genId,   "genid/F");
  
  fEcalInfoTree->Branch("ecal_subdet", &ecal_subdet, "ecal_subdet/F");
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


  fHcalInfoTree->Branch("hcal_subdet", &hcal_subdet, "hcal_subdet/F");
  fHcalInfoTree->Branch("hcal_ieta", &hcal_ieta, "hcal_ieta/F");
  fHcalInfoTree->Branch("hcal_iphi", &hcal_iphi, "hcal_iphi/F");
  fHcalInfoTree->Branch("hcal_TwrR", &hcal_TwrR, "hcal_TwrR/F");
  fHcalInfoTree->Branch("hcal_num", &hcal_num, "hcal_num/F");
  fHcalInfoTree->Branch("hcal_et", &hcal_et, "hcal_et/F");
  fHcalInfoTree->Branch("hcal_clust_et",  &hcal_clust_et , "hcal_clust_et/F");
  fHcalInfoTree->Branch("hcal_clust_eta", &hcal_clust_eta, "hcal_clust_eta/F");
  fHcalInfoTree->Branch("hcal_clust_phi", &hcal_clust_phi, "hcal_clust_phi/F");
  fHcalInfoTree->Branch("hcal_corr_et"  , &hcal_corr_et  , "hcal_corr_et/F");
  fHcalInfoTree->Branch("hcal_corr_emf" , &hcal_corr_emf , "hcal_corr_emf/F");
  fHcalInfoTree->Branch("ecal_et",      &hcal_ecal_et,     "hcal_ecal_et/F");
  fHcalInfoTree->Branch("ecal_etcorr",  &hcal_ecal_etcorr, "hcal_ecal_etcorr/F");
  fHcalInfoTree->Branch("ecal_eta", &hcal_ecal_eta, "hcal_ecal_eta/F");
  fHcalInfoTree->Branch("ecal_phi", &hcal_ecal_phi, "hcal_ecal_phi/F");
  fHcalInfoTree->Branch("genPt",   &hcal_genPt ,"hcal_genPt/F");
  fHcalInfoTree->Branch("genEta",  &hcal_genEta,"hcal_genEta/F");
  fHcalInfoTree->Branch("genPhi",  &hcal_genPhi,"hcal_genPhi/F");
  fHcalInfoTree->Branch("genid",   &hcal_genId, "hcal_genid/F");
  fHcalInfoTree->Branch("gendr",   &hcal_dr ,"hcal_dr/F");  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
NtupleProducer::endJob() {
  //
  // Save to ROOT file
  //
  
  connector_->write();
  metanalyzer_->write();
  jetanalyzer_->write();
  fOutputFile->cd();
  fTotalEvents->Write();
  fOutputFile->Write();
  fOutputFile->Close();
}

// ------------ method called when starting to processes a run  ------------

void
NtupleProducer::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
{
}

 
// ------------ method called when ending the processing of a run  ------------
/*
void
NtupleProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
NtupleProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
NtupleProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NtupleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NtupleProducer);
