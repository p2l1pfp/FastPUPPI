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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"     
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/PtEtaPhiMass.h"
#include "CondFormats/L1TObjects/interface/L1CaloEcalScale.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "FastPUPPI/NtupleProducer/interface/corrector.hh"
#include "FastPUPPI/NtupleProducer/interface/combiner.hh"

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

//--------------------------------------------------------------------------------------------------
class NtupleProducer : public edm::EDProducer {
public:
  explicit NtupleProducer(const edm::ParameterSet&);
  ~NtupleProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  const bool          zeroSuppress_;
  const edm::InputTag L1TrackTag_;
  const edm::InputTag EcalTPTag_;
  const edm::InputTag HcalTPTag_;
  const edm::InputTag GenParTag_;
  const edm::InputTag CorrectorTag_;
  const edm::InputTag ECorrectorTag_;
  const edm::InputTag TrackResTag_;
  const edm::InputTag EleResTag_;
  const edm::InputTag PionResTag_;
  const L1CaloEcalScale* ecalScale_;
  const L1CaloHcalScale* hcalScale_;
  std::auto_ptr<PFOutputCollection > corrCandidates_;
private:
  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void propagate(int iOption,std::vector<double> &iVars,const XYZTLorentzVector& iMom,const XYZTLorentzVector& iVtx,double iCharge,double iBField);
  void genMatch(std::vector<double> &iGenVars,int iType,double iEta,double iPhi,double iPt,const reco::GenParticleCollection &iGenParticles);
  TLorentzVector getVector(double iPt[][72],int iEta,int iPhi,int iEta0,int iPhi0,double iNSigma=2,double iENoise=1);
  void simpleCluster(std::vector<TLorentzVector> &iClusters,double  iEta,double iPhi,double iPt[][72],double iNSigma=2,double iENoise=1);
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  
  std::pair<float,float> towerEtaBounds(int ieta);
  float towerEta(int ieta);
  float towerPhi(int ieta, int iphi);
  float towerEtaSize(int ieta);
  float towerPhiSize(int ieta);
  int   translateIEta(int ieta,bool iInvert=false);
  int   translateIPhi(int iphi,bool iInvert=false);

  corrector* corrector_;
  corrector* ecorrector_;
  combiner * connector_;
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

  float ecal_subdet, ecal_ieta, ecal_iphi, ecal_curTwrEta, ecal_curTwrPhi, ecal_et, ecal_num;
  float hcal_subdet, hcal_ieta, hcal_iphi, hcal_TwrR, hcal_et, hcal_num, hcal_ecal_et,hcal_ecal_etcorr,hcal_ecal_eta,hcal_ecal_phi;
  float ecal_clust_et,ecal_clust_eta,ecal_clust_phi,ecal_corr_et;
  float hcal_clust_et,hcal_clust_eta,hcal_clust_phi,hcal_corr_et,hcal_corr_emf;
  float ecal_genPt, ecal_genEta, ecal_genPhi, ecal_genId;
  float hcal_genPt, hcal_genEta, hcal_genPhi, hcal_genId;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  double fBz;

  edm::ESHandle<CaloGeometry> calo;
  edm::ESHandle<HcalTrigTowerGeometry> theTrigTowerGeometry;

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
  HcalTPTag_            (iConfig.getParameter<edm::InputTag>("HcalTPTag")),
  GenParTag_            (iConfig.getParameter<edm::InputTag>("genParTag")),
  CorrectorTag_         (iConfig.getParameter<edm::InputTag>("corrector")),
  ECorrectorTag_        (iConfig.getParameter<edm::InputTag>("ecorrector")),
  TrackResTag_          (iConfig.getParameter<edm::InputTag>("trackres")),
  EleResTag_            (iConfig.getParameter<edm::InputTag>("eleres")),
  PionResTag_           (iConfig.getParameter<edm::InputTag>("pionres")),
  fOutputName           (iConfig.getUntrackedParameter<std::string>("outputName", "ntuple.root")),
  fOutputFile           (0),
  fTotalEvents          (0),
  fTrkInfoTree          (0),
  fEcalInfoTree         (0),
  fHcalInfoTree         (0)
{
  //now do what ever other initialization is needed
  // ECAL and HCAL LSB = 0.5
  ecalScale_  = new L1CaloEcalScale(0.5);
  hcalScale_  = new L1CaloHcalScale(0.5);
  corrector_  = new corrector(CorrectorTag_.label());
  ecorrector_ = new corrector(ECorrectorTag_.label(),1);
  connector_  = new combiner (TrackResTag_.label(),EleResTag_.label(),PionResTag_.label());
  produces<PFOutputCollection>();
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
NtupleProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  fTotalEvents->Fill(1);  
  using namespace edm;
  typedef std::vector<TTTrack<Ref_PixelDigi_> >         vec_track;
  
  iSetup.get<CaloGeometryRecord>().get(calo);
  const CaloGeometry* geo = calo.product(); 
  
  iSetup.get<CaloGeometryRecord>().get(theTrigTowerGeometry);
  const HcalTrigTowerGeometry* geoTrig = theTrigTowerGeometry.product();
  
  edm::Handle<reco::GenParticleCollection> hGenParProduct;
  iEvent.getByLabel(GenParTag_,hGenParProduct);
  const reco::GenParticleCollection genParticles = *(hGenParProduct.product());  
  /// ----------------TRACK INFO-------------------
  /// Stealing Jia Fu's code!
  /// https://github.com/jiafulow/SLHCL1TrackTriggerSimulations/blob/master/NTupleTools/src/NTupleTTTracks.cc
  edm::Handle< vec_track > pixelDigiTTTracks;
  iEvent.getByLabel(L1TrackTag_, pixelDigiTTTracks);
  connector_->clear();
  if (pixelDigiTTTracks.isValid()) {
    
    edm::LogInfo("NTupleTracks") << "Size: " << pixelDigiTTTracks->size();
    
    unsigned nPar = 4;
    unsigned n = 0;
    for (vec_track::const_iterator it = pixelDigiTTTracks->begin(); it != pixelDigiTTTracks->end(); ++it) {


      const GlobalVector&          momentum = it->getMomentum(nPar);
      const GlobalPoint&           poca     = it->getPOCA(nPar);  // point of closest approach
      bool iIsEle=false;
      double mass  = iIsEle ? 0.0005 : 0.139; 
      double energy=sqrt((mass*mass)+(momentum.mag()*momentum.mag()));
      const XYZTLorentzVector      tMom (momentum.x(),momentum.y(),momentum.z(),energy);
      const XYZTLorentzVector      tVtx (poca.x()     ,poca.y()     ,poca.z(),0.);
      
      std::vector<double> lVars;
      //Check below
      double charge = it->getRInv()/fabs(it->getRInv());
      propagate(1,lVars,tMom,tVtx,charge,fBz);
  
      runNum  = iEvent.id().run();
      lumiSec = iEvent.luminosityBlock();
      evtNum  = iEvent.id().event();
      trkNum  = n;
      
      trkEcalEta = float(lVars[4]);
      trkEcalPhi = float(lVars[5]);
      trkEcalR   = float(lVars[6]);
      
      trkPx  = momentum.x();
      trkPz  = momentum.y();
      trkPy  = momentum.z();
      trkPt  = momentum.perp();
      trkEta = momentum.eta();
      trkPhi = momentum.phi();
      trkz0  = poca.z();
      trkd0  = poca.perp();
      connector_->addTrack(trkPt,trkEta,trkPhi,trkz0,trkEcalEta,trkEcalPhi);      

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
      //std::cout << "=== Ecal surface eta : " << lVars[4] << " -- phi : " << lVars[5] << " -- R: " << lVars[6] << " -- eta/phi : " << momentum.eta() << " -- " << momentum.phi() << " -- pt : " << momentum.perp() << " -- " << sqrt(lVars[0]*lVars[0]+lVars[1]*lVars[1])   << std::endl;
      // Std::cout << "track info = " << momentum.perp() << "," << momentum.eta() << "," << momentum.phi() << "; poca z = " << poca.z() << std::endl;        
      n++;
      // if (n > 10) break; // just for testing so cut it off for now
    }  
  }
   
  
  /// ----------------ECAL INFO-------------------
  /// Stealing Jia Fu's code!
  /// https://github.com/jiafulow/SLHCL1TrackTriggerSimulations/blob/master/NTupleTools/src/NTupleTTTracks.cc
  /// 72 in phi and 56 in eta
  /// Ecal TPs
  edm::Handle< EcalTrigPrimDigiCollection > ecalTPs;
  iEvent.getByLabel(EcalTPTag_, ecalTPs);
  
  //std::cout << "ecalTPs size =  " << ecalTPs->size() << std::endl;
  double lEta[4032];
  double lPhi[4032];
  double lEt [4032][5];
  double lEEt[61][72];
  double lHEt[61][72];
  for(int i0 = 0; i0 < 61; i0++) { for(int i1 = 0; i1 < 72; i1++) {lEEt[i0][i1]=0; lHEt[i0][i1] = 0;}}
  std::vector<TLorentzVector> pClust;
  if (ecalTPs.isValid()){
    unsigned ne = 0;
    for (EcalTrigPrimDigiCollection::const_iterator it = ecalTPs->begin(); it != ecalTPs->end(); ++it) {
      int pIPhi = translateIPhi(it->id().iphi());
      int pIEta = translateIEta(it->id().ieta());
      double et = 0.;
      short sign = 1; if(it->id().ieta() < 0) sign=-1;
      if (ecalScale_!=0) et = ecalScale_->et( it->compressedEt(),(unsigned short)abs(it->id().ieta()), sign);
      if(et < 0.1) et = 0.;
      lEEt[pIEta][pIPhi] = et;
    }
    for (EcalTrigPrimDigiCollection::const_iterator it = ecalTPs->begin(); it != ecalTPs->end(); ++it) {
      short ieta = (short) it->id().ieta(); 
      unsigned short absIeta = (unsigned short) abs(ieta);
      short sign = ieta/absIeta;
      
      unsigned short compEt = it->compressedEt();
      double et = 0.;
      if (ecalScale_!=0) et = ecalScale_->et( compEt, absIeta, sign );
      
      float curTowerEta = towerEta(it->id().ieta());
      float curTowerPhi = towerPhi(it->id().ieta(),it->id().iphi());
      
      ecal_subdet = it->id().subDet();
      ecal_ieta = it->id().ieta();
      ecal_iphi = it->id().iphi();
      ecal_et = et;
      ecal_curTwrEta = curTowerEta;
      ecal_curTwrPhi = curTowerPhi;
      ecal_num = ne;
      if(fabs(ecal_ieta) < 31) simpleCluster(pClust,translateIEta(ecal_ieta),translateIPhi(ecal_iphi),lEEt);
      ecal_clust_et=0,ecal_clust_eta=0,ecal_clust_phi =0;
      if(pClust.size() > 0) ecal_clust_et =pClust[0].Pt();
      if(pClust.size() > 0) ecal_clust_eta=pClust[0].Eta();
      if(pClust.size() > 0) ecal_clust_phi=pClust[0].Phi();
      pClust.clear();
      ecal_corr_et = ecorrector_->correct(0.,double(ecal_clust_et),ecal_ieta);
      lEta[ne]  = it->id().ieta();
      lPhi[ne]  = it->id().iphi();
      lEt [ne][0] = ecal_clust_et;
      lEt [ne][1] = ecal_corr_et;
      lEt [ne][2] = ecal_clust_eta;
      lEt [ne][3] = ecal_clust_phi;
      std::vector<double> lGenVars;
      genMatch(lGenVars,1,double(curTowerEta),double(curTowerPhi),double(et),genParticles);
      ecal_genPt=0; ecal_genEta=0; ecal_genPhi=0; ecal_genId=0;
      if(lGenVars.size() > 3) { 
	ecal_genPt   = float(lGenVars[0]);
	ecal_genEta  = float(lGenVars[1]);
	ecal_genPhi  = float(lGenVars[2]);
	ecal_genId   = float(lGenVars[3]);
      }
      if(ecal_genPt > 1. && zeroSuppress_) fEcalInfoTree->Fill();      
      if(!zeroSuppress_) fEcalInfoTree->Fill();
      //std::cout << "ecal info: " << it->id().subDet() << "," << it->id().ieta() << "," << it->id().iphi() << "," << et << "," << curTowerPhi << "," << curTowerEta << std::endl;
      ne++;
      // if (ne > 20) break;
    }
  }
  
  
  // / ----------------HCAL INFO-------------------
  // / Stealing some other code!
  // / https://github.com/cms-sw/cmssw/blob/0397259dd747cee94b68928f17976224c037057a/L1Trigger/L1TNtuples/src/L1AnalysisCaloTP.cc#L40
  // / 72 in phi and 56 in eta = 4032
  // / 144 TP is HF: 4 (in eta) * 18 (in phi) * 2 (sides)
  // / Hcal TPs
  edm::Handle< HcalTrigPrimDigiCollection > hcalTPs;
  iEvent.getByLabel(HcalTPTag_, hcalTPs);
  if (hcalTPs.isValid()){
    unsigned nh = 0;
    for (HcalTrigPrimDigiCollection::const_iterator it = hcalTPs->begin(); it != hcalTPs->end(); ++it) {
      int pIPhi = translateIPhi(it->id().iphi());
      int pIEta = translateIEta(it->id().ieta());
      double et = 0.;      
      if (hcalScale_!=0) et = hcalScale_->et( it->SOI_compressedEt(), it->id().ietaAbs(), it->id().zside() );
      if(abs(it->id().ieta()) < 31) lHEt[pIEta][pIPhi] = et;
    }
    for (HcalTrigPrimDigiCollection::const_iterator it = hcalTPs->begin(); it != hcalTPs->end(); ++it) {
      std::vector<HcalDetId> detids = geoTrig->detIds(it->id());
      unsigned int ndetsPerTower = detids.size();
      float towerEta = 0.;
      float towerPhi = 0.;
      float towerR   = 0.;
      float curphimax = -9999.;
      float curphimin = 9999.;
      for (unsigned int i = 0; i < ndetsPerTower; i++){
	const CaloCellGeometry *cell = geo->getGeometry( detids[i] );       
	// std::cout << detids[i].det() << "," << detids[i].subdetId() << "," << detids[i].iphi() << "," << detids[i].ieta() << "," << cell->etaPos() << "," << cell->phiPos() << std::endl;
	towerEta += cell->etaPos();
	towerR   += cell->getPosition().mag();
	
	if (curphimax < cell->phiPos()) curphimax = cell->phiPos();
	if (curphimin > cell->phiPos()) curphimin = cell->phiPos();
      }
      // std::cout << "curphimax = " << curphimax << ", curphimin = " << curphimin << std::endl;
      if (curphimax - curphimin < PI){ towerPhi = (curphimax + curphimin)/2.; }
      else{ towerPhi = -3.14159;} // special case
      // else{ towerPhi = (curphimax - curphimin + 2*PI)/2.; }
      // std::cout << "towerPhi = " << towerPhi << std::endl;
      // if (fabs(towerPhi) > PI && towerPhi > 0) towerPhi -= PI;
      // if (fabs(towerPhi) > PI && towerPhi < 0) towerPhi += PI;
      towerPhi *= 1.;
      towerEta /= float(ndetsPerTower);
      towerR   /= float(ndetsPerTower);
      hcal_subdet = it->id().subdet();
      hcal_ieta   = it->id().ieta();
      hcal_iphi   = it->id().iphi();
      hcal_TwrR   = towerR;
      hcal_num    = nh;
      unsigned short compEt = it->SOI_compressedEt();
      double et = 0.;
      if (hcalScale_!=0) et = hcalScale_->et( compEt, it->id().ietaAbs(), it->id().zside() );
      hcal_et     = et;
      simpleCluster(pClust,translateIEta(hcal_ieta),translateIPhi(hcal_iphi),lHEt);
      hcal_clust_et=0,hcal_clust_eta=0,hcal_clust_phi=0;
      if(pClust.size() > 0) hcal_clust_et =pClust[0].Pt();
      if(pClust.size() > 0) hcal_clust_eta=pClust[0].Eta();
      if(pClust.size() > 0) hcal_clust_phi=pClust[0].Phi();
      pClust.clear();
     
      //std::cout << "hcal info: " << it->id().subdet() << "," << it->id().ieta() << "," << it->id().iphi() << "," << et << "," << towerPhi << "," << towerEta << "," << towerR << std::endl;  
      std::vector<double> lGenVars;
      genMatch(lGenVars,1,double(hcal_clust_eta),double(hcal_clust_phi),double(et),genParticles);
      hcal_genPt=0; hcal_genEta=0; hcal_genPhi=0; hcal_genId=0;
      if(lGenVars.size() > 3) { 
	hcal_genPt   = float(lGenVars[0]);
	hcal_genEta  = float(lGenVars[1]);
	hcal_genPhi  = float(lGenVars[2]);
	hcal_genId   = float(lGenVars[3]);
      }
      hcal_ecal_et   = -1;
      for(int i1 = 0; i1 < 4032; i1++) { 
	if(it->id().ieta() != lEta[i1]) continue;
	if(it->id().iphi() != lPhi[i1]) continue;
	hcal_ecal_et     = lEt[i1][0];
	hcal_ecal_etcorr = lEt[i1][1];
	hcal_ecal_eta    = lEt[i1][2];
	hcal_ecal_phi    = lEt[i1][3];
	break;
      }
      hcal_corr_et = 0;
      if(hcal_ecal_et > -1 || hcal_clust_et > 0) hcal_corr_et   = corrector_->correct(double(hcal_clust_et),double(hcal_ecal_et),hcal_ieta);
      if(hcal_corr_et) hcal_corr_emf  = corrector_->ecalFrac();
      connector_->addCalo(double(hcal_corr_et),double(hcal_ecal_etcorr),double(hcal_clust_eta),double(hcal_clust_phi),double(hcal_ecal_eta),double(hcal_ecal_phi));
      //if(hcal_genPt > 5.) std::cout << "===> Ecal Clust " << hcal_ecal_et << " -- " << hcal_clust_et << " -- Et " << hcal_et << " ---> " << hcal_genPt << std::endl;
      if(hcal_genPt > 1. && zeroSuppress_) fHcalInfoTree->Fill();      
      if(!zeroSuppress_) fHcalInfoTree->Fill();
      nh++;
      if (nh > 99999) break;
    }
  }
  connector_->link();
  std::vector<combiner::Particle> lCandidates = connector_->candidates();
  corrCandidates_.reset( new PFOutputCollection );
  for(unsigned int i0 = 0; i0 < lCandidates.size(); i0++) { 
    reco::PFCandidate::ParticleType id = reco::PFCandidate::ParticleType::X; 
    int pCharge=0; 
    if(lCandidates[i0].id == 0) id = reco::PFCandidate::ParticleType::h;
    if(lCandidates[i0].id == 1) id = reco::PFCandidate::ParticleType::e;
    if(lCandidates[i0].id == 2) id = reco::PFCandidate::ParticleType::h0;
    if(lCandidates[i0].id == 3) id = reco::PFCandidate::ParticleType::gamma;
    if(lCandidates[i0].id < 2)  pCharge = 1;
    TLorentzVector pVec;  pVec. SetPtEtaPhiM(lCandidates[i0].Et,lCandidates[i0].Eta,lCandidates[i0].Phi,lCandidates[i0].M);
    LorentzVector  pLVec; pLVec.SetPxPyPzE(pVec.Px(),pVec.Py(),pVec.Pz(),pVec.E());
    reco::PFCandidate pCand(pCharge,pLVec,id);
    corrCandidates_->push_back(pCand);
  }
  //Fill!
  iEvent.put(corrCandidates_);
}
// -- propagator 
void NtupleProducer::propagate(int iOption,std::vector<double> &iVars,const XYZTLorentzVector& iMom,const XYZTLorentzVector& iVtx,double iCharge,double iBField) { 
  BaseParticlePropagator particle = BaseParticlePropagator(RawParticle(iMom,iVtx),0.,0.,iBField);
  particle.setCharge(iCharge);
  double ecalShowerDepth=0;
  if(iOption == 0 || iOption == 1) particle.propagateToEcalEntrance(false);
  if(iOption == 2 || iOption == 3) particle.propagateToHcalEntrance(false);
  if(iOption == 4) particle.propagateToHcalExit    (false);
  if(iOption == 1) ecalShowerDepth = reco::PFCluster::getDepthCorrection(particle.momentum().E(),false,false);
  math::XYZVector point = math::XYZVector(particle.vertex())+math::XYZTLorentzVector(particle.momentum()).Vect().Unit()*ecalShowerDepth;
  iVars.push_back(particle.momentum().px());
  iVars.push_back(particle.momentum().py());
  iVars.push_back(particle.momentum().pz());
  iVars.push_back(particle.momentum().energy());
  iVars.push_back(point.eta());
  iVars.push_back(point.phi());
  iVars.push_back(point.rho());
}
// -- helpers for ECAL
std::pair<float,float> NtupleProducer::towerEtaBounds(int ieta)
{
  // if(ieta==0) ieta = 1;
  // if(ieta>kHFEnd) ieta = kHFEnd;
  // if(ieta<(-1*kHFEnd)) ieta = -1*kHFEnd;
  const float towerEtas[33] = {0,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696,0.783,0.870,0.957,1.044,1.131,1.218,1.305,1.392,1.479,1.566,1.653,1.740,1.830,1.930,2.043,2.172,2.322,2.5,2.650,3.000,3.5,4.0,4.5,5.0}; 
  // const float towerEtas[42] = {0,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696,0.783,0.870,0.957,1.044,1.131,1.218,1.305,1.392,1.479,1.566,1.653,1.740,1.830,1.930,2.043,2.172,2.322,2.5,2.650,2.853,3.139,3.314,3.489,3.664,3.839,4.013,4.191,4.363,4.538,4.716,4.889,5.191,5.191};
  return std::make_pair( towerEtas[abs(ieta)-1],towerEtas[abs(ieta)] );
}

float NtupleProducer::towerEta(int ieta)
{
  std::pair<float,float> bounds = towerEtaBounds(ieta);
  float eta = (bounds.second+bounds.first)/2.;
  float sign = ieta>0 ? 1. : -1.;
  return sign*eta; 
}

float NtupleProducer::towerPhi(int ieta, int iphi)
{
  float phi = (float(iphi)-0.5)*towerPhiSize(ieta);
  if (phi > M_PI) phi = phi - (2*M_PI);
  return phi;
}

int NtupleProducer::translateIEta(int ieta,bool iInvert) {
  int lEta = ieta+30;
  if(iInvert) lEta = ieta-30;
  return lEta;
}
int NtupleProducer::translateIPhi(int iphi,bool iInvert) {
  int lPhi = iphi-1;
  if(iInvert) lPhi = iphi+1;
  return lPhi;
}
float NtupleProducer::towerEtaSize(int ieta)
{
  std::pair<float,float> bounds = towerEtaBounds(ieta);
  float size = (bounds.second-bounds.first);
  return size;
}
float NtupleProducer::towerPhiSize(int ieta)
{
  const int kNphi = 72;
  return 2.*M_PI/kNphi;
}
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
    if(deltaR > 0.1) continue;
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
TLorentzVector NtupleProducer::getVector(double iPt[][72],int iEta,int iPhi,int iEta0,int iPhi0,double iNSigma,double iENoise) { //iEtaC,iPhiC
  double lPtTot = 0;
  std::vector<std::pair<int,int> > lGrow;
  for(int i0=-1; i0 < 2; i0++) {
    for(int i1=-1; i1 < 2; i1++) {
      if(i0 == 0 && i1 == 0) continue;
      if(iEta+i0 < 0 || iEta+i0 > 60) continue; 
      int pPhi = iPhi+i1;
      if(pPhi < 0)  pPhi=72+pPhi;
      if(pPhi > 71) pPhi=pPhi-72;
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
  if(lPtTot > 0 && iPt[iEta0][iPhi0]/lPtTot > 1) std::cout << "!!!!!!!!!!!!!====> ptfrac: " << iPt[iEta0][iPhi0]/lPtTot << std::endl;
  if(lPtTot >  iNSigma * iENoise && !(iEta0 == iEta && iPhi0 == iPhi)) lPt = iPt[iEta0][iPhi0]/lPtTot * lPt;
  if(lPtTot >  iNSigma * iENoise &&  iEta0 == iEta && iPhi0 == iPhi  ) lPt = 0;
  TLorentzVector lVec;
  lVec.SetPtEtaPhiM(lPt,towerEta(translateIEta(iEta,true)),towerPhi(translateIEta(iEta,true),translateIPhi(iPhi,true)),0);
  //if(lPt > 0) for(unsigned int i0 = 0; i0 < lGrow.size(); i0++) { lVec += getVector(iPt,lGrow[i0].first,lGrow[i0].second,iEta,iPhi,iEta,iPhiC);}
  return lVec;
}
//--- Simple Clustering Start with Local maxima (in 3x3) keep adding neighbors 2sigma above threshold require bottom left to be equal or greater (to avoid points with the same value)
void NtupleProducer::simpleCluster(std::vector<TLorentzVector> &iClusters,double  iEta,double iPhi,double iPt[][72],double iNSigma,double iENoise) { 
  for(int i0 = 0; i0 < 61; i0++) { 
    for(int i1 = 0; i1 < 72; i1++) { 
      if(i0 != iEta || i1 != iPhi) continue;
      if (iPt[i0][i1] < iNSigma * iENoise)        continue;
      //Max requirement
      TLorentzVector pVec = getVector(iPt,i0,i1,i0,i1,iNSigma,iENoise);
      if(pVec.Pt() == 0) continue;
      for(int i2=-1; i2 < 2; i2++) {
	for(int i3=-1; i3 < 2; i3++) {
	  if(i2 == 0 && i3 == 0) continue;
	  if(i0+i2 < 0 || i0+i2 > 60) continue; 
	  int pPhi = i1+i3;
	  if(pPhi < 0)  pPhi=72+pPhi;
	  if(pPhi > 71) pPhi=pPhi-72;
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
  fEcalInfoTree->Branch("genid",   &ecal_genId, "ecal_enid/F");


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
  fHcalInfoTree->Branch("genid",   &hcal_genId, "hcal_enid/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
NtupleProducer::endJob() {
  //
  // Save to ROOT file
  //
  fOutputFile->cd();
  fTotalEvents->Write();
  fOutputFile->Write();
  fOutputFile->Close();
  
}

// ------------ method called when starting to processes a run  ------------

void
NtupleProducer::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
{
  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  fBz = magneticField->inTesla(GlobalPoint(0,0,0)).z();
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
