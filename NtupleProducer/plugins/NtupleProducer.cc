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
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
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

const double PI = 3.1415926535897;

// ROOT classes
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObject.h>

//--------------------------------------------------------------------------------------------------
class NtupleProducer : public edm::EDProducer {
public:
  explicit NtupleProducer(const edm::ParameterSet&);
  ~NtupleProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  const edm::InputTag L1TrackTag_;
  const edm::InputTag EcalTPTag_;
  const edm::InputTag HcalTPTag_;
  const edm::InputTag GenParTag_;
  const L1CaloEcalScale* ecalScale_;
  const L1CaloHcalScale* hcalScale_;
private:
  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void propagate(int iOption,std::vector<double> &iVars,const XYZTLorentzVector& iMom,const XYZTLorentzVector& iVtx,double iCharge,double iBField);
  void genMatch(std::vector<double> &iGenVars,int iType,double iEta,double iPhi,double iPt,const reco::GenParticleCollection &iGenParticles);
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;

  std::pair<float,float> towerEtaBounds(int ieta);
  float towerEta(int ieta);
  float towerPhi(int ieta, int iphi);
  float towerEtaSize(int ieta);
  float towerPhiSize(int ieta);

  // declare variables for output file
  std::string           fOutputName;
  TFile                 *fOutputFile;
  TH1D                  *fTotalEvents;
  TTree                 *fEventTree;
  TTree                 *fTrkInfoTree;
  
  float runNum, lumiSec, evtNum;
  float trkNum;
  float trkPx, trkPz, trkPy, trkPt, trkEta, trkPhi, trkz0, trkd0;    
  float trkEcalEta, trkEcalPhi, trkEcalR;
  float genPt, genEta, genPhi, genId;

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
  L1TrackTag_           (iConfig.getParameter<edm::InputTag>("L1TrackTag")),
  EcalTPTag_            (iConfig.getParameter<edm::InputTag>("EcalTPTag")),
  HcalTPTag_            (iConfig.getParameter<edm::InputTag>("HcalTPTag")),
  GenParTag_            (iConfig.getParameter<edm::InputTag>("genParTag")),
  fOutputName           (iConfig.getUntrackedParameter<std::string>("outputName", "ntuple.root")),
  fOutputFile           (0),
  fTotalEvents          (0),
  fTrkInfoTree          (0)
{
  //now do what ever other initialization is needed
  // ECAL and HCAL LSB = 0.5
  ecalScale_ = new L1CaloEcalScale(0.5);
  hcalScale_ = new L1CaloHcalScale(0.5);
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
  
  if (ecalTPs.isValid()){
    unsigned ne = 0;
    for (EcalTrigPrimDigiCollection::const_iterator it = ecalTPs->begin(); it != ecalTPs->end(); ++it) {
      
      short ieta = (short) it->id().ieta(); 
      unsigned short absIeta = (unsigned short) abs(ieta);
      short sign = ieta/absIeta;
      
      unsigned short compEt = it->compressedEt();
      double et = 0.;
      if (ecalScale_!=0) et = ecalScale_->et( compEt, absIeta, sign );
      
      float curTowerEta = towerEta(it->id().ieta());
      float curTowerPhi = towerPhi(it->id().ieta(),it->id().iphi());
      
      //std::cout << "ecal info: " << it->id().subDet() << "," << it->id().ieta() << "," << it->id().iphi() << "," << et << "," << curTowerPhi << "," << curTowerEta << std::endl;
      curTowerPhi++;
      curTowerEta++;
      et++;
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
  //std::cout << "hcalTPs size =  " << hcalTPs->size() << std::endl;
  
  if (hcalTPs.isValid()){
    
    unsigned nh = 0;
    
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
      towerEta /= float(ndetsPerTower);
      towerR   /= float(ndetsPerTower);
      unsigned short compEt = it->SOI_compressedEt();
      double et = 0.;
      if (hcalScale_!=0) et = hcalScale_->et( compEt, it->id().ietaAbs(), it->id().zside() );
      
      //std::cout << "hcal info: " << it->id().subdet() << "," << it->id().ieta() << "," << it->id().iphi() << "," << et << "," << towerPhi << "," << towerEta << "," << towerR << std::endl;  
      towerEta++;
      towerPhi++;
      towerR++;
      et++;
      nh++;
      if (nh > 99999) break;
      
    }
  }

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
    if(iType == 1 && itGenP->charge() != 0) continue;
    double deltaEta = itGenP->eta()-iEta;
    double deltaPhi = fabs(itGenP->phi()-iPhi); if(deltaPhi > 2.*TMath::Pi()-deltaPhi) deltaPhi = 2.*TMath::Pi()-deltaPhi;
    double deltaR   = sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);
    if(itGenP->pt() > 20 && iPt > 20) std::cout << "==ieta=> " << iEta << " -iphi- " << iPhi << " -eta- " << itGenP->eta() << " -phi- " << itGenP->phi() << " -ipt- " << iPt << " -pt- " << itGenP->pt() << " -- " << deltaR << " - " << deltaEta << " -- " << deltaPhi << std::endl;
    if(deltaR > 0.1 ) continue;
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
    if(itGenP->pt() > 20 && iPt > 20) std::cout << "==> Matched " << lVec.Pt() << std::endl;
  }
  if(lVec.Pt() > 0) iGenVars.push_back(lVec.Pt());
  if(lVec.Pt() > 0) iGenVars.push_back(lVec.Eta());
  if(lVec.Pt() > 0) iGenVars.push_back(lVec.Phi());
  if(lVec.Pt() > 0) iGenVars.push_back(lId);
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
