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
#include "CondFormats/L1TObjects/interface/L1CaloEcalScale.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
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

//
// class declaration
//

class NtupleProducer : public edm::EDProducer {
public:
  explicit NtupleProducer(const edm::ParameterSet&);
  ~NtupleProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  const edm::InputTag L1TrackTag_;
  const edm::InputTag EcalTPTag_;
  const edm::InputTag HcalTPTag_;
  const L1CaloEcalScale* ecalScale_;
  const L1CaloHcalScale* hcalScale_;
private:
  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void propagate(int iOption,std::vector<double> &iVars,const XYZTLorentzVector& iMom,const XYZTLorentzVector& iVtx,double iCharge,double iBField);
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;

  std::pair<float,float> towerEtaBounds(int ieta);
  float towerEta(int ieta);
  float towerPhi(int ieta, int iphi);
  float towerEtaSize(int ieta);
  float towerPhiSize(int ieta);

  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  typedef std::vector<double> HCalCollection;
  typedef std::vector<double> ECalCollection;
  typedef std::vector<double> TrkCollection;
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
  L1TrackTag_(iConfig.getParameter<edm::InputTag>("L1TrackTag")),
  EcalTPTag_(iConfig.getParameter<edm::InputTag>("EcalTPTag")),
  HcalTPTag_(iConfig.getParameter<edm::InputTag>("HcalTPTag"))
{
  //now do what ever other initialization is needed
  // ECAL and HCAL LSB = 0.5
  ecalScale_ = new L1CaloEcalScale(0.5);
  hcalScale_ = new L1CaloHcalScale(0.5);
  produces<HCalCollection>( "hcalET" ).setBranchAlias( "hcalET" );
  produces<ECalCollection>( "ecalET" ).setBranchAlias( "ecalET" ); 
  produces<TrkCollection>( "trkPtPerp" ).setBranchAlias( "trkPtPerp" ); 
  produces<TrkCollection>( "trkPtEta" ).setBranchAlias( "trkPtEta" ); 
  produces<TrkCollection>( "trkPtPhi" ).setBranchAlias( "trkPtPhi" ); 
  produces<TrkCollection>( "trkPOCAz" ).setBranchAlias( "trkPOCAz" );  
  produces<TrkCollection>( "trkEcalEta" ).setBranchAlias( "trkEcalEta" );  
  produces<TrkCollection>( "trkEcalPhi" ).setBranchAlias( "trkEcalPhi" );  
  produces<TrkCollection>( "trkEcalR"   ).setBranchAlias( "trkEcalR" );  
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

   using namespace edm;
   typedef std::vector<TTTrack<Ref_PixelDigi_> >         vec_track;

   iSetup.get<CaloGeometryRecord>().get(calo);
   const CaloGeometry* geo = calo.product(); 

   iSetup.get<CaloGeometryRecord>().get(theTrigTowerGeometry);
   const HcalTrigTowerGeometry* geoTrig = theTrigTowerGeometry.product();


   /// ----------------TRACK INFO-------------------
   /// Stealing Jia Fu's code!
   /// https://github.com/jiafulow/SLHCL1TrackTriggerSimulations/blob/master/NTupleTools/src/NTupleTTTracks.cc
   edm::Handle< vec_track > pixelDigiTTTracks;
   iEvent.getByLabel(L1TrackTag_, pixelDigiTTTracks);

   std::auto_ptr<TrkCollection> trkPtPerp ( new TrkCollection );
   std::auto_ptr<TrkCollection> trkPtEta ( new TrkCollection );
   std::auto_ptr<TrkCollection> trkPtPhi ( new TrkCollection );
   std::auto_ptr<TrkCollection> trkPOCAz( new TrkCollection );
   std::auto_ptr<TrkCollection> trkEcalEta ( new TrkCollection );
   std::auto_ptr<TrkCollection> trkEcalPhi ( new TrkCollection );
   std::auto_ptr<TrkCollection> trkEcalR   ( new TrkCollection );

   const int trksize = pixelDigiTTTracks->size();      
   trkPtPerp->reserve( trksize );
   trkPtEta->reserve( trksize );
   trkPtPhi->reserve( trksize );
   trkPOCAz->reserve( trksize );
   trkEcalEta->reserve(trksize);
   trkEcalPhi->reserve(trksize);
   trkEcalR  ->reserve(trksize);
   
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

       trkPtPerp->push_back(momentum.perp());
       trkPtEta->push_back(momentum.eta());
       trkPtPhi->push_back(momentum.phi());
       trkPOCAz->push_back(poca.z());
       std::vector<double> lVars;
       //Check below
       double charge = it->getRInv()/fabs(it->getRInv());
       propagate(1,lVars,tMom,tVtx,charge,fBz);
       trkEcalEta->push_back(lVars[4]);
       trkEcalPhi->push_back(lVars[5]);
       trkEcalR  ->push_back(lVars[6]);
       //std::cout << "=== Ecal surface eta : " << lVars[4] << " -- phi : " << lVars[5] << " -- R: " << lVars[6] << " -- eta/phi : " << momentum.eta() << " -- " << momentum.phi() << " -- pt : " << momentum.perp() << " -- " << sqrt(lVars[0]*lVars[0]+lVars[1]*lVars[1])   << std::endl;
       // Std::cout << "track info = " << momentum.perp() << "," << momentum.eta() << "," << momentum.phi() << "; poca z = " << poca.z() << std::endl;        
       n++;
       if (n > 10) break; // just for testing so cut it off for now
     }  
   }
   
   /// ----------------ECAL INFO-------------------
   /// Stealing Jia Fu's code!
   /// https://github.com/jiafulow/SLHCL1TrackTriggerSimulations/blob/master/NTupleTools/src/NTupleTTTracks.cc
   /// 72 in phi and 56 in eta
   /// Ecal TPs
   edm::Handle< EcalTrigPrimDigiCollection > ecalTPs;
   iEvent.getByLabel(EcalTPTag_, ecalTPs);
   
   std::auto_ptr<ECalCollection> ecalET ( new ECalCollection );
   const int size = ecalTPs->size();
   ecalET->reserve( size );
   std::cout << "ecalTPs size =  " << ecalTPs->size() << std::endl;

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

       std::cout << "ecal info: " << it->id().subDet() << "," << it->id().ieta() << "," << it->id().iphi() << "," << et << "," << curTowerPhi << "," << curTowerEta << std::endl;

       ecalET->push_back( et );
       // if (et > 0) std::cout << "ecal info: " << it->id().ieta() << "," << it->id().iphi() << "," << et << std::endl;
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
    std::cout << "hcalTPs size =  " << hcalTPs->size() << std::endl;

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

        std::cout << "hcal info: " << it->id().subdet() << "," << it->id().ieta() << "," << it->id().iphi() << "," << et << "," << towerPhi << "," << towerEta << "," << towerR << std::endl;  

        nh++;
        if (nh > 99999) break;

      }
    }

   iEvent.put ( ecalET, "ecalET" );
   iEvent.put ( trkPtPerp, "trkPtPerp" );
   iEvent.put ( trkPtEta, "trkPtEta" );
   iEvent.put ( trkPtPhi, "trkPtPhi" );
   iEvent.put ( trkPOCAz, "trkPOCAz" );
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

// ------------ method called once each job just before starting event loop  ------------
void 
NtupleProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
NtupleProducer::endJob() {
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
