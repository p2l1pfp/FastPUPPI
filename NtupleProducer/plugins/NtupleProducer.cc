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
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  typedef std::vector<double> HCalCollection;
  typedef std::vector<double> ECalCollection;
  typedef std::vector<double> TrkCollection;
  double fBz;
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
   /// Ecal TPs
   edm::Handle< EcalTrigPrimDigiCollection > ecalTPs;
   iEvent.getByLabel(EcalTPTag_, ecalTPs);
   
   std::auto_ptr<ECalCollection> ecalET ( new ECalCollection );
   const int size = ecalTPs->size();
   ecalET->reserve( size );

   if (ecalTPs.isValid()){
     unsigned ne = 0;
     for (EcalTrigPrimDigiCollection::const_iterator it = ecalTPs->begin(); it != ecalTPs->end(); ++it) {

       short ieta = (short) it->id().ieta(); 
       unsigned short absIeta = (unsigned short) abs(ieta);
       short sign = ieta/absIeta;
      
       // unsigned short cal_iphi = (unsigned short) it->id().iphi(); 
       // unsigned short iphi = (72 + 18 - cal_iphi) % 72; // transform TOWERS (not regions) into local rct (intuitive) phi bins
      
       unsigned short compEt = it->compressedEt();
       double et = 0.;
       if (ecalScale_!=0) et = ecalScale_->et( compEt, absIeta, sign );

       ecalET->push_back( et );
       // if (et > 0) std::cout << "ecal info: " << it->id().ieta() << "," << it->id().iphi() << "," << et << std::endl;
       ne++;
       // if (ne > 20) break;
     }
   }
   iEvent.put ( ecalET, "ecalET" );
   iEvent.put ( trkPtPerp, "trkPtPerp" );
   iEvent.put ( trkPtEta, "trkPtEta" );
   iEvent.put ( trkPtPhi, "trkPtPhi" );
   iEvent.put ( trkPOCAz, "trkPOCAz" );
}
// -- propa gator 
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
