// -*- C++ -*-
//
// Package:    InputBuilder
// Class:      InputBuilder
// 
/**\class InputBuilder InputBuilder.cc FastPUPPI/InputBuilder/plugins/InputBuilder.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nhan Viet Tran
//         Created:  Fri, 24 Jun 2016 07:55:13 GMT
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#ifdef HASL1TK
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"     
#endif
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CondFormats/L1TObjects/interface/L1CaloEcalScale.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
// class declaration
//

using namespace std;


class InputBuilder : public edm::EDAnalyzer {
   public:
      explicit InputBuilder(const edm::ParameterSet&);
      ~InputBuilder();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      const edm::InputTag L1TrackTag_;
      const edm::InputTag EcalTPTag_;
      const edm::InputTag HcalTPTag_;
      const L1CaloEcalScale* ecalScale_;
      const L1CaloHcalScale* hcalScale_;

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
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
InputBuilder::InputBuilder(const edm::ParameterSet& iConfig):
  L1TrackTag_(iConfig.getParameter<edm::InputTag>("L1TrackTag")),
  EcalTPTag_(iConfig.getParameter<edm::InputTag>("EcalTPTag")),
  HcalTPTag_(iConfig.getParameter<edm::InputTag>("HcalTPTag"))
{
  //now do what ever initialization is needed
  // ECAL and HCAL LSB = 0.5
  ecalScale_ = new L1CaloEcalScale(0.5);
  hcalScale_ = new L1CaloHcalScale(0.5);

}


InputBuilder::~InputBuilder()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
InputBuilder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
#ifdef HASL1TK
  typedef std::vector<TTTrack<Ref_PixelDigi_> >         vec_track;

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

      std::cout << "track info = " << momentum.perp() << "," << momentum.eta() << "," << momentum.phi() << "; poca z = " << poca.z() << std::endl;        

      n++;
      if (n > 10) break; // just for testing so cut it off for now
    }  
  }
#endif
  /// ----------------ECAL INFO-------------------
  /// Stealing Jia Fu's code!
  /// https://github.com/jiafulow/SLHCL1TrackTriggerSimulations/blob/master/NTupleTools/src/NTupleTTTracks.cc
  /// Ecal TPs
  edm::Handle< EcalTrigPrimDigiCollection > ecalTPs;
  iEvent.getByLabel(EcalTPTag_, ecalTPs);
  std::cout << "ecalTPs size =  " << ecalTPs->size() << std::endl;
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

      if (et > 0) std::cout << "ecal info: " << it->id().ieta() << "," << it->id().iphi() << "," << et << std::endl;

      ne++;
      // if (ne > 20) break;
    }
  }

  /// ----------------HCAL INFO-------------------
  /// Stealing some other code!
  /// https://github.com/cms-sw/cmssw/blob/0397259dd747cee94b68928f17976224c037057a/L1Trigger/L1TNtuples/src/L1AnalysisCaloTP.cc#L40
  /// Hcal TPs
  edm::Handle< HcalTrigPrimDigiCollection > hcalTPs;
  iEvent.getByLabel(HcalTPTag_, hcalTPs);
  std::cout << "hcalTPs size =  " << hcalTPs->size() << std::endl;
  if (hcalTPs.isValid()){
    
    unsigned nh = 0;

    for (HcalTrigPrimDigiCollection::const_iterator it = hcalTPs->begin(); it != hcalTPs->end(); ++it) {

      short ieta = (short) it->id().ieta(); 
      unsigned short absIeta = (unsigned short) abs(ieta);
      short sign = ieta/absIeta;

      unsigned short cal_iphi = (unsigned short) it->id().iphi();
      unsigned short iphi = (72 + 18 - cal_iphi) % 72;
      if (absIeta >= 29) {  // special treatment for HF
        iphi = iphi/4;
      }

      unsigned short compEt = it->SOI_compressedEt();
      double et = 0.;
      if (hcalScale_!=0) et = hcalScale_->et( compEt, absIeta, sign );

      std::cout << "hcal info: " << it->id().ieta() << "," << it->id().iphi() << "," << et << std::endl;

      nh++;
      if (nh > 20) break;

    }
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
InputBuilder::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
InputBuilder::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
InputBuilder::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
InputBuilder::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
InputBuilder::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
InputBuilder::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
InputBuilder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(InputBuilder);
