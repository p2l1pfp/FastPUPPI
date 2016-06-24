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

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"     
//
// class declaration
//

using namespace std;


class InputBuilder : public edm::EDAnalyzer {
   public:
      explicit InputBuilder(const edm::ParameterSet&);
      ~InputBuilder();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      const edm::InputTag L1TrackTag_;

      // edm::InputTag L1TrackTag_;
      // edm::EDGetTokenT< vector<TTTrack<edm::Ref<edm::DetSetVector<PixelDigi>,PixelDigi,edm::refhelper::FindForDetSetVector<PixelDigi> > > >> > L1TrackTok_;

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
  L1TrackTag_(iConfig.getParameter<edm::InputTag>("L1TrackTag"))
{
   //now do what ever initialization is needed

  // L1TrackTag_ = iConfig.getParameter<edm::InputTag>("JetTag");
  // L1TrackTok_ = consumes< vector<TTTrack<edm::Ref<edm::DetSetVector<PixelDigi>,PixelDigi,edm::refhelper::FindForDetSetVector<PixelDigi> > > >> >(JetTag_);

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

   std::cout << "hello world!" << std::endl;

   // edm::Handle< vector<TTTrack<edm::Ref<edm::DetSetVector<PixelDigi>,PixelDigi,edm::refhelper::FindForDetSetVector<PixelDigi> > > >>  > L1Tracks;
   // iEvent.getByToken(L1TrackTok_,L1Tracks);
   // std::cout << "L1Tracks size = " << L1Tracks->size() << std::endl;

  // L1 tracks & stubs
  // edm::Handle<L1TkTrack_PixelDigi_Collection> L1TrackHandle;
  // edm::InputTag L1TrackTag; 
  // iEvent.getByLabel("L1Tracks", "Level1TkTracks", L1TrackHandle);

    // typedef std::vector<TTTrack<Ref_PixelDigi_> >         vec_track;


  /// TTTrack
  edm::Handle< std::vector<TTTrack<Ref_PixelDigi_> > > pixelDigiTTTracks;
  iEvent.getByLabel(L1TrackTag_, pixelDigiTTTracks);
  std::cout << "pixelDigiTTracks size =  " << pixelDigiTTTracks->size() << std::endl;

// #ifdef THIS_IS_AN_EVENT_EXAMPLE
//    Handle<ExampleData> pIn;
//    iEvent.getByLabel("example",pIn);
// #endif
   
// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//    ESHandle<SetupData> pSetup;
//    iSetup.get<SetupRecord>().get(pSetup);
// #endif

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
