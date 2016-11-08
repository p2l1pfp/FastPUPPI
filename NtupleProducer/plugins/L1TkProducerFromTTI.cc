#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#ifdef HASL1TK
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"     
#endif

#include "FastPUPPI/NtupleProducer/interface/L1TPFParticle.h"
#include "FastPUPPI/NtupleProducer/interface/L1TPFUtils.h"

namespace l1tpf {
    class TkProducerFromTTI : public edm::EDProducer {
        public:
            explicit TkProducerFromTTI(const edm::ParameterSet&) ;
            ~TkProducerFromTTI() {}

        private:
            edm::EDGetToken L1TrackTag_;
            float fBz_;

            virtual void produce(edm::Event&, const edm::EventSetup&) override;

            virtual void beginRun(edm::Run const&, edm::EventSetup const& iSetup) override {
                edm::ESHandle<MagneticField> magneticField;
                iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
                fBz_ = magneticField->inTesla(GlobalPoint(0,0,0)).z();
            }

    }; // class
} // namespace

l1tpf::TkProducerFromTTI::TkProducerFromTTI(const edm::ParameterSet&iConfig)
{
#ifdef HASL1TK
    L1TrackTag_ = consumes<std::vector<TTTrack<Ref_PixelDigi_>>>(iConfig.getParameter<edm::InputTag>("L1TrackTag"));
#endif
    produces<std::vector<l1tpf::Particle>>();
}


void 
l1tpf::TkProducerFromTTI::produce(edm::Event & iEvent, const edm::EventSetup &) 
{
  std::unique_ptr<std::vector<l1tpf::Particle>> out(new std::vector<l1tpf::Particle>());
   
#ifdef HASL1TK
  /// ----------------TRACK INFO-------------------
  /// Stealing Jia Fu's code!
  /// https://github.com/jiafulow/SLHCL1TrackTriggerSimulations/blob/master/NTupleTools/src/NTupleTTTracks.cc
  typedef std::vector<TTTrack<Ref_PixelDigi_> >         vec_track;

  edm::Handle< vec_track > pixelDigiTTTracks;
  iEvent.getByToken(L1TrackTag_, pixelDigiTTTracks);
  if (pixelDigiTTTracks.isValid()) {
    
    edm::LogInfo("NTupleTracks") << "Size: " << pixelDigiTTTracks->size();
    
    unsigned nPar = 4;
    for (vec_track::const_iterator it = pixelDigiTTTracks->begin(); it != pixelDigiTTTracks->end(); ++it) {
      const GlobalVector&          momentum = it->getMomentum(nPar);
      const GlobalPoint&           poca     = it->getPOCA(nPar);  // point of closest approach
      bool iIsEle=false;
      double mass  = iIsEle ? 0.0005 : 0.139; 
      double energy=sqrt((mass*mass)+(momentum.mag()*momentum.mag()));
      const math::XYZTLorentzVector tMom(momentum.x(),momentum.y(),momentum.z(),energy);
      const math::XYZTLorentzVector tVtx(poca.x()     ,poca.y()     ,poca.z(),0.);
      
      double charge = it->getRInv()/fabs(it->getRInv());

      //Check below
      std::vector<double> lVars;
      l1tpf::propagate(1,lVars,tMom,tVtx,charge,fBz_);
      float trkEcalEta = lVars[4];
      float trkEcalPhi = lVars[5];
 
      int pQuality = 0;
      float nPS = 0.;     // number of stubs in PS modules
      float nstubs = 0;
      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > >  theStubs = it-> getStubRefs() ;
      // loop over the stubs
      for (unsigned int istub=0; istub<(unsigned int)theStubs.size(); istub++) {
	nstubs ++;
	StackedTrackerDetId detIdStub( theStubs.at(istub)->getDetId() );
	bool isPS = true;//theStackedGeometry -> isPSModule( detIdStub );
	if (isPS) nPS ++;
      } // end loop over stubs
      if(nPS           > 2) pQuality++;
      if(nstubs        > 3) pQuality+=2;
      if(it->getChi2() < 100.) pQuality += 4; 

      out->emplace_back(momentum.perp(), momentum.eta(),momentum.phi(), 0.137, 
                        0, 0.0, poca.z(), trkEcalEta, trkEcalPhi, charge, pQuality);
      out->back().setVertex(reco::Particle::Point(poca.x(),poca.y(),poca.z()));

      //std::cout << "=== Ecal surface eta : " << lVars[4] << " -- phi : " << lVars[5] << " -- R: " << lVars[6] << " -- eta/phi : " << momentum.eta() << " -- " << momentum.phi() << " -- pt : " << momentum.perp() << " -- " << sqrt(lVars[0]*lVars[0]+lVars[1]*lVars[1])   << std::endl;
      // Std::cout << "track info = " << momentum.perp() << "," << momentum.eta() << "," << momentum.phi() << "; poca z = " << poca.z() << std::endl;        
      // if (it - pixelDigiTTTracks->begin() > 10) break; // just for testing so cut it off for now
    }  
  }
#endif
  iEvent.put(std::move(out));
}
using l1tpf::TkProducerFromTTI;
DEFINE_FWK_MODULE(TkProducerFromTTI);
