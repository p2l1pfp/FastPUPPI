////////////////////////////// L1TrackTableProducer ///////////////////////////
// Original author: G. Karathanasis, CERN gkaratha@cern.ch
// Saves all L1 tracks in FastPUPPI. Has hardcoded the basic variables and assoc
// iations with MC particles. Additional variables can be added in runtime like
// in nano. Same with cuts
 

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/Candidate.h"


#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "SimDataFormats/Associations/interface/TTTrackAssociationMap.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include <algorithm>




class L1TrackTableProducer : public edm::global::EDProducer<>  {
    public:
        explicit L1TrackTableProducer(const edm::ParameterSet&);
        ~L1TrackTableProducer();

    private:
        virtual void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

        typedef TTTrack<Ref_Phase2TrackerDigi_> l1Track;
        typedef TTTrackAssociationMap<Ref_Phase2TrackerDigi_> l1TrackingParticleMap;
        edm::EDGetTokenT< std::vector<l1Track> > l1TrkToken_;
        edm::EDGetTokenT< l1TrackingParticleMap > l1TrackingPartToken_;
        StringCutObjectSelector<l1Track> tk_sel_;
        std::string name_;

        struct ExtraVar {
            std::string name, expr;
            StringObjectFunction<l1Track> func;
            ExtraVar(const std::string & n, const std::string & expr) : name(n), expr(expr), func(expr, true) {}
        };
        std::vector<ExtraVar> extraVars_;
        

};



L1TrackTableProducer::L1TrackTableProducer(const edm::ParameterSet& iConfig) :
    l1TrkToken_{ consumes<std::vector<l1Track>> ( iConfig.getParameter <edm::InputTag> ("tracks") ) },
    l1TrackingPartToken_{ consumes<l1TrackingParticleMap>( iConfig.getParameter <edm::InputTag> ("trackingParticleMap") )},
    tk_sel_{iConfig.getParameter<std::string>("selection"), true},
    name_{iConfig.getParameter<std::string>("name")}
{
   if (iConfig.existsAs<edm::ParameterSet>("variables")) {
        edm::ParameterSet vars = iConfig.getParameter<edm::ParameterSet>("variables");
        auto morenames = vars.getParameterNamesForType<std::string>();
        for (const std::string & name : morenames) {
            extraVars_.emplace_back(name, vars.getParameter<std::string>(name));
        }
    }
    produces<nanoaod::FlatTable>(name_);    
 }



L1TrackTableProducer::~L1TrackTableProducer() { }

// ------------ method called for each event  ------------
    void
L1TrackTableProducer::produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{
    // get input
    edm::Handle<std::vector<l1Track>> l1_tracks;
    edm::Handle<l1TrackingParticleMap> l1_tracking_parts_map;
    iEvent.getByToken(l1TrkToken_, l1_tracks);
    iEvent.getByToken(l1TrackingPartToken_, l1_tracking_parts_map);
    std::vector<bool> isLooseGenuineTP,isGenuineTP,isUnknownTP,isCombinatoricTP;
    std::vector<float> trk_pt,trk_eta,trk_phi,trk_q,tp_pt,tp_eta,tp_phi;
    std::vector<int> tp_pdgid;
    std::vector<float> extra_var_values[extraVars_.size()];
  
 
    for (size_t itrk=0; itrk<l1_tracks->size(); itrk++) {
      // get and select
      edm::Ptr<l1Track> l1trk_ptr(l1_tracks, itrk);
      if (!tk_sel_(*l1trk_ptr))
         continue;
      trk_pt.emplace_back(l1trk_ptr->momentum().perp());
      trk_eta.emplace_back(l1trk_ptr->momentum().eta());
      trk_phi.emplace_back(l1trk_ptr->momentum().phi());
      trk_q.emplace_back(l1trk_ptr->rInv()/fabs(l1trk_ptr->rInv()));
      for (size_t ivar =0; ivar< extraVars_.size(); ivar++) {
         extra_var_values[ivar].emplace_back( extraVars_[ivar].func(*l1trk_ptr) ); 
      }
      isLooseGenuineTP.emplace_back(l1_tracking_parts_map->isLooselyGenuine(l1trk_ptr));
      isGenuineTP.emplace_back(l1_tracking_parts_map->isGenuine(l1trk_ptr));
      isUnknownTP.emplace_back(l1_tracking_parts_map->isUnknown(l1trk_ptr));
      isCombinatoricTP.emplace_back(l1_tracking_parts_map->isCombinatoric(l1trk_ptr));
      if ( l1_tracking_parts_map->isLooselyGenuine(l1trk_ptr) || l1_tracking_parts_map->isGenuine(l1trk_ptr) ){
         edm::Ptr<TrackingParticle> trackingpart_ptr = l1_tracking_parts_map->findTrackingParticlePtr(l1trk_ptr);
         tp_pt.emplace_back( trackingpart_ptr.isNull() ? -99 : trackingpart_ptr->p4().pt() );
         tp_eta.emplace_back( trackingpart_ptr.isNull() ? -99 : trackingpart_ptr->p4().eta());
         tp_phi.emplace_back( trackingpart_ptr.isNull() ? -99 : trackingpart_ptr->p4().phi());
         tp_pdgid.emplace_back( trackingpart_ptr.isNull() ? -99 : trackingpart_ptr->pdgId());     
      } else{
         tp_pt.emplace_back(-99);
         tp_eta.emplace_back(-99);
         tp_phi.emplace_back(-99);
         tp_pdgid.emplace_back(-99);
      }
          
    }
    // create the table
    unsigned int ncands = isLooseGenuineTP.size();
    auto out = std::make_unique<nanoaod::FlatTable>(ncands, name_, false);
    out->addColumn<float>("pt", trk_pt, "l1 track pt");
    out->addColumn<float>("eta", trk_eta, "l1 track eta");
    out->addColumn<float>("phi", trk_phi, "l1 track phi");
    out->addColumn<int>("charge", trk_q, "l1 track charge");
    out->addColumn<bool>("isLooseGenuineTP", isLooseGenuineTP, "track loosely matched to TP");
    out->addColumn<bool>("isGenuineTP", isGenuineTP, "track matched to TP");
    out->addColumn<bool>("isUnknownTP", isUnknownTP, "track not matched to TP");
    out->addColumn<bool>("isCombinatoricTP", isCombinatoricTP, "combinatoric track");
    out->addColumn<float>("tp_pt", tp_pt, "pt of tracking particle");
    out->addColumn<float>("tp_eta", tp_eta, "eta of tracking particle");
    out->addColumn<float>("tp_phi", tp_phi, "phi of tracking particle");
    out->addColumn<int>("tp_pdgid", tp_pdgid, "pdgId of tracking particle");

    // write arbitrary number of variables
    for (size_t ivar =0; ivar< extraVars_.size(); ivar++) {
       out->addColumn<float>(extraVars_[ivar].name, extra_var_values[ivar], extraVars_[ivar].expr);
    }

    iEvent.put(std::move(out), name_);
}



//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TrackTableProducer);
