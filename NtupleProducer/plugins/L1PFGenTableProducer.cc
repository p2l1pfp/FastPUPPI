// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TPFUtils.h"

#include <algorithm>

class L1PFGenTableProducer : public edm::global::EDProducer<>  {
    public:
        explicit L1PFGenTableProducer(const edm::ParameterSet&);
        ~L1PFGenTableProducer();

    private:
        virtual void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;


        std::string name_;
        edm::EDGetTokenT<std::vector<reco::GenParticle>> genparticles_;
        StringCutObjectSelector<reco::GenParticle> sel_;

};

L1PFGenTableProducer::L1PFGenTableProducer(const edm::ParameterSet& iConfig) :
    name_(iConfig.getParameter<std::string>("name")),
    genparticles_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("src"))),
    sel_(iConfig.getParameter<std::string>("cut"), true)
{
    produces<nanoaod::FlatTable>();

}

L1PFGenTableProducer::~L1PFGenTableProducer() { }

// ------------ method called for each event  ------------
    void
L1PFGenTableProducer::produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{
    edm::Handle<std::vector<reco::GenParticle>> genparticles;
    iEvent.getByToken(genparticles_, genparticles);

    std::vector<const reco::GenParticle *> selected;


    for (const reco::GenParticle &gen : *genparticles)
        if(sel_(gen)) selected.push_back(&gen);

    const float bz = 3.8112;

    // create the table
    unsigned int ncands = selected.size();
    auto out = std::make_unique<nanoaod::FlatTable>(ncands, name_, false, true);

    std::vector<float> vals_caloeta, vals_calophi;
    vals_caloeta.resize(ncands);
    vals_calophi.resize(ncands);

    for (unsigned int i = 0; i < ncands; ++i) {
        math::XYZTLorentzVector vertex(selected[i]->vx(),selected[i]->vy(),selected[i]->vz(),0.);
        auto caloetaphi = l1tpf::propagateToCalo(selected[i]->p4(),vertex,selected[i]->charge(),bz);
        vals_caloeta[i] = caloetaphi.first;
        vals_calophi[i] = caloetaphi.second;
    }

    out->addColumn<float>("caloeta", vals_caloeta, "");
    out->addColumn<float>("calophi", vals_calophi, "");

    // save to the event branches
    iEvent.put(std::move(out));

   
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PFGenTableProducer);
