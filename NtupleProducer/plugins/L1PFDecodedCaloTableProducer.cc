// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/L1TParticleFlow/interface/PFCluster.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloCrystalCluster.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"


#include <algorithm>

class L1PFDecodedCaloTableProducer : public edm::global::EDProducer<>  {
    public:
        explicit L1PFDecodedCaloTableProducer(const edm::ParameterSet&);
        ~L1PFDecodedCaloTableProducer();

    private:
        virtual void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;


        std::string name_;
        edm::EDGetTokenT<l1t::PFClusterCollection> clusters_;
        StringCutObjectSelector<l1t::PFCluster> sel_;

};

L1PFDecodedCaloTableProducer::L1PFDecodedCaloTableProducer(const edm::ParameterSet& iConfig) :
    name_(iConfig.getParameter<std::string>("name")),
    clusters_(consumes<l1t::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("src"))),
    sel_(iConfig.getParameter<std::string>("cut"), true)
{
    produces<nanoaod::FlatTable>();

}

L1PFDecodedCaloTableProducer::~L1PFDecodedCaloTableProducer() { }

// ------------ method called for each event  ------------
    void
L1PFDecodedCaloTableProducer::produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{
    edm::Handle<l1t::PFClusterCollection> clusters;
    iEvent.getByToken(clusters_, clusters);

    std::vector<const l1t::PFCluster *> selected;


    for (const l1t::PFCluster &cl : *clusters)
        if(sel_(cl)) selected.push_back(&cl);

    
    // create the table
    unsigned int ncands = selected.size();
    auto out = std::make_unique<nanoaod::FlatTable>(ncands, name_, false, true);

    std::vector<float> vals_empt, vals_srrTot, vals_hwSrrTot, vals_meanz, vals_hwMeanZ, vals_hoe, vals_piIdProb, vals_PuIdProb, vals_EmIdProb, vals_caloIso, vals_showerShape; 
    vals_empt.resize(ncands);
    vals_srrTot.resize(ncands);
    vals_hwSrrTot.resize(ncands);
    vals_meanz.resize(ncands);
    vals_hwMeanZ.resize(ncands);
    vals_hoe.resize(ncands);
    vals_piIdProb.resize(ncands);
    vals_PuIdProb.resize(ncands);
    vals_EmIdProb.resize(ncands);
    vals_caloIso.resize(ncands);
    vals_showerShape.resize(ncands);
    
    for (unsigned int i = 0; i < ncands; ++i) {
        const auto cand = selected[i];
        auto obj = cand->caloDigiObj();

        if(auto digi = std::get_if<l1ct::EmCaloObj>(&obj)){
            vals_srrTot[i] =   digi->floatSrrTot();
            vals_hwSrrTot[i] = digi->hwSrrTot.to_float();
            vals_meanz[i] =    digi->floatMeanZ();
            vals_hwMeanZ[i] =  digi->hwMeanZ.to_float();
            vals_hoe[i] =      digi->floatHoe();
            vals_piIdProb[i] = digi->floatPiProb();
            vals_PuIdProb[i] = digi->floatPuProb();
            vals_EmIdProb[i] = digi->floatEmProb();

            const l1tp2::CaloCrystalCluster *crycl = dynamic_cast<const l1tp2::CaloCrystalCluster *>(cand->constituentsAndFractions().front().first.get());
            if(crycl) {
                vals_caloIso[i] = crycl->isolation();
                vals_showerShape[i] = crycl->e2x5() / crycl->e5x5();
            }
        } else if(auto digi = std::get_if<l1ct::HadCaloObj>(&obj)){
            vals_empt[i] = digi->floatEmPt();

            vals_srrTot[i] = digi->floatSrrTot();
            vals_hwSrrTot[i] = digi->hwSrrTot.to_float();
            vals_meanz[i] = digi->floatMeanZ();
            vals_hwMeanZ[i] = digi->hwMeanZ.to_float();
            vals_hoe[i] = digi->floatHoe();
            vals_piIdProb[i] = digi->floatPiProb();
            vals_PuIdProb[i] = digi->floatPuProb();
            vals_EmIdProb[i] = digi->floatEmProb();

        }
    }


    out->addColumn<float>("empt", vals_empt, "");
    out->addColumn<float>("srrTot", vals_srrTot, "");
    out->addColumn<float>("hwSrrTot", vals_hwSrrTot, "");
    out->addColumn<float>("meanz", vals_meanz, "");
    out->addColumn<float>("hwMeanZ", vals_hwMeanZ, "");
    out->addColumn<float>("hoe", vals_hoe, "");
    out->addColumn<float>("piIdProb", vals_piIdProb, "");
    out->addColumn<float>("PuIdProb", vals_PuIdProb, "");
    out->addColumn<float>("EmIdProb", vals_EmIdProb, "");
    out->addColumn<float>("caloIso", vals_caloIso, "");
    out->addColumn<float>("showerShape", vals_showerShape, "");

    // save to the event branches
    iEvent.put(std::move(out));

   
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PFDecodedCaloTableProducer);
