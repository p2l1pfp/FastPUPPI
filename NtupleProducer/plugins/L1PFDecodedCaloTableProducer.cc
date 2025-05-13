// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/L1TParticleFlow/interface/PFCluster.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloCrystalCluster.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"


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

    std::vector<float> vals_empt, vals_srrTot, vals_hwSrrTot, vals_meanz, 
        vals_hwMeanZ, vals_hoe, vals_piIdProb, vals_PuIdProb, vals_EmIdProb, 
        vals_caloIso, vals_showerShape; 
     std::vector<float> vals_showerlength,
        vals_coreshowerlength, vals_emf, vals_hw_emf, vals_abseta, 
        vals_hw_abseta, 
        vals_hw_meanz, 
        vals_sigmaetaeta, vals_hw_sigmaetaeta, 
        vals_sigmaphiphi, vals_hw_sigmaphiphi, vals_sigmazz, vals_hw_sigmazz;

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
    vals_showerlength.resize(ncands);
    vals_coreshowerlength.resize(ncands);
    vals_emf.resize(ncands);
    vals_hw_emf.resize(ncands);
    vals_abseta.resize(ncands);
    vals_hw_abseta.resize(ncands);
    vals_hw_meanz.resize(ncands);
    vals_sigmaetaeta.resize(ncands);
    vals_hw_sigmaetaeta.resize(ncands);
    vals_sigmaphiphi.resize(ncands);
    vals_hw_sigmaphiphi.resize(ncands);
    vals_sigmazz.resize(ncands);
    vals_hw_sigmazz.resize(ncands);


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

            const l1t::HGCalMulticluster *hgcalcl = dynamic_cast<const l1t::HGCalMulticluster *>(cand->constituentsAndFractions().front().first.get());
            if(hgcalcl) {
                static constexpr float ETAPHI_LSB = M_PI / 720;
                static constexpr float SIGMAZZ_LSB = 778.098 / (1 << 7);
                static constexpr float SIGMAPHIPHI_LSB = 0.12822 / (1 << 7);
                static constexpr float SIGMAETAETA_LSB = 0.148922 / (1 << 5);

                ap_uint<6> w_showerlenght = hgcalcl->showerLength();
                ap_uint<6> w_coreshowerlenght = hgcalcl->coreShowerLength();
                ap_uint<8> w_emf = std::min(round(hgcalcl->eot() * 256), float(255.));
                ap_uint<10> w_abseta = round(fabs(hgcalcl->eta()) / ETAPHI_LSB);
                ap_ufixed<12, 11, AP_RND_CONV, AP_SAT> w_meanz_f = fabs(hgcalcl->zBarycenter()) - 320;  // LSB = 0.5cm
                ap_uint<12> w_meanz = w_meanz_f.range();
                ap_uint<5> w_sigmaetaeta = round(hgcalcl->sigmaEtaEtaTot() / SIGMAETAETA_LSB);
                ap_uint<7> w_sigmaphiphi = round(hgcalcl->sigmaPhiPhiTot() / SIGMAPHIPHI_LSB);
                ap_uint<7> w_sigmazz = round(hgcalcl->sigmaZZ() / SIGMAZZ_LSB);

                vals_showerlength[i] = w_showerlenght.to_int();
                vals_coreshowerlength[i] = w_coreshowerlenght.to_int();
                vals_emf[i] = w_emf / 256.;
                vals_hw_emf[i] = w_emf.to_float();
                vals_abseta[i] = w_abseta * ETAPHI_LSB;
                vals_hw_abseta[i] = w_abseta.to_float();
                vals_hw_meanz[i] = w_meanz*0.5;
                vals_sigmaetaeta[i] = w_sigmaetaeta * SIGMAETAETA_LSB;
                vals_hw_sigmaetaeta[i] = w_sigmaetaeta.to_float();
                vals_sigmaphiphi[i] = w_sigmaphiphi * SIGMAPHIPHI_LSB;
                vals_hw_sigmaphiphi[i] = w_sigmaphiphi.to_float();
                vals_sigmazz[i] = w_sigmazz * SIGMAZZ_LSB;
                vals_hw_sigmazz[i] = w_sigmazz.to_float();
            }
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

    out->addColumn<float>("showerlength", vals_showerlength, "");
    out->addColumn<float>("coreshowerlength", vals_coreshowerlength, "");
    out->addColumn<float>("emf", vals_emf, "");
    out->addColumn<float>("hwEmf", vals_hw_emf, "");
    out->addColumn<float>("abseta", vals_abseta, "");
    out->addColumn<float>("hwAbseta", vals_hw_abseta, "");
    out->addColumn<float>("hwFPMeanz", vals_hw_meanz, "");    
    out->addColumn<float>("sigmaetaeta", vals_sigmaetaeta, "");
    out->addColumn<float>("hwSigmaetaeta", vals_hw_sigmaetaeta, "");
    out->addColumn<float>("sigmaphiphi", vals_sigmaphiphi, "");
    out->addColumn<float>("hwSigmaphiphi", vals_hw_sigmaphiphi, "");
    out->addColumn<float>("sigmazz", vals_sigmazz, "");
    out->addColumn<float>("hwSigmazz", vals_hw_sigmazz, "");

    // save to the event branches
    iEvent.put(std::move(out));

   
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PFDecodedCaloTableProducer);
