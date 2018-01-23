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
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FastPUPPI/NtupleProducer/interface/corrector.hh"
#include "FastPUPPI/NtupleProducer/interface/metanalyzer.hh"
#include "FastPUPPI/NtupleProducer/interface/jetanalyzer.hh"
#include "FastPUPPI/NtupleProducer/interface/isoanalyzer.hh"
#include "FastPUPPI/NtupleProducer/interface/SimpleCalibrations.h"
#include "FastPUPPI/NtupleProducer/interface/DiscretePF.h"
#include "FastPUPPI/NtupleProducer/interface/AlternativePF.h"
#include "FastPUPPI/NtupleProducer/interface/BitwisePF.h"

#include "DataFormats/Math/interface/deltaR.h"


// ROOT classes
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObject.h>

typedef math::XYZTLorentzVector                        LorentzVector;
typedef std::vector< reco::PFCandidate >               PFOutputCollection;
typedef std::vector<l1tpf::Particle>                   L1PFCollection;

//--------------------------------------------------------------------------------------------------
class NtupleProducer : public edm::EDProducer {
public:
  explicit NtupleProducer(const edm::ParameterSet&);
  ~NtupleProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  const edm::InputTag L1TrackTag_;
  const std::vector<edm::InputTag> CaloClusterTags_;
  const std::vector<edm::InputTag> EmClusterTags_;
  const bool correctCaloEnergies_;
  const float rawCaloPtMin_, rawEmCaloPtMin_;
  const edm::InputTag MuonTPTag_;
  const edm::InputTag GenParTag_, GenOriginTag_;
  const std::string CorrectorTag_;
  const unsigned int CorrectorEmfBins_;
  const double CorrectorEmfMax_;
  const std::string ECorrectorTag_;
  const std::string TrackResTag_;
  const std::string EleResTag_;
  const std::string PionResTag_;
 
private:
  inline std::string getFilePath(const edm::ParameterSet & pset, const std::string & name) const {
    std::string ret = pset.getParameter<std::string>(name);
    if (ret[0] != '/') ret = edm::FileInPath(ret).fullPath();
    return ret;
  }

  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void genMatch(std::vector<double> &iGenVars, int iType, double iEta, double iPhi, double iPt, const reco::GenParticleCollection &iGenParticles);
  void addPF(const std::vector<l1tpf::Particle> &iCandidates, const std::string &iLabel, edm::Event& iEvent);
  void addUInt(unsigned int value,std::string iLabel,edm::Event& iEvent);

  edm::EDGetTokenT<reco::GenParticleCollection>   TokGenPar_;
  edm::EDGetTokenT<math::XYZPointF>               TokGenOrigin_;
  edm::EDGetTokenT<L1PFCollection>                TokL1TrackTPTag_;
  std::vector<edm::EDGetTokenT<L1PFCollection>>   TokCaloClusterTags_, TokEmClusterTags_;
  edm::EDGetTokenT<L1PFCollection>                TokMuonTPTag_;
  double trkPt_, trkMaxChi2_;
  unsigned trkMinStubs_;
  bool   metRate_;
  double etaCharged_;
  double puppiPtCut_, puppiDr_;
  double vtxRes_;
  bool   vtxAdaptiveCut_;
  l1tpf_int::PFAlgo::VertexAlgo vtxAlgo_;  
  corrector* corrector_;
  corrector* ecorrector_;
  // simplified corrections
  l1tpf::SimpleCorrEm simpleCorrEm_;
  l1tpf::SimpleCorrHad simpleCorrHad_;
  metanalyzer* metanalyzer_;
  jetanalyzer* jetanalyzer_;
  isoanalyzer* isoanalyzer_;
  // simplified resolutions
  l1tpf::SimpleResol simpleResolEm_, simpleResolHad_, simpleResolTrk_;
  // discretized version
  l1tpf_int::RegionMapper l1regions_;
  std::unique_ptr<l1tpf_int::PFAlgo> l1pfalgo_;
  // fill track tree to TFileService
  int fTrackTree;
  // debug flag
  int fDebug;
  float fDebugEta, fDebugPhi, fDebugR;

  // Region dump
  FILE *fRegionDump;
   
  // declare variables for output file
  std::string           fOutputName;
  TFile                 *fOutputFile;
  TH1D                  *fTotalEvents;

  TTree                 *fTrkInfoTree;
  float runNum, lumiSec, evtNum;
  float trkNum;
  float trkPx, trkPz, trkPy, trkPt, trkEta, trkPhi, trkz0, trkdz, trkd0;    
  float trkEcalEta, trkEcalPhi, trkEcalR;
  float genPt, genEta, genPhi, genId;
};

//
// constructors and destructor
//
NtupleProducer::NtupleProducer(const edm::ParameterSet& iConfig):
  L1TrackTag_           (iConfig.getParameter<edm::InputTag>("L1TrackTag")),
  CaloClusterTags_      (iConfig.getParameter<std::vector<edm::InputTag>>("CaloClusterTags")),
  EmClusterTags_        (iConfig.getParameter<std::vector<edm::InputTag>>("EmClusterTags")),
  correctCaloEnergies_  (iConfig.getParameter<bool>("correctCaloEnergies")),
  rawCaloPtMin_         (correctCaloEnergies_   ? iConfig.getParameter<double>("rawCaloPtMin") : 0.0), 
  rawEmCaloPtMin_       (correctCaloEnergies_ && !EmClusterTags_.empty() ? iConfig.getParameter<double>("rawEmCaloPtMin") : 0.0), 
  MuonTPTag_            (iConfig.getParameter<edm::InputTag>("MuonTPTag")),
  GenParTag_            (iConfig.getParameter<edm::InputTag>("genParTag")),
  GenOriginTag_         (iConfig.getParameter<edm::InputTag>("genOriginTag")),
  CorrectorTag_         (getFilePath(iConfig,"corrector")),
  CorrectorEmfBins_     (iConfig.getParameter<uint32_t>("correctorEmfBins")),
  CorrectorEmfMax_      (iConfig.getParameter<double>("correctorEmfMax")),
  ECorrectorTag_        (getFilePath(iConfig,"ecorrector")),
  TrackResTag_          (getFilePath(iConfig,"trackres")),
  EleResTag_            (getFilePath(iConfig,"eleres")),
  PionResTag_           (getFilePath(iConfig,"pionres")),
  trkPt_                (iConfig.getParameter<double>       ("trkPtCut")),
  trkMaxChi2_           (iConfig.getParameter<double>       ("trkMaxChi2")),
  trkMinStubs_          (iConfig.getParameter<unsigned>     ("trkMinStubs")),
  metRate_              (iConfig.getParameter<bool>         ("metRate")),
  etaCharged_           (iConfig.getParameter<double>       ("etaCharged")),
  puppiPtCut_           (iConfig.getParameter<double>       ("puppiPtCut")),
  puppiDr_              (iConfig.getParameter<double>       ("puppiDr")),
  vtxRes_               (iConfig.getParameter<double>       ("vtxRes")),
  vtxAdaptiveCut_       (iConfig.getParameter<bool>         ("vtxAdaptiveCut")),
  vtxAlgo_              (l1tpf_int::PFAlgo::OldVtxAlgo), 
  simpleCorrEm_         (iConfig, "simpleCorrEm"),
  simpleCorrHad_        (iConfig, "simpleCorrHad"),
  l1regions_            (iConfig),
  l1pfalgo_             (nullptr),
  fTrackTree            (iConfig.getUntrackedParameter<int>("fillTrackTree",0)),
  fDebug                (iConfig.getUntrackedParameter<int>("debug",0)),
  fDebugEta             (iConfig.getUntrackedParameter<double>("debugEta",0)),
  fDebugPhi             (iConfig.getUntrackedParameter<double>("debugPhi",0)),
  fDebugR               (iConfig.getUntrackedParameter<double>("debugR",-1)),
  fRegionDump           (nullptr),
  fOutputName           (iConfig.getUntrackedParameter<std::string>("outputName", "ntuple.root")),
  fOutputFile           (0),
  fTotalEvents          (0),
  fTrkInfoTree          (0)
{
  std::string vtxAlgo = iConfig.getParameter<std::string>("vtxAlgo");
  if      (vtxAlgo == "TP") vtxAlgo_ = l1tpf_int::PFAlgo::TPVtxAlgo;
  else if (vtxAlgo != "old") throw cms::Exception("Configuration") << "Unsupported vtxAlgo " << vtxAlgo << "\n";
  //now do what ever other initialization is needed
  corrector_  = new corrector(CorrectorTag_,CorrectorEmfBins_,CorrectorEmfMax_,fDebug);
  ecorrector_ = new corrector(ECorrectorTag_,1,1.0,fDebug);
  simpleResolEm_ = l1tpf::SimpleResol(iConfig,"simpleResolEm", true);
  simpleResolHad_ = l1tpf::SimpleResol(iConfig,"simpleResolHad", true);
  simpleResolTrk_ = l1tpf::SimpleResol(iConfig,"simpleResolTrk", true);
  edm::ParameterSet linkcfg = iConfig.getParameter<edm::ParameterSet>("linking");
  auto algo = linkcfg.getParameter<std::string>("algo");
  if (algo == "PFAlgo3") {
      l1pfalgo_.reset(new l1tpf_int::PFAlgo3(iConfig));
  } else if (algo == "PFAlgoOld") {
      l1pfalgo_.reset(new l1tpf_int::PFAlgo(iConfig));
  } else if (algo == "BitwisePF") {
      l1pfalgo_.reset(new l1tpf_int::BitwisePF(iConfig));
  } else throw cms::Exception("Configuration", "Unsupported PFAlgo");

  std::string regionDumpFile = iConfig.getUntrackedParameter<std::string>("regionDumpFileName", "");
  if (!regionDumpFile.empty()) fRegionDump = fopen(regionDumpFile.c_str(), "wb");

  if (fOutputName.empty()) {
      metanalyzer_ = nullptr;
      jetanalyzer_ = nullptr;
      isoanalyzer_ = nullptr;
  } else {
      metanalyzer_ = new metanalyzer("MetFile.root");
      jetanalyzer_ = new jetanalyzer("JetFile.root");
      isoanalyzer_ = new isoanalyzer("IsoFile.root");
  }
  produces<PFOutputCollection>("RawEmCalo");
  produces<PFOutputCollection>("RawCalo");
  produces<PFOutputCollection>("EmCalo");
  produces<PFOutputCollection>("Calo");
  produces<PFOutputCollection>("TK");
  produces<PFOutputCollection>("TKVtx");
  produces<PFOutputCollection>("PF");
  produces<PFOutputCollection>("PFDiscarded");
  produces<PFOutputCollection>("Puppi");
  produces<float>("z0");
  produces<float>("alphaCMed"); produces<float>("alphaCRms"); produces<float>("alphaFMed"); produces<float>("alphaFRms");

  TokGenPar_       = consumes<reco::GenParticleCollection>( GenParTag_    );
  TokGenOrigin_    = consumes<math::XYZPointF>( GenOriginTag_    );
  TokL1TrackTPTag_ = consumes<L1PFCollection>( L1TrackTag_  );
  for (const edm::InputTag &tag : CaloClusterTags_) {
    TokCaloClusterTags_.push_back(consumes<L1PFCollection>(tag));
  }
  for (const edm::InputTag &tag : EmClusterTags_) {
    TokEmClusterTags_.push_back(consumes<L1PFCollection>(tag));
  }
  TokMuonTPTag_    = consumes<L1PFCollection>( MuonTPTag_  );
  for (int tot = 0; tot <= 1; ++tot) {
      for (int i = 0; i < l1tpf_int::Region::n_input_types; ++i) {
          produces<unsigned int>(std::string(tot ? "totNL1" : "maxNL1")+l1tpf_int::Region::inputTypeName(i));
      }
      for (int i = 0; i < l1tpf_int::Region::n_output_types; ++i) {
          produces<unsigned int>(std::string(tot ? "totNL1PF"    : "maxNL1PF"   )+l1tpf_int::Region::outputTypeName(i));
          produces<unsigned int>(std::string(tot ? "totNL1Puppi" : "maxNL1Puppi")+l1tpf_int::Region::outputTypeName(i));
      }
  }
}

NtupleProducer::~NtupleProducer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   if (fRegionDump) fclose(fRegionDump);
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
NtupleProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //fOutputFile->cd();
  if (fTotalEvents) fTotalEvents->Fill(1);  
  using namespace edm;

  edm::Handle<reco::GenParticleCollection> hGenParProduct;
  iEvent.getByToken(TokGenPar_,hGenParProduct);
  const reco::GenParticleCollection & genParticles = *hGenParProduct;  
  edm::Handle<math::XYZPointF> hGenOrigin;
  iEvent.getByToken(TokGenOrigin_,hGenOrigin);
  const math::XYZPointF & genOrigin = *hGenOrigin;
  

  l1regions_.clear(); 
  /// ----------------TRACK INFO-------------------
  edm::Handle<std::vector<l1tpf::Particle>> hl1tks;
  iEvent.getByLabel(L1TrackTag_, hl1tks);
  std::vector<l1tpf::Particle> l1tks(*hl1tks);
  trkNum = 0;
  for (l1tpf::Particle & tk : l1tks) { // no const & since we modify it to set the sigma
      // the combiner knows the sigma, the track producer doesn't
      tk.setSigma(simpleResolTrk_(tk.pt(), std::abs(tk.eta())));
      tk.setCaloSigma(simpleResolHad_(tk.pt(), std::abs(tk.eta())));
      if (fDebug > 1) printf("tk %7.2f eta %+4.2f has sigma %4.2f  sigma/pt %5.3f, calo sigma/pt %5.3f, stubs %2d, chi2 %7.1f\n", tk.pt(), tk.eta(), tk.sigma(), tk.sigma()/tk.pt(), tk.caloSigma()/tk.pt(), int(tk.quality()),tk.normalizedChi2());
      // adding objects to PF
      if (fDebugR > 0 && deltaR(tk.eta(),tk.phi(),fDebugEta,fDebugPhi) > fDebugR) continue;
      if(tk.pt() > trkPt_ && unsigned(tk.quality()) >= trkMinStubs_ && tk.normalizedChi2() < trkMaxChi2_ ) l1regions_.addTrack(tk);      
      /// filling the tree    
      if (!fTrkInfoTree) continue;
      trkPx  = tk.px();
      trkPz  = tk.py();
      trkPy  = tk.pz();
      trkPt  = tk.pt();
      trkEta = tk.eta();
      trkPhi = tk.phi();
      trkz0  = tk.vertex().Z();
      trkdz  = tk.vertex().Z() - genOrigin.Z();
      trkd0  = tk.vertex().Rho();
      trkEcalEta = tk.caloEta();
      trkEcalPhi = tk.caloPhi();
      trkEcalR   = -999.; //FIXME: this wasn't in the Particle class. Do we need it?
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
      trkNum++;
  }
  
  /// ----------------Muon INFO-------------------
  edm::Handle<std::vector<l1tpf::Particle>> hl1mus;
  iEvent.getByToken(TokMuonTPTag_, hl1mus);
  std::vector<l1tpf::Particle> l1mus(*hl1mus);
  for (l1tpf::Particle & mu : l1mus) { // no const & since we modify it to set the sigma
      // possible debug filtering of inputs in some area
      if (fDebugR > 0 && deltaR(mu.eta(),mu.phi(),fDebugEta,fDebugPhi) > fDebugR) continue;
      // set resolution
      mu.setSigma(simpleResolTrk_(mu.pt(), std::abs(mu.eta())));  // this needs to be updated with the muon resolutions!
      // add to inputs
      l1regions_.addMuon(mu);
  }


  /// ----------------Clustered Calo INFO-------------------
  L1PFCollection calos, emcalos;
  edm::Handle<L1PFCollection> hcals;
  for (const auto & token : TokCaloClusterTags_) {
      iEvent.getByToken(token, hcals);
      for (const auto & c : *hcals) { 
          if (c.pt() > rawCaloPtMin_) calos.push_back(c); 
      }
  }
  for (const auto & token : TokEmClusterTags_) {
      iEvent.getByToken(token, hcals);
      for (const auto & c : *hcals) { 
          if (c.pt() > rawEmCaloPtMin_) emcalos.push_back(c); 
      }
  }

  // make a copy of the uncorrected ones
  L1PFCollection lRawCalo(calos), lRawEmCalo(emcalos);
  addPF(lRawCalo,   "RawCalo",   iEvent);
  addPF(lRawEmCalo, "RawEmCalo", iEvent);

  // calibrate and do the rest
  for (l1tpf::Particle &calo : calos) {
      if (correctCaloEnergies_) {
          float ptcorr = simpleCorrHad_.empty() ? corrector_->correct(calo.pt(), calo.emEt(), calo.iEta(), calo.iPhi()) :
                                                  simpleCorrHad_(calo.pt(), std::abs(calo.eta()), calo.emEt()/calo.pt());
          calo.setPt(ptcorr);
      }
      calo.setSigma(calo.pdgId() == l1tpf::Particle::GAMMA ?
              simpleResolEm_(calo.pt(), std::abs(calo.eta())) :
              simpleResolHad_(calo.pt(), std::abs(calo.eta())) );
      if (fDebug > 1) printf("calo %7.2f eta %+4.2f has sigma %4.2f  sigma/pt %5.3f\n", calo.pt(), calo.eta(), calo.sigma(), calo.sigma()/calo.pt());
  }
  for (l1tpf::Particle &calo : emcalos) {
      if (correctCaloEnergies_) {
          float ptcorr = simpleCorrEm_.empty() ? ecorrector_->correct(0, calo.pt(), calo.iEta(), calo.iPhi()) :
                                                 simpleCorrEm_(calo.pt(), std::abs(calo.eta()));
          calo.setPt(ptcorr);
      }
      calo.setSigma(simpleResolEm_(calo.pt(), std::abs(calo.eta())));
      if (fDebug > 1) printf("emcalo %7.2f eta %+4.2f has sigma %4.2f  sigma/pt %5.3f\n", calo.pt(), calo.eta(), calo.sigma(), calo.sigma()/calo.pt());
  }

  // pass to the PF algo
  for (const l1tpf::Particle & calo : calos) {
      if (fDebugR > 0 && deltaR(calo.eta(),calo.phi(),fDebugEta,fDebugPhi) > fDebugR) continue;
      l1regions_.addCalo(calo); 
  }
  for (const l1tpf::Particle & calo : emcalos) {
      if (fDebugR > 0 && deltaR(calo.eta(),calo.phi(),fDebugEta,fDebugPhi) > fDebugR) continue;
      l1regions_.addEmCalo(calo); 
  }

  // First, get a copy of the discretized and corrected inputs, and write them out
  std::vector<l1tpf::Particle> lEmCalo  = l1regions_.fetchCalo(/*ptmin=*/0.1, /*em=*/true);
  std::vector<l1tpf::Particle> lCorrCalo = l1regions_.fetchCalo(/*ptmin=*/0.1);
  std::vector<l1tpf::Particle> lTKCands  = l1regions_.fetchTracks();
  addPF(lEmCalo ,"EmCalo"  ,iEvent);
  addPF(lCorrCalo,"Calo"    ,iEvent);
  addPF(lTKCands ,"TK"      ,iEvent);
  if (fRegionDump) {
    uint32_t run = iEvent.id().run(), lumi = iEvent.id().luminosityBlock(); uint64_t event = iEvent.id().event();
    fwrite(&run, sizeof(uint32_t), 1, fRegionDump);
    fwrite(&lumi, sizeof(uint32_t), 1, fRegionDump);
    fwrite(&event, sizeof(uint64_t), 1, fRegionDump);
    l1tpf_int::writeManyToFile(l1regions_.regions(), fRegionDump);
  } 


  // then do the vertexing, and save it out
  float z0, genZ0 = genOrigin.Z();  
  l1pfalgo_->doVertexing(l1regions_.regions(), vtxAlgo_, z0);
  iEvent.put(std::make_unique<float>(z0), "z0");
  if (fRegionDump) {
    fwrite(&z0, sizeof(float), 1, fRegionDump);
    fwrite(&genZ0, sizeof(float), 1, fRegionDump);
  } 

  // then also save the tracks with a vertex cut
  std::vector<l1tpf::Particle> lTKVtxCands = l1regions_.fetchTracks(/*ptmin=*/0.0, /*frompv=*/true);
  addPF(lTKVtxCands, "TKVtx", iEvent);

  // then run PF in each region
  for (auto & l1region : l1regions_.regions()) {
      l1pfalgo_->runPF(l1region);
      l1pfalgo_->runChargedPV(l1region, z0);
  }
  std::vector<l1tpf::Particle> lCands = l1regions_.fetch(false);
  addPF(lCands, "PF", iEvent);

  // then get our alphas (globally)
  float alphaCMed, alphaCRms, alphaFMed, alphaFRms;
  l1pfalgo_->computePuppiMedRMS(l1regions_.regions(), alphaCMed, alphaCRms, alphaFMed, alphaFRms);
  iEvent.put(std::make_unique<float>(alphaCMed), "alphaCMed"); iEvent.put(std::make_unique<float>(alphaCRms), "alphaCRms");
  iEvent.put(std::make_unique<float>(alphaFMed), "alphaFMed"); iEvent.put(std::make_unique<float>(alphaFRms), "alphaFRms");
  if (fRegionDump) {
    fwrite(&alphaCMed, sizeof(float), 1, fRegionDump);
    fwrite(&alphaCRms, sizeof(float), 1, fRegionDump);
    fwrite(&alphaFMed, sizeof(float), 1, fRegionDump);
    fwrite(&alphaFRms, sizeof(float), 1, fRegionDump);
  } 

  // then run puppi (regionally)
  for (auto & l1region : l1regions_.regions()) {
      l1pfalgo_->runPuppi(l1region, -1., alphaCMed, alphaCRms, alphaFMed, alphaFRms);
  }
  // and save puppi
  std::vector<l1tpf::Particle> lPupCands = l1regions_.fetch(true);
  addPF(lPupCands,"Puppi"   ,iEvent);

  // then get and save objects discarted in the PF algo, for debugging
  std::vector<l1tpf::Particle> lPFDisc = l1regions_.fetch(false,0.01,true);
  produces<PFOutputCollection>("PFDiscarded");

  // then go do the multiplicities
  for (int i = 0; i < l1tpf_int::Region::n_input_types; ++i) {
      auto totAndMax = l1regions_.totAndMaxInput(i);
      addUInt(totAndMax.first,  std::string("totNL1")+l1tpf_int::Region::inputTypeName(i), iEvent);
      addUInt(totAndMax.second, std::string("maxNL1")+l1tpf_int::Region::inputTypeName(i), iEvent);
  }
  for (int i = 0; i < l1tpf_int::Region::n_output_types; ++i) {
      auto totAndMaxPF = l1regions_.totAndMaxOutput(i,false);
      auto totAndMaxPuppi = l1regions_.totAndMaxOutput(i,true);
      addUInt(totAndMaxPF.first,  std::string("totNL1PF")+l1tpf_int::Region::outputTypeName(i), iEvent);
      addUInt(totAndMaxPF.second, std::string("maxNL1PF")+l1tpf_int::Region::outputTypeName(i), iEvent);
      addUInt(totAndMaxPuppi.first,  std::string("totNL1Puppi")+l1tpf_int::Region::outputTypeName(i), iEvent);
      addUInt(totAndMaxPuppi.second, std::string("maxNL1Puppi")+l1tpf_int::Region::outputTypeName(i), iEvent);
  }

 
  if (metanalyzer_) {
      metanalyzer_->clear();
      metanalyzer_->setZ(lTKCands,z0);
      metanalyzer_->setMETRecoil(2,lTKCands ,false);
      metanalyzer_->setMETRecoil(0,lCands,false);
      metanalyzer_->setMETRecoil(3,lRawCalo ,true);
      metanalyzer_->setMETRecoil(1,lCorrCalo,true);
      metanalyzer_->setMETRecoil(5,lPupCands,false);
      metanalyzer_->setMETRecoil(4,lTKVtxCands ,false);
      metanalyzer_->setGenMET(genParticles);
      metanalyzer_->fill();
  }

  if (jetanalyzer_) {
      jetanalyzer_->clear();
      jetanalyzer_->setZ(lTKCands,z0);
      jetanalyzer_->setGenJets(genParticles,1);
      jetanalyzer_->setJets(lCands,0);
      jetanalyzer_->setJets(lTKCands,2);
      jetanalyzer_->setJets(lTKVtxCands,3);
      jetanalyzer_->setJets(lRawCalo ,4);
      jetanalyzer_->setJets(lCorrCalo,5);
      jetanalyzer_->setJets(lPupCands,6);
      jetanalyzer_->fill();
  }
  
  if (isoanalyzer_) {
      isoanalyzer_->clear();
      isoanalyzer_->setGenMuons(genParticles,1);
      isoanalyzer_->matchMuons(lCands);
      isoanalyzer_->computeIso(lCands     , 0.3, "pf");
      isoanalyzer_->computeIso(lTKCands   , 0.3, "tk");
      isoanalyzer_->computeIso(lTKVtxCands, 0.3, "tkvtx");
      isoanalyzer_->computeIso(lPupCands  , 0.3, "pup");
      isoanalyzer_->fill();
  }

}

void NtupleProducer::addPF(const L1PFCollection &iCandidates, const std::string &iLabel, edm::Event& iEvent) { 
  std::unique_ptr<PFOutputCollection > corrCandidates( new PFOutputCollection );
  for(unsigned int i0 = 0; i0 < iCandidates.size(); i0++) { 
    reco::PFCandidate::ParticleType id = reco::PFCandidate::ParticleType::X; 
    int pCharge=0; 
    if(iCandidates[i0].pdgId() == l1tpf::Particle::CH) id = reco::PFCandidate::ParticleType::h;
    if(iCandidates[i0].pdgId() == l1tpf::Particle::EL) id = reco::PFCandidate::ParticleType::e;
    if(iCandidates[i0].pdgId() == l1tpf::Particle::NH) id = reco::PFCandidate::ParticleType::h0;
    if(iCandidates[i0].pdgId() == l1tpf::Particle::GAMMA) id = reco::PFCandidate::ParticleType::gamma;
    if(iCandidates[i0].pdgId() == l1tpf::Particle::MU) id = reco::PFCandidate::ParticleType::mu;
    if(iCandidates[i0].pdgId() == l1tpf::Particle::CH || iCandidates[i0].pdgId() == l1tpf::Particle::EL)  pCharge = 1;
    if(iCandidates[i0].pdgId() == l1tpf::Particle::MU)  pCharge = iCandidates[i0].charge();
    if (pCharge == 0 && (id !=  reco::PFCandidate::ParticleType::h0 && id != reco::PFCandidate::ParticleType::gamma)) {
        std::cout << "ERROR for " << iLabel << " candidate id " << iCandidates[i0].pdgId()  << ", pt " << iCandidates[i0].pt()  << ", eta " << iCandidates[i0].eta() << " has charge zero" << std::endl;
    }
    reco::PFCandidate pCand(pCharge,iCandidates[i0].p4(),id);
    pCand.setStatus(iCandidates[i0].status());
    if (pCharge != 0) {
        pCand.setVertex(reco::Particle::Point(0,0,iCandidates[i0].dz()));
    }
    pCand.set_mva_Isolated(iCandidates[i0].puppiWeight());
    corrCandidates->push_back(pCand);
  }
  //Fill!
  iEvent.put(std::move(corrCandidates),iLabel);
}
void NtupleProducer::addUInt(unsigned int value,std::string iLabel,edm::Event& iEvent) { 
  iEvent.put(std::make_unique<unsigned>(value), iLabel);
}

//////////////////////////// ------------------------------------------------------
//////////////////////////// ------------------------------------------------------

//--- Gen Matching
void NtupleProducer::genMatch(std::vector<double> &iGenVars,int iType,double iEta,double iPhi,double iPt,const reco::GenParticleCollection &iGenParticles) { 
  int lId = -999;
  double lPt = -1;
  double lDeltaRMin = 100;
  TLorentzVector lVec; lVec.SetPtEtaPhiM(0,0,0,0);
  for (reco::GenParticleCollection::const_iterator itGenP = iGenParticles.begin(); itGenP!=iGenParticles.end(); ++itGenP) {
    if(iType == 0 && itGenP->charge() == 0) continue;
    double deltaEta = itGenP->eta()-iEta;
    double deltaPhi = fabs(itGenP->phi()-iPhi); if(deltaPhi > 2.*TMath::Pi()-deltaPhi) deltaPhi = 2.*TMath::Pi()-deltaPhi;
    double deltaR   = sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);
    if(deltaR > 0.2) continue;
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
  if (fOutputName.empty() && !fTrackTree) return;
  //
  // Create output file, trees, and histograms
  //
  if (!fOutputName.empty()) {
      fOutputFile = new TFile(fOutputName.c_str(), "RECREATE");
      fTotalEvents = new TH1D("TotalEvents","TotalEvents",1,-10,10);
      fTrkInfoTree = new TTree("TrkInfo",   "TrkInfo");
  } else {
      edm::Service<TFileService> fs;
      fTotalEvents = fs->make<TH1D>("TotalEvents","TotalEvents",1,-10,10);
      fTrkInfoTree = fs->make<TTree>("TrkInfo",   "TrkInfo");
  }

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
  fTrkInfoTree->Branch("trkdz",   &trkdz,   "trkdz/F");
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
  if (fOutputName.empty()) return;
  //
  // Save to ROOT file
  //
  metanalyzer_->write();
  jetanalyzer_->write();
  isoanalyzer_->write();
  
  fOutputFile->cd();
  fTotalEvents->Write();
  fOutputFile->Write();
  fOutputFile->Close();

}

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
