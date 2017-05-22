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

#include "FastPUPPI/NtupleProducer/interface/corrector.hh"
#include "FastPUPPI/NtupleProducer/interface/combiner.hh"
#include "FastPUPPI/NtupleProducer/interface/metanalyzer.hh"
#include "FastPUPPI/NtupleProducer/interface/jetanalyzer.hh"
#include "FastPUPPI/NtupleProducer/interface/isoanalyzer.hh"
#include "FastPUPPI/NtupleProducer/interface/SimpleCalibrations.h"
#include "FastPUPPI/NtupleProducer/interface/DiscretePF.h"

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
  const edm::InputTag MuonTPTag_;
  const edm::InputTag GenParTag_;
  const std::string CorrectorTag_;
  const std::string Corrector2Tag_;
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
  void addPF(const std::vector<combiner::Particle> &iCandidates, const std::string &iLabel, edm::Event& iEvent);
  void addUInt(unsigned int value,std::string iLabel,edm::Event& iEvent);

  edm::EDGetTokenT<reco::GenParticleCollection>   TokGenPar_;
  edm::EDGetTokenT<L1PFCollection>                TokL1TrackTPTag_;
  std::vector<edm::EDGetTokenT<L1PFCollection>>   TokCaloClusterTags_, TokEmClusterTags_;
  edm::EDGetTokenT<L1PFCollection>                TokMuonTPTag_;
  double trkPt_;
  bool   metRate_;
  double etaCharged_;
  double puppiPtCut_, puppiDr_;
  double vtxRes_;
  corrector* corrector_;
  corrector* ecorrector_;
  combiner * connector_;
  combiner * rawconnector_;
  metanalyzer* metanalyzer_;
  jetanalyzer* jetanalyzer_;
  isoanalyzer* isoanalyzer_;
  // simplified resolutions
  l1tpf::SimpleResol simpleResolEm_, simpleResolHad_, simpleResolTrk_;
  // discretized version
  l1tpf_int::RegionMapper l1regions_;
  l1tpf_int::PFAlgo       l1pfalgo_;
  // debug flag
  int fDebug;
     
  // declare variables for output file
  std::string           fOutputName;
  TFile                 *fOutputFile;
  TH1D                  *fTotalEvents;

  TTree                 *fTrkInfoTree;
  TTree                 *fEcalInfoTree;
  TTree                 *fHcalInfoTree;
  float runNum, lumiSec, evtNum;
  float trkNum;
  float trkPx, trkPz, trkPy, trkPt, trkEta, trkPhi, trkz0, trkd0;    
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
  MuonTPTag_            (iConfig.getParameter<edm::InputTag>("MuonTPTag")),
  GenParTag_            (iConfig.getParameter<edm::InputTag>("genParTag")),
  CorrectorTag_         (getFilePath(iConfig,"corrector")),
  ECorrectorTag_        (getFilePath(iConfig,"ecorrector")),
  TrackResTag_          (getFilePath(iConfig,"trackres")),
  EleResTag_            (getFilePath(iConfig,"eleres")),
  PionResTag_           (getFilePath(iConfig,"pionres")),
  trkPt_                (iConfig.getParameter<double>       ("trkPtCut")),
  metRate_              (iConfig.getParameter<bool>         ("metRate")),
  etaCharged_           (iConfig.getParameter<double>       ("etaCharged")),
  puppiPtCut_           (iConfig.getParameter<double>       ("puppiPtCut")),
  puppiDr_              (iConfig.getParameter<double>       ("puppiDr")),
  vtxRes_               (iConfig.getParameter<double>       ("vtxRes")),
  l1regions_            (iConfig),
  l1pfalgo_             (iConfig),
  fDebug                (iConfig.getUntrackedParameter<int>("debug",0)),
  fOutputName           (iConfig.getUntrackedParameter<std::string>("outputName", "ntuple.root")),
  fOutputFile           (0),
  fTotalEvents          (0),
  fTrkInfoTree          (0)
{
  //now do what ever other initialization is needed
  corrector_  = new corrector(CorrectorTag_,11,fDebug);
  ecorrector_ = new corrector(ECorrectorTag_,1,fDebug);
  connector_  = new combiner (PionResTag_,EleResTag_,TrackResTag_,!fOutputName.empty() ? "puppi.root" : "",etaCharged_,puppiPtCut_,puppiDr_,vtxRes_,fDebug);
  rawconnector_  = new combiner (PionResTag_,EleResTag_,TrackResTag_,!fOutputName.empty() ? "puppiraw.root" : "",etaCharged_,puppiPtCut_,puppiDr_,vtxRes_);
  simpleResolEm_ = l1tpf::SimpleResol(iConfig,"simpleResolEm");
  simpleResolHad_ = l1tpf::SimpleResol(iConfig,"simpleResolHad");
  simpleResolTrk_ = l1tpf::SimpleResol(iConfig,"simpleResolTrk");
  if (iConfig.existsAs<edm::ParameterSet>("linking")) {
    edm::ParameterSet linkcfg = iConfig.getParameter<edm::ParameterSet>("linking");
    connector_->configureLinking(linkcfg.getParameter<double>("trackCaloDR"),
                                linkcfg.getParameter<double>("trackCaloNSigmaLow"),
                                linkcfg.getParameter<double>("trackCaloNSigmaHigh"),
                                linkcfg.getParameter<bool>("useTrackCaloSigma"),
                                linkcfg.getParameter<bool>("rescaleUnmatchedTrack"),
                                linkcfg.getParameter<double>("maxInvisiblePt"));
  }
  if (fOutputName.empty()) {
      metanalyzer_ = nullptr;
      jetanalyzer_ = nullptr;
      isoanalyzer_ = nullptr;
  } else {
      metanalyzer_ = new metanalyzer("MetFile.root");
      jetanalyzer_ = new jetanalyzer("JetFile.root");
      isoanalyzer_ = new isoanalyzer("IsoFile.root");
  }
  produces<PFOutputCollection>("TK");
  produces<PFOutputCollection>("TKVtx");
  produces<PFOutputCollection>("RawCalo");
  produces<PFOutputCollection>("Calo");
  produces<PFOutputCollection>("PF");
  produces<PFOutputCollection>("Puppi");
  produces<PFOutputCollection>("L1TK");
  produces<PFOutputCollection>("L1Calo");
  produces<PFOutputCollection>("L1PF");
  produces<PFOutputCollection>("L1Puppi");
  TokGenPar_       = consumes<reco::GenParticleCollection>( GenParTag_    );
  TokL1TrackTPTag_ = consumes<L1PFCollection>( L1TrackTag_  );
  for (const edm::InputTag &tag : CaloClusterTags_) {
    TokCaloClusterTags_.push_back(consumes<L1PFCollection>(tag));
  }
  for (const edm::InputTag &tag : EmClusterTags_) {
    TokEmClusterTags_.push_back(consumes<L1PFCollection>(tag));
  }
  TokMuonTPTag_    = consumes<L1PFCollection>( MuonTPTag_  );
  produces<unsigned int>("totNL1TK");
  produces<unsigned int>("totNL1Mu");
  produces<unsigned int>("totNL1Calo");
  produces<unsigned int>("totNL1PF");
  produces<unsigned int>("totNL1Puppi");
  produces<unsigned int>("maxNL1TK");
  produces<unsigned int>("maxNL1Mu");
  produces<unsigned int>("maxNL1Calo");
  produces<unsigned int>("maxNL1PF");
  produces<unsigned int>("maxNL1Puppi");
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
NtupleProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //fOutputFile->cd();
  if (!fOutputName.empty()) fTotalEvents->Fill(1);  
  using namespace edm;

  edm::Handle<reco::GenParticleCollection> hGenParProduct;
  iEvent.getByToken(TokGenPar_,hGenParProduct);
  const reco::GenParticleCollection genParticles = *(hGenParProduct.product());  

  connector_->clear();
  rawconnector_->clear();
  l1regions_.clear(); 
  /// ----------------TRACK INFO-------------------
  edm::Handle<std::vector<l1tpf::Particle>> l1tks;
  iEvent.getByLabel(L1TrackTag_, l1tks);
  trkNum = 0;
  for (l1tpf::Particle tk : *l1tks) { // no const & since we modify it to set the sigma
      // the combiner knows the sigma, the track producer doesn't
      if (simpleResolTrk_.empty()) {
          tk.setSigma(connector_->getTrkRes(tk.pt(), tk.eta(), tk.phi())); 
      } else {
          tk.setSigma(simpleResolTrk_(tk.pt(), std::abs(tk.eta())));
      }
      if (simpleResolHad_.empty()) {
          tk.setCaloSigma(connector_->getPionRes(tk.pt(), tk.eta(), tk.phi()));
      } else {
          tk.setCaloSigma(simpleResolHad_(tk.pt(), std::abs(tk.eta())));
      }
      if (fDebug > 1) printf("tk %7.2f eta %+4.2f has sigma %4.2f  sigma/pt %5.3f, calo sigma/pt %5.3f\n", tk.pt(), tk.eta(), tk.sigma(), tk.sigma()/tk.pt(), tk.caloSigma()/tk.pt());
      // adding objects to PF
      if(tk.pt() > trkPt_) l1regions_.addTrack(tk);      
      if(tk.pt() > trkPt_) connector_->addTrack(tk);      
      if(tk.pt() > trkPt_) rawconnector_->addTrack(tk);
      /// filling the tree    
      if (fOutputName.empty()) continue;
      trkPx  = tk.px();
      trkPz  = tk.py();
      trkPy  = tk.pz();
      trkPt  = tk.pt();
      trkEta = tk.eta();
      trkPhi = tk.phi();
      trkz0  = tk.vertex().Z();
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
  edm::Handle<std::vector<l1tpf::Particle>> l1mus;
  iEvent.getByToken(TokMuonTPTag_, l1mus);
  for (l1tpf::Particle mu : *l1mus) { // no const & since we modify it to set the sigma
      mu.setSigma(connector_->getTrkRes(mu.pt(), mu.eta(), mu.phi()));  // this needs to be updated with the muon resolutions!      
      connector_->addMuon(mu);
      l1regions_.addMuon(mu);
  }


  /// ----------------Clustered Calo INFO-------------------
  L1PFCollection calos, emcalos;
  edm::Handle<L1PFCollection> hcals;
  for (const auto & token : TokCaloClusterTags_) {
      iEvent.getByToken(token, hcals);
      calos.insert(calos.end(), hcals->begin(), hcals->end());
  }
  for (const auto & token : TokEmClusterTags_) {
      iEvent.getByToken(token, hcals);
      emcalos.insert(emcalos.end(), hcals->begin(), hcals->end());
  }


  // add uncalibrated (or at least un-recalibrated) calos to rawconnector
  for (const l1tpf::Particle & calo : calos) {
      rawconnector_->addCalo(calo);
  }

  // ccalibrate and do the rest
  for (l1tpf::Particle &calo : calos) {
      if (correctCaloEnergies_) {
          float ptcorr = corrector_->correct(calo.pt(), calo.emEt(), calo.iEta(), calo.iPhi());
          calo.setPt(ptcorr);
      }
      if (simpleResolEm_.empty() || simpleResolHad_.empty()) {
          calo.setSigma(calo.pdgId() == combiner::Particle::GAMMA ?
                  connector_->getEleRes(calo.pt(), calo.eta(), calo.phi()) :
                  connector_->getPionRes(calo.pt(), calo.eta(), calo.phi()));
      } else {
          calo.setSigma(calo.pdgId() == combiner::Particle::GAMMA ?
                  simpleResolEm_(calo.pt(), std::abs(calo.eta())) :
                  simpleResolHad_(calo.pt(), std::abs(calo.eta())) );
          if (fDebug > 1) printf("calo %7.2f eta %+4.2f has sigma %4.2f  sigma/pt %5.3f\n", calo.pt(), calo.eta(), calo.sigma(), calo.sigma()/calo.pt());
      }
  }

  // pass to the PF algo
  for (const l1tpf::Particle & calo : calos) {
      connector_->addCalo(calo);
      l1regions_.addCalo(calo); 
  }
  for (const l1tpf::Particle & calo : emcalos) {
      l1regions_.addEmCalo(calo); 
  }


  std::vector<combiner::Particle> lRawCalo      = rawconnector_->candidates();
  std::vector<combiner::Particle> lCorrCalo     = connector_->candidates();
  connector_->link(metRate_);
  std::vector<combiner::Particle> lCands        = connector_->candidates();
  std::vector<combiner::Particle> lTKCands      = connector_->tkcandidates();
  connector_->doVertexing();
  std::vector<combiner::Particle> lTKVtxCands   = connector_->tkvtxcandidates();
  connector_->fetchPuppi();
  connector_->fill();
  std::vector<combiner::Particle> lPupCands     = connector_->puppiFetch();

  // Now we run the discretized version
  // First, get a copy of the discretized inputs (for reference later)
  std::vector<combiner::Particle> ll1CaloCands = l1regions_.fetchCalo();
  std::vector<combiner::Particle> ll1TkCands   = l1regions_.fetchTracks();
  // then get global inputs
  float z0 = connector_->dZ();
  std::pair<float,float> alphaC = connector_->alphaCMedRms(), alphaF = connector_->alphaFMedRms();
  if (fDebug) std::cout << " z0 = " << z0 << "\t alphaC = " << alphaC.first << " +/- " << alphaC.second << "\t alphaF = " << alphaF.first << " +/- " << alphaF.second << std::endl;

  // run PF in each region
  for (auto & l1region : l1regions_.regions()) {
      l1pfalgo_.runPF(l1region);
      l1pfalgo_.runPuppi(l1region, z0, -1., alphaC.first, alphaC.second, alphaF.first, alphaF.second);
  }

  addPF(lRawCalo ,"RawCalo" ,iEvent);
  addPF(lCorrCalo,"Calo"    ,iEvent);
  addPF(lTKCands ,"TK"      ,iEvent);
  addPF(lTKVtxCands,"TKVtx" ,iEvent);
  addPF(lCands   ,"PF"      ,iEvent);
  addPF(lPupCands,"Puppi"   ,iEvent);

  std::vector<combiner::Particle> ll1PFCands   = l1regions_.fetch(false);
  std::vector<combiner::Particle> ll1PupCands  = l1regions_.fetch(true);
  addPF(ll1CaloCands,"L1Calo" ,iEvent);
  addPF(ll1TkCands,  "L1TK"   ,iEvent);
  addPF(ll1PFCands,  "L1PF"   ,iEvent);
  addPF(ll1PupCands, "L1Puppi",iEvent);
  unsigned int totNL1Calo = 0, totNL1TK = 0, totNL1Mu = 0, totNL1PF = 0, totNL1Puppi = 0;
  unsigned int maxNL1Calo = 0, maxNL1TK = 0, maxNL1Mu = 0, maxNL1PF = 0, maxNL1Puppi = 0;
  for (const auto & r : l1regions_.regions()) {
      totNL1Calo += r.calo.size();
      totNL1TK += r.track.size();
      totNL1PF += r.pf.size();
      totNL1Mu += r.muon.size();
      totNL1Puppi += r.puppi.size();
      maxNL1Calo = std::max<unsigned>( maxNL1Calo, r.calo.size() );
      maxNL1TK = std::max<unsigned>( maxNL1TK, r.track.size() );
      maxNL1PF = std::max<unsigned>( maxNL1PF, r.pf.size() );
      maxNL1Mu = std::max<unsigned>( maxNL1Mu, r.muon.size() );
      maxNL1Puppi = std::max<unsigned>( maxNL1Puppi, r.puppi.size() );
  }
  addUInt(totNL1Calo, "totNL1Calo", iEvent); addUInt(maxNL1Calo, "maxNL1Calo", iEvent);
  addUInt(totNL1TK, "totNL1TK", iEvent); addUInt(maxNL1TK, "maxNL1TK", iEvent);
  addUInt(totNL1Mu, "totNL1Mu", iEvent); addUInt(maxNL1Mu, "maxNL1Mu", iEvent);
  addUInt(totNL1PF, "totNL1PF", iEvent); addUInt(maxNL1PF, "maxNL1PF", iEvent);
  addUInt(totNL1Puppi, "totNL1Puppi", iEvent); addUInt(maxNL1Puppi, "maxNL1Puppi", iEvent);
  
  if (metanalyzer_) {
      metanalyzer_->clear();
      metanalyzer_->setZ(lTKCands,connector_->dZ());
      metanalyzer_->setMETRecoil(2,lTKCands ,false);
      metanalyzer_->setMETRecoil(0,lCands,false);
      metanalyzer_->setMETRecoil(3,lRawCalo ,true);
      metanalyzer_->setMETRecoil(1,lCorrCalo,true);
      metanalyzer_->setMETRecoil(5,lPupCands,false);
      metanalyzer_->setMETRecoil(4,lTKVtxCands ,false);
      metanalyzer_->setMETRecoil(6,ll1PFCands ,false);
      metanalyzer_->setMETRecoil(7,ll1PupCands ,false);      
      metanalyzer_->setGenMET(genParticles);
      metanalyzer_->fill();
  }

  if (jetanalyzer_) {
      jetanalyzer_->clear();
      jetanalyzer_->setZ(lTKCands,connector_->dZ());
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
    if(iCandidates[i0].pdgId() == combiner::Particle::CH) id = reco::PFCandidate::ParticleType::h;
    if(iCandidates[i0].pdgId() == combiner::Particle::EL) id = reco::PFCandidate::ParticleType::e;
    if(iCandidates[i0].pdgId() == combiner::Particle::NH) id = reco::PFCandidate::ParticleType::h0;
    if(iCandidates[i0].pdgId() == combiner::Particle::GAMMA) id = reco::PFCandidate::ParticleType::gamma;
    if(iCandidates[i0].pdgId() == combiner::Particle::MU) id = reco::PFCandidate::ParticleType::mu;
    if(iCandidates[i0].pdgId() == combiner::Particle::CH || iCandidates[i0].pdgId() == combiner::Particle::EL)  pCharge = 1;
    if(iCandidates[i0].pdgId() == combiner::Particle::MU)  pCharge = iCandidates[i0].charge();
    reco::PFCandidate pCand(pCharge,iCandidates[i0].p4(),id);
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
  if (fOutputName.empty()) return;
  //
  // Create output file, trees, and histograms
  //
  fOutputFile = new TFile(fOutputName.c_str(), "RECREATE");
  fTotalEvents = new TH1D("TotalEvents","TotalEvents",1,-10,10);
  fTrkInfoTree     = new TTree("TrkInfo",   "TrkInfo");

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
  
  connector_->write();
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
