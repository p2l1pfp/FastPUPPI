#ifndef FASTPUPPI_NTUPLERPRODUCER_COMBINER_HH
#define FASTPUPPI_NTUPLERPRODUCER_COMBINER_HH
#include "TGraph.h"
#include <algorithm>  
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "L1TPFParticle.h"
#include "L1TPFUtils.h"

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

class combiner { 
public:
  enum MyParticleId { CH=0, EL=1, NH=2, GAMMA=3, MU=4 };
  typedef l1tpf::Particle Particle;

  combiner(const std::string &iPionFile,const std::string & iElectronFile,const std::string &iTrackFile,const std::string &iFile,double iEtaCharged,double iPuppiPt,double iVtxRes);
  l1tpf::Particle makeCalo(double iCalo,double iEcal,double iCaloEta,double iCaloPhi,double iEcalEta,double iEcalPhi) const ;
  void addCalo(const l1tpf::Particle & particle);
  void addMuon(const l1tpf::Particle & particle); 
  void loadFile(TGraph** &iF1, std::string iFile);
  double correct(double iHcal,double iEcal,int iEta);
  void addTrack(const l1tpf::Particle & particle); 
  void link(bool iMetRate);
  void doVertexing();  
  void fetchPuppi();
  void fill();
  void write(){fFile->cd(); fTree->Write(); fFile->Write(); fFile->Close();}
  void merge(Particle &iTkParticle,Particle &iParticle1,std::vector<Particle> &iCollection);
  inline void    clear() {fTkParticles.clear(); fParticles.clear(); fMuParticles.clear(); fTkParticlesWVertexing.clear(); fParticlesPuppi.clear(); }
  inline std::vector<Particle> candidates()   { return fParticles;}
  inline std::vector<Particle> tkcandidates() { return fTkParticles;}
  inline std::vector<Particle> tkvtxcandidates() { return fTkParticlesWVertexing;}
  inline std::vector<Particle> mucandidates() { return fMuParticles;}
  inline std::vector<Particle> puppiFetch() { return fParticlesPuppi; }
  inline double dZ() { return fDZ;} 
  inline std::pair<float,float> alphaCMedRms() { return std::make_pair(alphaCMed,alphaCRms); } 
  inline std::pair<float,float> alphaFMedRms() { return std::make_pair(alphaFMed,alphaFRms); } 
  inline double  getTrkRes (double iPt,double iEta,double iPhi) const { return fTrackRes   [l1tpf::translateAEta(l1tpf::translateIEta(iEta))]->Eval(iPt);}
  inline double  getEleRes (double iPt,double iEta,double iPhi) const { return fElectronRes[l1tpf::translateAEta(l1tpf::translateIEta(iEta))]->Eval(iPt);}
  inline double  getPionRes(double iPt,double iEta,double iPhi) const { 
    //double lPt30 = fPionRes    [l1tpf::translateAEta(l1tpf::translateIEta(iEta))]->Eval(30.);
    double lPt   = fPionRes    [l1tpf::translateAEta(l1tpf::translateIEta(iEta))]->Eval(iPt);
    return lPt*0.5;}//(sqrt(lPt*lPt+lPt30*lPt30));}
private:
  void insert(const Particle &iPartcle,std::vector<Particle> &iParticles);
  inline int     translateIEtaOld(double iEta) { return int(10*std::max(std::min(iEta,3.0),-3.0))+30;}  
  double deltaR(Particle &iParticle1,Particle &iParticle2);
  double deltaRraw(Particle &iParticle1,Particle &iParticle2);
  TGraph **fTrackRes;
  TGraph **fElectronRes;
  TGraph **fPionRes;
  int   fNEta;
  double fEta;
  double fDRMatch;
  double fPuppiPt;
  double fVtxRes;
  std::vector<Particle> fTkParticles;
  std::vector<Particle> fTkParticlesWVertexing;  
  std::vector<Particle> fMuParticles;
  std::vector<Particle> fParticles;
  std::vector<Particle> fParticlesPuppi;
  TFile  *fFile;
  TTree  *fTree;
  double  fDZ;

  void computeAlphas(std::vector<Particle> neighborParticles, int isCentral);
  void computeMedRMS();
  void computeWeights();
  float alphaFMed;
  float alphaFRms;
  float alphaCMed;
  float alphaCRms;
  std::vector<float> b_pt;
  std::vector<float> b_alphaF;
  std::vector<float> b_alphaC;
  std::vector<float> b_Eta;
  std::vector<float> b_Phi;
  std::vector<float> b_Et;
  std::vector<float> b_PuppiWeight;
  int b_nParticles; 
  float b_alphaFMed;
  float b_alphaFRms;
  float b_alphaCMed;
  float b_alphaCRms;



};
#endif
