#ifndef FASTPUPPI_NTUPLERPRODUCER_COMBINER_HH
#define FASTPUPPI_NTUPLERPRODUCER_COMBINER_HH
#include "TGraph.h"
#include <algorithm>  
#include <iostream>
#include "TFile.h"
#include "TTree.h"

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

class combiner { 
public:
  struct Particle { 
    Particle(double iEt,double iEta,double iPhi,double iM,int iId,double iSigma,double iDZ,double iCaloEta=0,double iCaloPhi=0, double iCharge = 0, double iQuality = -999, double iIsPV = 0, float alphaF = -999, float alphaC = -999, float puppiWeight = -99) 
    {
      Et=iEt; Eta=iEta; Phi=iPhi; M=iM; id=iId; sigma=iSigma; dZ=iDZ; caloEta=iCaloEta; caloPhi=iCaloPhi; charge = iCharge; quality = iQuality; isPV = iIsPV;
    }
    double Et;
    double Eta;
    double Phi;
    double M;
    double dZ;
    double sigma;
    double caloEta;
    double caloPhi;
    double charge;
    double quality;
    int id;
    int isPV; 

    float alphaF;
    float alphaC;
    float puppiWeight;

  };

  combiner(const std::string iPionFile,const std::string iElectronFile,const std::string iTrackFile,std::string iFile);
  void addCalo(double iCalo,double iEcal,double iCaloEta,double iCaloPhi,double iEcalEta,double iEcalPhi);
  void loadFile(TGraph** &iF1, std::string iFile);
  double correct(double iHcal,double iEcal,int iEta);
  void addTrack(double iPt,double iEta,double iPhi,double idZ,double iCaloEta,double iCaloPhi, double iCharge,int iQuality);
  void addMuon(double iPt, double iEta, double iPhi, double charge, double quality);
  void link();
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
private:
  void insert(Particle &iPartcle,std::vector<Particle> &iParticles);
  inline double  getTrkRes (double iPt,double iEta,double iPhi) { return fTrackRes   [translateIEta(iEta)]->Eval(iPt);}
  inline double  getPionRes(double iPt,double iEta,double iPhi) { return fElectronRes[translateIEta(iEta)]->Eval(iPt);}
  inline double  getEleRes (double iPt,double iEta,double iPhi) { return fPionRes    [translateIEta(iEta)]->Eval(iPt);}
  inline int     translateIEta(double iEta) { return int(10*std::max(std::min(iEta,3.0),-3.0))+30;}  
  double deltaR(Particle &iParticle1,Particle &iParticle2);
  double deltaRraw(Particle &iParticle1,Particle &iParticle2);
  TGraph **fTrackRes;
  TGraph **fElectronRes;
  TGraph **fPionRes;
  int   fNEta;
  double fDRMatch;
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
