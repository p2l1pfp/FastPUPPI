#include "FastPUPPI/NtupleProducer/interface/PuppiContainer.h"
#include "Math/ProbFunc.h"
#include "TMath.h"
#include <iostream>
#include <math.h>

using namespace std;

PuppiContainer::PuppiContainer(const edm::ParameterSet &iConfig) {
    fPuppiDiagnostics = iConfig.getParameter<bool>("puppiDiagnostics");
    fApplyCHS        = iConfig.getParameter<bool>("applyCHS");
    fInvert          = iConfig.getParameter<bool>("invertPuppi");    
    fUseExp          = iConfig.getParameter<bool>("useExp");
    fPuppiWeightCut  = iConfig.getParameter<double>("MinPuppiWeight");
    std::vector<edm::ParameterSet> lAlgos = iConfig.getParameter<std::vector<edm::ParameterSet> >("algos");
    fNAlgos = lAlgos.size();
    for(unsigned int i0 = 0; i0 < lAlgos.size(); i0++) {
        PuppiAlgo pPuppiConfig(lAlgos[i0]);
        fPuppiAlgo.push_back(pPuppiConfig);
    }
}

void PuppiContainer::initialize(const std::vector<combiner::Particle> &iRecoObjects) {
    //Clear everything
    fRecoParticles.resize(0);
    fPFParticles  .resize(0);
    fChargedPV    .resize(0);
    fPupParticles .resize(0);
    fMuons        .resize(0);
    fWeights      .resize(0);
    fVals.resize(0);
    fRawAlphas.resize(0);
    fAlphaMed     .resize(0);
    fAlphaRMS     .resize(0);
    //Link to the RecoObjects
    fPVFrac = 0.;
    fNPV    = 1.;
    fRecoParticles = iRecoObjects;
    for (unsigned int i = 0; i < fRecoParticles.size(); i++){
      combiner::Particle fRecoParticle = fRecoParticles[i];
      if(fRecoParticle.id == 4) fMuons.push_back(fRecoParticle);
      if(fRecoParticle.id == 4) continue;
      fPFParticles.push_back(fRecoParticle);
      if(std::abs(fRecoParticle.pvid) == 1) fChargedPV.push_back(fRecoParticle);
      if(std::abs(fRecoParticle.pvid) >= 1 ) fPVFrac+=1.;
    }
    if (fPVFrac != 0) fPVFrac = double(fChargedPV.size())/fPVFrac;
    else fPVFrac = 0;
}
PuppiContainer::~PuppiContainer(){}
double PuppiContainer::goodVar(const combiner::Particle  &iPart,const std::vector<combiner::Particle> &iParts, int iOpt,double iRCone) {
    double lPup = 0;
    lPup = var_within_R(iOpt,iParts,iPart,iRCone);
    return lPup;
}
void PuppiContainer::sel(double iR,vector<combiner::Particle> &iFill,const vector<combiner::Particle> iParticles,const combiner::Particle &iPart) { 
  for(unsigned int i0 = 0; i0 < iParticles.size(); i0++) { 
    double pDEta = iParticles[i0].Eta-iPart.Eta;
    double pDPhi = fabs(iParticles[i0].Phi-iPart.Phi);
    if(pDPhi > 2.*TMath::Pi()-pDPhi) pDPhi =  2.*TMath::Pi()-fabs(iParticles[i0].Phi-iPart.Phi);
    if(pDEta*pDEta+pDPhi*pDPhi < iR*iR) iFill.push_back(iParticles[i0]); 
  }
}
double PuppiContainer::var_within_R(int iId,const vector<combiner::Particle> & particles,const combiner::Particle& centre, double R){
    if(iId == -1) return 1;
    vector<combiner::Particle> near_particles; 
    sel(R,near_particles,particles,centre);
    double var = 0;
    for(unsigned int i=0; i<near_particles.size(); i++){
        double pDEta = near_particles[i].Eta-centre.Eta;
        double pDPhi = std::abs(near_particles[i].Phi-centre.Phi);
        if(pDPhi > 2.*M_PI-pDPhi) pDPhi =  2.*M_PI-pDPhi;
        double pDR2 = pDEta*pDEta+pDPhi*pDPhi;
        if(std::abs(pDR2)  <  0.0001) continue;
        if(iId == 0) var += (near_particles[i].Et/pDR2);
        if(iId == 1) var += near_particles[i].Et;
        if(iId == 2) var += (1./pDR2);
        if(iId == 3) var += (1./pDR2);
        if(iId == 4) var += near_particles[i].Et;
        if(iId == 5) var += (near_particles[i].Et * near_particles[i].Et/pDR2);
    }
    if(iId == 1) var += centre.Et; //Sum in a cone
    if(iId == 0 && var != 0) var = log(var);
    if(iId == 3 && var != 0) var = log(var);
    if(iId == 5 && var != 0) var = log(var);
    return var;
}
//In fact takes the median not the average
void PuppiContainer::getRMSAvg(int iOpt,std::vector<combiner::Particle> const &iConstits,std::vector<combiner::Particle> const &iParticles,std::vector<combiner::Particle> const &iChargedParticles) {
    for(unsigned int i0 = 0; i0 < iConstits.size(); i0++ ) {
        double pVal = -1;
        //Calculate the Puppi Algo to use
        int  pPupId   = getPuppiId(iConstits[i0].Et,iConstits[i0].Eta);
        if(pPupId == -1 || fPuppiAlgo[pPupId].numAlgos() <= iOpt){
            fVals.push_back(-1);
            continue;
        }
        //Get the Puppi Sub Algo (given iteration)
        int  pAlgo    = fPuppiAlgo[pPupId].algoId   (iOpt);
        bool pCharged = fPuppiAlgo[pPupId].isCharged(iOpt);
        double pCone  = fPuppiAlgo[pPupId].coneSize (iOpt);
        //Compute the Puppi Metric
        if(!pCharged) pVal = goodVar(iConstits[i0],iParticles       ,pAlgo,pCone);
        if( pCharged) pVal = goodVar(iConstits[i0],iChargedParticles,pAlgo,pCone);
        fVals.push_back(pVal);

        fPuppiAlgo[pPupId].add(iConstits[i0],pVal,iOpt);
        /*
	//code added by Nhan, now instead for every algorithm give it all the particles
        for(int i1 = 0; i1 < fNAlgos; i1++){
            pAlgo    = fPuppiAlgo[i1].algoId   (iOpt);
            pCharged = fPuppiAlgo[i1].isCharged(iOpt);
            pCone    = fPuppiAlgo[i1].coneSize (iOpt);
            double curVal = -1; 
            if(!pCharged) curVal = goodVar(iConstits[i0],iParticles       ,pAlgo,pCone);
            if( pCharged) curVal = goodVar(iConstits[i0],iChargedParticles,pAlgo,pCone);
            //std::cout << "i1 = " << i1 << ", curVal = " << curVal << ", eta = " << iConstits[i0].eta() << ", pupID = " << pPupId << std::endl;
            fPuppiAlgo[i1].add(iConstits[i0],curVal,iOpt);
        }
	*/
    }
    for(int i0 = 0; i0 < fNAlgos; i0++) fPuppiAlgo[i0].computeMedRMS(iOpt,fPVFrac);
}
//In fact takes the median not the average
void PuppiContainer::getRawAlphas(int iOpt,std::vector<combiner::Particle> const &iConstits,std::vector<combiner::Particle> const &iParticles,std::vector<combiner::Particle> const &iChargedParticles) {
    for(int j0 = 0; j0 < fNAlgos; j0++){
        for(unsigned int i0 = 0; i0 < iConstits.size(); i0++ ) {
            double pVal = -1;
            //Get the Puppi Sub Algo (given iteration)
            int  pAlgo    = fPuppiAlgo[j0].algoId   (iOpt);
            bool pCharged = fPuppiAlgo[j0].isCharged(iOpt);
            double pCone  = fPuppiAlgo[j0].coneSize (iOpt);
            //Compute the Puppi Metric
            if(!pCharged) pVal = goodVar(iConstits[i0],iParticles       ,pAlgo,pCone);
            if( pCharged) pVal = goodVar(iConstits[i0],iChargedParticles,pAlgo,pCone);
            fRawAlphas.push_back(pVal);
        }
    }
}
int    PuppiContainer::getPuppiId( float iPt, float iEta) {
    int lId = -1;
    for(int i0 = 0; i0 < fNAlgos; i0++) {
        int nEtaBinsPerAlgo = fPuppiAlgo[i0].etaBins();
        for (int i1 = 0; i1 < nEtaBinsPerAlgo; i1++){
            if ( (std::abs(iEta) > fPuppiAlgo[i0].etaMin(i1)) && (std::abs(iEta) < fPuppiAlgo[i0].etaMax(i1)) ){ 
                fPuppiAlgo[i0].fixAlgoEtaBin( i1 );
                if(iPt > fPuppiAlgo[i0].ptMin()){
                    lId = i0; 
                    break;
                }
            }
        }
    }
    return lId;
}
std::vector<double> const & PuppiContainer::puppiFetch() {
    fPupParticles .resize(0);
    fWeights      .resize(0);
    fVals         .resize(0);
    for(int i0 = 0; i0 < fNAlgos; i0++) fPuppiAlgo[i0].reset();
    
    int lNMaxAlgo = 1;
    for(int i0 = 0; i0 < fNAlgos; i0++) lNMaxAlgo = std::max(fPuppiAlgo[i0].numAlgos(),lNMaxAlgo);
    //Run through all compute mean and RMS
    int lNParticles    = fPFParticles.size();
    for(int i0 = 0; i0 < lNMaxAlgo; i0++) {
        getRMSAvg(i0,fPFParticles,fPFParticles,fChargedPV);
    }
    if (fPuppiDiagnostics) getRawAlphas(0,fPFParticles,fPFParticles,fChargedPV);

    std::vector<double> pVals;
    for(int i0 = 0; i0 < lNParticles; i0++) {
        //Refresh
        pVals.clear();
        double pWeight = 1;
        //Get the Puppi Id and if ill defined move on
        int  pPupId   = getPuppiId(fPFParticles[i0].Et,fPFParticles[i0].Eta);
        if(pPupId == -1) {
            fWeights .push_back(pWeight);
            fAlphaMed.push_back(-10);
            fAlphaRMS.push_back(-10);            
            continue;
        }
        // fill the p-values
        double pChi2   = 0;
        //Fill and compute the PuppiWeight
        int lNAlgos = fPuppiAlgo[pPupId].numAlgos();
        for(int i1 = 0; i1 < lNAlgos; i1++) pVals.push_back(fVals[lNParticles*i1+i0]);
        pWeight = fPuppiAlgo[pPupId].compute(pVals,pChi2);
        //Apply the CHS weights
        if(fPFParticles[i0].id == 1 && fApplyCHS ) pWeight = 1;
        if(fPFParticles[i0].id == 2 && fApplyCHS ) pWeight = 0;
        //Basic Cuts
        if(pWeight                         < fPuppiWeightCut) pWeight = 0;  //==> Elminate the low Weight stuff
        if(pWeight*fPFParticles[i0].Et     < fPuppiAlgo[pPupId].neutralPt(fNPV) && fPFParticles[i0].id == 0 ) pWeight = 0;  //threshold cut on the neutral Pt
        if(fInvert) pWeight = 1.-pWeight;
        fWeights .push_back(pWeight);
        fAlphaMed.push_back(fPuppiAlgo[pPupId].median());
        fAlphaRMS.push_back(fPuppiAlgo[pPupId].rms());        
        //Produce
	combiner::Particle curjet(pWeight*fPFParticles[i0].Et,fPFParticles[i0].Eta,fPFParticles[i0].Phi,fPFParticles[i0].M,fPFParticles[i0].id,
				  sqrt(pWeight)*fPFParticles[i0].sigma,fPFParticles[i0].dZ,fPFParticles[i0].pvid,
				  fPFParticles[i0].caloEta,fPFParticles[i0].caloPhi,fPFParticles[i0].charge,fPFParticles[i0].quality); 
        fPupParticles.push_back(curjet);
    }
    for(unsigned int i0 = 0; i0 < fMuons.size(); i0++) fPupParticles.push_back(fMuons[i0]);
    return fWeights;
}


