#ifndef COMMONTOOLS_PUPPI_PUPPICONTAINER_H_
#define COMMONTOOLS_PUPPI_PUPPICONTAINER_H_
#include "FastPUPPI/NtupleProducer/interface/combiner.hh"
#include "FastPUPPI/NtupleProducer/interface/PuppiAlgo.h"

class PuppiContainer{
public:
    PuppiContainer(const edm::ParameterSet &iConfig);
    ~PuppiContainer(); 
    void initialize(const std::vector<combiner::Particle> &iRecoObjects);
    void setNPV(int iNPV){ fNPV = iNPV; }

    std::vector<combiner::Particle> const & pfParticles() const { return fPFParticles; }    
    std::vector<combiner::Particle> const & pvParticles() const { return fChargedPV; }        
    std::vector<double> const & puppiFetch();
    const std::vector<double> & puppiRawAlphas(){ return fRawAlphas; }
    const std::vector<double> & puppiAlphas(){ return fVals; }
    // const std::vector<double> puppiAlpha   () {return fAlpha;}
    const std::vector<double> & puppiAlphasMed() {return fAlphaMed;}
    const std::vector<double> & puppiAlphasRMS() {return fAlphaRMS;}

    int puppiNAlgos(){ return fNAlgos; }
    std::vector<combiner::Particle> const & puppiParticles() const { return fPupParticles;}

protected:
    double  goodVar      (const combiner::Particle &iPart,const std::vector<combiner::Particle> &iParts, int iOpt,double iRCone);
    void    sel(double iR,std::vector<combiner::Particle> &iFill,const std::vector<combiner::Particle> iParticles,const combiner::Particle &iPart);
    void    getRMSAvg    (int iOpt,std::vector<combiner::Particle> const &iConstits,std::vector<combiner::Particle> const &iParticles,std::vector<combiner::Particle> const &iChargeParticles);
    void    getRawAlphas    (int iOpt,std::vector<combiner::Particle> const &iConstits,std::vector<combiner::Particle> const &iParticles,std::vector<combiner::Particle> const &iChargeParticles);
    double  getChi2FromdZ(double iDZ);
    int     getPuppiId   ( float iPt, float iEta);
    double  var_within_R (int iId,const std::vector<combiner::Particle> & particles, const combiner::Particle& centre, double R);  
    
    bool      fPuppiDiagnostics;
    std::vector<combiner::Particle> fRecoParticles;
    std::vector<combiner::Particle> fPFParticles;
    std::vector<combiner::Particle> fChargedPV;
    std::vector<combiner::Particle> fPupParticles;
    std::vector<combiner::Particle> fMuons;
    std::vector<double>    fWeights;
    std::vector<double>    fVals;
    std::vector<double>    fRawAlphas;
    std::vector<double>    fAlphaMed;
    std::vector<double>    fAlphaRMS;

    bool   fApplyCHS;
    bool   fInvert;
    bool   fUseExp;
    double fNeutralMinPt;
    double fNeutralSlope;
    double fPuppiWeightCut;
    int    fNAlgos;
    int    fNPV;
    double fPVFrac;
    std::vector<PuppiAlgo> fPuppiAlgo;
};
#endif

