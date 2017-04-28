#ifndef FASTPUPPI_NTUPLERPRODUCER_L1TPFCandidate_HH
#define FASTPUPPI_NTUPLERPRODUCER_L1TPFCandidate_HH

#include "FastPUPPI/NtupleProducer/interface/L1TPFUtils.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include <TLorentzVector.h>

namespace l1tpf {
    class Particle : public reco::LeafCandidate {
        public:
            Particle() {}
            Particle(double iEt,double iEta,double iPhi,double iM,int iId,double iSigma=0,double iDZ=0,double iCaloEta=0,double iCaloPhi=0, double iCharge = 0, double iQuality = -999, double iIsPV = 0, float alphaF = -999, float alphaC = -999, float puppiWeight = -99) :
                LeafCandidate(iCharge, reco::LeafCandidate::PolarLorentzVector(iEt,iEta,iPhi,iM), reco::LeafCandidate::Point(), iId),
                dZ_(iDZ),
                sigma_(iSigma),
                caloEta_(iCaloEta),
                caloPhi_(iCaloPhi),
                eta_(iEta),phi_(iPhi),
                quality_(iQuality),
                isPV_(iIsPV),
                alphaF_(alphaF),
                alphaC_(alphaC),
                puppiWeight_(puppiWeight) {}

            float dz() const { return dZ_; }
            float sigma() const { return sigma_; }  
            float caloEta() const { return caloEta_; }
            float caloPhi() const { return caloPhi_; }
	    //iEta,iPhi (usuals)
	    int   iEta() const { return l1tpf::translateIEta(eta_);}
	    int   iPhi() const { return l1tpf::translateIPhi(phi_,eta_);}
	    //iEta,iPhi as they are stored in Arrays
	    int   aEta() const { return l1tpf::translateAEta(iEta());}
	    int   aPhi() const { return l1tpf::translateAPhi(iPhi());}
	    //Center of the trigger tower
	    float dEta() const { return l1tpf::towerEta(l1tpf::translateIEta(eta_));}
	    float dPhi() const { return l1tpf::towerPhi(l1tpf::translateIEta(eta_),l1tpf::translateIPhi(phi_,eta_));}
            //Other stuff
	    float quality() const { return quality_; }
            float alphaF() const { return alphaF_; }
            float alphaC() const { return alphaC_; }
            float puppiWeight() const { return puppiWeight_; }
            int isPV() const { return isPV_; }
            float hOverE() const { return hOverE_; }
            // for L1Tk
            float normalizedChi2() const { return chi2n_; }

	    
            void setCaloEta(float caloEta) { caloEta_ = caloEta; }
            void setCaloPhi(float caloPhi) { caloPhi_ = caloPhi; }
            void setEta(float iEta) { eta_ = iEta; }
            void setPhi(float iPhi) { phi_ = iPhi; }
            void setCaloEtaPhi(float caloEta, float caloPhi) { caloEta_ = caloEta; caloPhi_ = caloPhi; }
            //void setIEtaIPhi(int iEta, int iPhi) { iEta_ = iEta; iPhi_ = iPhi; }

            void setDz(float dz) { dZ_ = dz; }
            void setSigma(float sigma) { sigma_ = sigma; }  
            void setQuality(float quality) { quality_ = quality; }
            void setAlphaF(float alphaF) { alphaF_ = alphaF; }
            void setAlphaC(float alphaC) { alphaC_ = alphaC; }
            void setPuppiWeight(float puppiWeight) { puppiWeight_ = puppiWeight; }
            void setIsPV(int isPV) { isPV_ = isPV; }

            // for HGC 3D clusters, or our own linked ecal+hcal clusters
            void setHOverE(float hOverE) { hOverE_ = hOverE; }

            // for L1Tk
            void setNormalizedChi2(float normalizedChi2) { chi2n_ = normalizedChi2; }

            TLorentzVector tp4() const { 
                TLorentzVector ret;
                ret.SetPtEtaPhiM(pt(),eta(),phi(),mass());
                return ret;
            }
            void setPtEtaPhiM(float pt, float eta, float phi, float mass) { setP4(PolarLorentzVector(pt,eta,phi,mass)); }
            void setPt(float pt) { setP4(PolarLorentzVector(pt,eta(),phi(),mass())); }
    
            void addToP4(const l1tpf::Particle &other) {
                setP4(p4() + other.p4());
            }

        protected:
            float dZ_;
            float sigma_;
            float caloEta_, caloPhi_;
	    float eta_, phi_;
            float quality_;
            int isPV_;
            float hOverE_, chi2n_;
            float alphaF_, alphaC_, puppiWeight_;
    }; // class
} // namespace

#endif
