#ifndef FASTPUPPI_NTUPLERPRODUCER_L1TPFCandidate_HH
#define FASTPUPPI_NTUPLERPRODUCER_L1TPFCandidate_HH

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include <TLorentzVector.h>

namespace l1tpf {
    class Particle : public reco::LeafCandidate {
        public:
            Particle() {}
            Particle(double iEt,double iEta,double iPhi,double iM,int iId,double iSigma,double iDZ,double iCaloEta=0,double iCaloPhi=0, double iCharge = 0, double iQuality = -999, double iIsPV = 0, float alphaF = -999, float alphaC = -999, float puppiWeight = -99) :
                LeafCandidate(iCharge, reco::LeafCandidate::PolarLorentzVector(iEt,iEta,iPhi,iM), reco::LeafCandidate::Point(), iId),
                dZ_(iDZ),
                sigma_(iSigma),
                caloEta_(iCaloEta),
                caloPhi_(iCaloPhi),
                iEta_(0),iPhi_(0),
                quality_(iQuality),
                isPV_(iIsPV),
                alphaF_(alphaF),
                alphaC_(alphaC),
                puppiWeight_(puppiWeight) {}

            float dz() const { return dZ_; }
            float sigma() const { return sigma_; }  
            float caloEta() const { return caloEta_; }
            float caloPhi() const { return caloPhi_; }
            int iEta() const { return iEta_; }
            int iPhi() const { return iPhi_; }
            float quality() const { return quality_; }
            float alphaF() const { return alphaF_; }
            float alphaC() const { return alphaC_; }
            float puppiWeight() const { return puppiWeight_; }
            int isPV() const { return isPV_; }

            void setCaloEta(float caloEta) { caloEta_ = caloEta; }
            void setCaloPhi(float caloPhi) { caloPhi_ = caloPhi; }
            void setIEta(int iEta) { iEta_ = iEta; }
            void setIPhi(int iPhi) { iPhi_ = iPhi; }
            void setCaloEtaPhi(float caloEta, float caloPhi) { caloEta_ = caloEta; caloPhi_ = caloPhi; }
            void setIEtaIPhi(int iEta, int iPhi) { iEta_ = iEta; iPhi_ = iPhi; }

            void setDz(float dz) { dZ_ = dz; }
            void setSigma(float sigma) { sigma_ = sigma; }  
            void setQuality(float quality) { quality_ = quality; }
            void setAlphaF(float alphaF) { alphaF_ = alphaF; }
            void setAlphaC(float alphaC) { alphaC_ = alphaC; }
            void setPuppiWeight(float puppiWeight) { puppiWeight_ = puppiWeight; }
            void setIsPV(int isPV) { isPV_ = isPV; }

            TLorentzVector tp4() const { 
                TLorentzVector ret;
                ret.SetPtEtaPhiM(pt(),eta(),phi(),mass());
                return ret;
            }
            void setPtEtaPhiM(float pt, float eta, float phi, float mass) { setP4(PolarLorentzVector(pt,eta,phi,mass)); }
            void setPt(float pt) { setP4(PolarLorentzVector(pt,eta(),phi(),mass())); }

        protected:
            float dZ_;
            float sigma_;
            float caloEta_, caloPhi_;
            int   iEta_, iPhi_;
            float quality_;
            int isPV_;
            float alphaF_, alphaC_, puppiWeight_;
    }; // class
} // namespace

#endif
