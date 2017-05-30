#ifndef FASTPUPPI_NTUPLERPRODUCER_ALTERNATIVEPF_H
#define FASTPUPPI_NTUPLERPRODUCER_ALTERNATIVEPF_H

// more PF algos
#include "FastPUPPI/NtupleProducer/interface/DiscretePF.h"

namespace l1tpf_int { 
class PFAlgo3 : public PFAlgo {
    public:
        PFAlgo3( const edm::ParameterSet& ) ;
        virtual void runPF(Region &r) const override;
    protected:
        float drMatchEm_, ptMinFracMatchEm_, drMatchEmHad_;
        bool caloReLinkStep_; float caloReLinkDr_, caloReLinkThreshold_;
        bool sumTkCaloErr2_, ecalPriority_;
        unsigned int tightTrackMinStubs_; float tightTrackMaxChi2_, tightTrackMaxInvisiblePt_;
        enum GoodTrackStatus { GoodTK_Calo_TkPt=0, GoodTK_Calo_TkCaloPt=1, GoodTk_Calo_CaloPt=2, GoodTK_NoCalo=3 };
        enum BadTrackStatus { BadTK_NoCalo=1 };
  };

} // end namespace

#endif
