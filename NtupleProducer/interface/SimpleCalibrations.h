#ifndef FastPUPPI_NtupleProducer_SimpleCalibrations_h
#define FastPUPPI_NtupleProducer_SimpleCalibrations_h
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <vector>

namespace l1tpf {
    class SimpleCorrEm {
        public:
            SimpleCorrEm() {}
            SimpleCorrEm(const edm::ParameterSet &iConfig, const std::string &name) {
                if (iConfig.existsAs<edm::ParameterSet>(name)) {
                    edm::ParameterSet cpset = iConfig.getParameter<edm::ParameterSet>(name);
                    std::vector<double> etaBins = cpset.getParameter<std::vector<double>>("etaBins");
                    std::vector<double> offset = cpset.getParameter<std::vector<double>>("offset");
                    std::vector<double> scale = cpset.getParameter<std::vector<double>>("scale");
                    etas.insert(etas.end(), etaBins.begin(), etaBins.end());
                    scales.insert(scales.end(), scale.begin(), scale.end());
                    offsets.insert(offsets.end(), offset.begin(), offset.end());
                }
            }
            
            bool empty() const { return etas.empty(); }
            float operator()(float pt, float abseta) const {
                for (unsigned int i = 0, n = etas.size(); i < n; ++i) {
                    if (abseta < etas[i]) return (pt - offsets[i])/scales[i];
                }
                return pt;
            }
        protected:
            std::vector<float> etas, scales, offsets;
    };
    class SimpleCorrHad {
        public:
            SimpleCorrHad() {}
            SimpleCorrHad(const edm::ParameterSet &iConfig, const std::string &name) {
                if (iConfig.existsAs<edm::ParameterSet>(name)) {
                    edm::ParameterSet cpset = iConfig.getParameter<edm::ParameterSet>(name);
                    std::vector<double> etaBins = cpset.getParameter<std::vector<double>>("etaBins");
                    std::vector<double> emfBins = cpset.getParameter<std::vector<double>>("emfBins");
                    std::vector<double> offset = cpset.getParameter<std::vector<double>>("offset");
                    std::vector<double> scale = cpset.getParameter<std::vector<double>>("scale");
                    etas.insert(etas.end(), etaBins.begin(), etaBins.end());
                    emfs.insert(emfs.end(), emfBins.begin(), emfBins.end());
                    scales.insert(scales.end(), scale.begin(), scale.end());
                    offsets.insert(offsets.end(), offset.begin(), offset.end());
                }
            }
            
            bool empty() const { return etas.empty(); }
            float operator()(float pt, float abseta, float emf) const {
                unsigned int i = 0, n = emfs.size();
                // step 1: find the first emf bin for this eta
                while (i < n && abseta > etas[i]) i++;
                if (i == n) {
                    //printf("for pt %7.2f eta %4.2f emf %4.2f will not apply any correction (no eta bin found)\n", pt, abseta, emf);
                    return pt;
                }
                unsigned int i2 = i;
                while (i < n && etas[i] == etas[i2] && emf > emfs[i]) i++;
                if (i == n || etas[i] != etas[i2]) {
                    //printf("for pt %7.2f eta %4.2f emf %4.2f will not apply any correction (no emf bin found)\n", pt, abseta, emf);
                    return pt;
                } 
                //printf("for pt %7.2f eta %4.2f emf %4.2f will use bin %d eta [ * , %4.2f ] emf [ * , %4.2f ] offset = %+5.2f scale=%.3f -> corr pt %7.2f\n",
                //        pt, abseta, emf, i, etas[i], emfs[i], offsets[i], scales[i], (pt-offsets[i])/scales[i]);
                return (pt-offsets[i])/scales[i]; 
            } 
        protected:
            std::vector<float> etas, emfs, scales, offsets;
    };
    class SimpleResol { 
        public:
            SimpleResol() {}
            SimpleResol(const edm::ParameterSet &iConfig, const std::string &name) {
                if (iConfig.existsAs<edm::ParameterSet>(name)) {
                    edm::ParameterSet cpset = iConfig.getParameter<edm::ParameterSet>(name);
                    std::vector<double> etaBins = cpset.getParameter<std::vector<double>>("etaBins");
                    std::vector<double> offset = cpset.getParameter<std::vector<double>>("offset");
                    std::vector<double> scale = cpset.getParameter<std::vector<double>>("scale");
                    etas.insert(etas.end(), etaBins.begin(), etaBins.end());
                    scales.insert(scales.end(), scale.begin(), scale.end());
                    offsets.insert(offsets.end(), offset.begin(), offset.end());
                    std::string skind = cpset.getParameter<std::string>("kind");
                    if      (skind == "track") kind = Track;
                    else if (skind == "calo")  kind = Calo;
                    else throw cms::Exception("Configuration", "Bad kind of resolution");
                }
            }
            float operator()(const float pt, const float abseta) const { 
                for (unsigned int i = 0, n = etas.size(); i < n; ++i) {
                    if (abseta < etas[i]) {
                        switch(kind) {
                            case Track: return pt * std::min<float>(1.f, std::hypot(pt*scales[i]*0.001,offsets[i])); 
                            case Calo:  return std::min<float>(pt,pt*scales[i]+offsets[i]); 
                        }
                    }
                }
                return std::min<float>(pt,0.3*pt+7); // saturate to 100% at 10 GeV, and to 30% at high pt
            } 
            bool empty() const { return etas.empty(); }
        protected:
            std::vector<double> etas, offsets, scales; 
            enum Kind { Calo, Track };
            Kind kind;
    };

};

#endif
