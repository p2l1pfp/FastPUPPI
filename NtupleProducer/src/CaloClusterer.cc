#include "../interface/CaloClusterer.h"
#include <DataFormats/Math/interface/deltaPhi.h>
#include "FWCore/ParameterSet/interface/ParameterSet.h"

const float l1pf_calo::Stage1Grid::towerEtas_[l1pf_calo::Stage1Grid::nEta_] = {0,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696,0.783,0.870,0.957,1.044,1.131,1.218,1.305,1.392,1.479,1.566,1.653,1.740,1.830,1.930,2.043,2.172,2.322,2.5,2.650,2.853,3.139,3.314,3.489,3.664,3.839,4.013,4.191,4.363,4.538,4.716,4.889,5.191};

l1pf_calo::Stage1Grid::Stage1Grid() :
            Grid(2*((ietaCoarse_-1)*nPhi_ + (ietaVeryCoarse_-ietaCoarse_)*(nPhi_/2) + (nEta_-ietaVeryCoarse_+1) * (nPhi_/4))),
            cell_map_(2*nEta_*nPhi_, -1)
{
    int icell = 0;
    for (int ie = -nEta_; ie <= nEta_; ++ie) {
        int absie = std::abs(ie);
        for (int iph = 1; iph <= nPhi_; ++iph) {
            if (!valid_ieta_iphi(ie,iph)) continue;
            ieta_[icell] = ie;
            iphi_[icell] = iph;
            eta_[icell]      = (ie > 0 ? 0.5 : -0.5)*(towerEtas_[absie-1] + towerEtas_[absie]);
            etaWidth_[icell] = (towerEtas_[absie] - towerEtas_[absie-1]);
            phiWidth_[icell] = 2*M_PI/nPhi_;
            if (absie >= ietaVeryCoarse_)  phiWidth_[icell] *= 4;
            else if (absie >= ietaCoarse_) phiWidth_[icell] *= 2;
            phi_[icell]      = (iph-1)*2*M_PI/nPhi_ + 0.5*phiWidth_[icell];
            if (phi_[icell] > M_PI) phi_[icell] -= 2*M_PI;
            std::fill(neighbours_[icell].begin(), neighbours_[icell].end(), -1);
            cell_map_[(ie+nEta_) + 2 * nEta_*(iph-1)] = icell;
            icell++;
        }
    }
    assert(unsigned(icell) == ncells_);
    // now link the cells
    for (icell = 0; icell < int(ncells_); ++icell) {
        int ie = ieta_[icell], iph = iphi_[icell];
        int ineigh = 0;
        for (int deta = -1; deta <= +1; ++deta) {
            for (int dphi = -1; dphi <= +1; ++dphi) {
                if (deta == 0 && dphi == 0) continue;
                neighbours_[icell][ineigh++] = imove(ie, iph, deta, dphi);
            }
        }
    } 
    //// consistency check 1: check that find_cell works
    //for (float teta = 0; teta <= 5.0; teta += 0.02) {
    //    for (float tphi = -M_PI; tphi <= M_PI; tphi += 0.02) {
    //        find_cell(+teta, tphi);
    //        find_cell(-teta, tphi);
    //    }
    //}
}

int l1pf_calo::Stage1Grid::find_cell(float eta, float phi) const {
    int ieta = (eta != 0) ? std::distance(towerEtas_, std::lower_bound(towerEtas_, towerEtas_+nEta_, std::abs(eta))) : 1;
    assert(ieta > 0 && ieta < nEta_);
    if (ieta > nEta_) ieta = nEta_;
    if (eta < 0) ieta = -ieta;
    if (phi > 2*M_PI) phi -= 2*M_PI;
    if (phi < 0) phi += 2*M_PI;
    int iphi = std::floor(phi * nPhi_/(2*M_PI));
    if (phi >= 2*M_PI) iphi = nPhi_-1; // fix corner case due to roundings etc
    assert(iphi < nPhi_);
    if      (std::abs(ieta) >= ietaVeryCoarse_) iphi -= (iphi%4);
    else if (std::abs(ieta) >= ietaCoarse_)     iphi -= (iphi%2);
    iphi += 1;
    //if (!valid_ieta_iphi(ieta,iphi)) {
    //    printf("Error in finding cell for eta %+7.4f phi %+7.4f, got ieta = %+3d iphi %2d which is not valid\n",
    //        eta, phi, ieta, iphi);
    //}
    assert(valid_ieta_iphi(ieta,iphi));
    int icell = ifind_cell(ieta,iphi);
    assert(icell != -1);

    //if (std::abs(eta - eta_[icell]) > 0.501*etaWidth_[icell] || std::abs(deltaPhi(phi, phi_[icell])) > 0.501*phiWidth_[icell]) {
    //    printf("Mismatch in finding cell for eta %+7.4f phi %+7.4f, got ieta = %+3d iphi %2d which has eta %+7.4f +- %.4f phi %+7.4f +- %.4f ; deta = %+7.4f dphi = %+7.4f\n",
    //        eta, phi, ieta, iphi, eta_[icell], etaWidth_[icell], phi_[icell], phiWidth_[icell], eta - eta_[icell], deltaPhi(phi, phi_[icell]));
    //}
    //assert(std::abs(eta - eta_[icell]) <= 0.5*etaWidth_[icell]);  
    //assert(std::abs(deltaPhi(phi, phi_[icell])) <= 0.5*phiWidth_[icell]);  
    return icell;
}

int l1pf_calo::Stage1Grid::imove(int ieta, int iphi, int deta, int dphi) {
    int ie = ieta, iph = iphi;
    switch (deta) {
        case -1: ie = (ie == -nEta_ ? 0 : (ie == +1 ? -1 : ie-1)); break;
        case +1: ie = (ie == +nEta_ ? 0 : (ie == -1 ? +1 : ie+1)); break;
        case 0: break;
        default: assert(false);
    };
    if (ie == 0) return -1;
    switch (dphi) {
        case -1: iph = (iph ==   1   ? nPhi_ : iph-1); break;
        case +1: iph = (iph == nPhi_ ?   1   : iph+1); break;
        case 0: break;
        default: assert(false);
    };
    if (!valid_ieta_iphi(ie,iph)) return -1;
    int icell = ifind_cell(ie,iph);
    assert(!(ie == ieta && iph == iphi));
    assert(icell != -1);
    assert(icell != ifind_cell(ieta,iphi));
    return icell;
}

l1pf_calo::SingleCaloClusterer::SingleCaloClusterer(const edm::ParameterSet &pset) :
    grid_(new Stage1Grid()),
    rawet_(*grid_),
    precluster_(*grid_),
    cluster_(*grid_),
    zsEt_(pset.getParameter<double>("zsEt")),
    seedEt_(pset.getParameter<double>("seedEt")),
    minClusterEt_(pset.getParameter<double>("minClusterEt"))
{
}

l1pf_calo::SingleCaloClusterer::~SingleCaloClusterer() 
{
}

void l1pf_calo::SingleCaloClusterer::run() 
{
    unsigned int i, ncells = grid_->size();

    // kill zeros
    for (i = 0; i < ncells; ++i) {
        if (rawet_[i] < zsEt_) rawet_[i] = 0;
    }

    precluster_.clear();
    // pre-cluster step 1: at each cell, set the value equal to itself if it's a local maxima, zero otherwise
    // can be done in parallel on all cells
    for (i = 0; i < ncells; ++i) {
        if (rawet_[i] > seedEt_) {
            precluster_[i].ptLocalMax = rawet_[i];
            //printf("   candidate precluster pt %7.2f at %4d (ieta %+3d iphi %2d)\n",  rawet_[i], i, grid_->ieta(i), grid_->iphi(i));
            for (int ineigh = 0; ineigh <= 3; ++ineigh) {
                if (rawet_.neigh(i,ineigh) >  rawet_[i]) precluster_[i].ptLocalMax = 0;
                //int ncell = grid_->neighbour(i,ineigh);
                //if (ncell == -1) printf("   \t neigh %d is null\n", ineigh);
                //else printf("   \t neigh %d at %4d (ieta %+3d iphi %2d) has pt %7.2f: comparison %1d \n", ineigh, ncell, grid_->ieta(ncell), grid_->iphi(ncell), rawet_[ncell], precluster_[i].ptLocalMax > 0);
            }
            for (int ineigh = 4; ineigh <  8; ++ineigh) {
                if (rawet_.neigh(i,ineigh) >= rawet_[i]) precluster_[i].ptLocalMax = 0;
                //int ncell = grid_->neighbour(i,ineigh);
                //if (ncell == -1) printf("   \t neigh %d is null\n", ineigh);
                //else printf("   \t neigh %d at %4d (ieta %+3d iphi %2d) has pt %7.2f: comparison %1d \n", ineigh, ncell, grid_->ieta(ncell), grid_->iphi(ncell), rawet_[ncell], precluster_[i].ptLocalMax > 0);
            }
        }
    }
    // pre-cluster step 2: at each non-max cell, add up ptLocalMax of neighbours
    //                     could be done for all cells, but it's not needed for max cells since they can't have ptLocalMax around them
    for (i = 0; i < ncells; ++i) {
        if (precluster_[i].ptLocalMax == 0) {
            float tot = 0;
            for (int ineigh = 0; ineigh < 8; ++ineigh) {
                tot += precluster_.neigh(i,ineigh).ptLocalMax;
            }
            precluster_[i].ptOverNeighLocalMaxSum = tot ? rawet_[i]/tot : 0;
        }
    }

    cluster_.clear();
    // cluster: at each localMax cell, take itself plus the weighted contributions of the neighbours
    for (i = 0; i < ncells; ++i) {
        if (precluster_[i].ptLocalMax > 0) {
            float myet = rawet_[i];
            float tot  = myet;
            float avg_eta = 0;
            float avg_phi = 0;
            for (int ineigh = 0; ineigh < 8; ++ineigh) {
                int ineighcell = grid_->neighbour(i, ineigh);
                if (ineighcell == -1) continue; // skip dummy cells
                float fracet = myet * precluster_.neigh(i,ineigh).ptOverNeighLocalMaxSum;
                tot  += fracet;
                avg_eta += fracet * (grid_->eta(ineighcell) - grid_->eta(i));
                avg_phi += fracet * deltaPhi(grid_->phi(ineighcell), grid_->phi(i));
            }
            if (tot > minClusterEt_) {
                cluster_[i].et  = tot;
                cluster_[i].eta = grid_->eta(i) + avg_eta / tot;
                cluster_[i].phi = grid_->phi(i) + avg_phi / tot;
                // wrap around phi
                if (cluster_[i].phi >  M_PI) cluster_[i].phi -= 2*M_PI;
                if (cluster_[i].phi < -M_PI) cluster_[i].phi += 2*M_PI;
            }
        }
    }


}

l1pf_calo::SimpleCaloLinker::SimpleCaloLinker(const edm::ParameterSet &pset, const SingleCaloClusterer & ecal,  const SingleCaloClusterer & hcal) :
    grid_(new Stage1Grid()),
    ecal_(ecal), hcal_(hcal),
    ecalToHCal_(*grid_),
    cluster_(*grid_),
    hoeCut_(pset.getParameter<double>("hoeCut")),
    minPhotonEt_(pset.getParameter<double>("minPhotonEt")),
    minHadronEt_(pset.getParameter<double>("minHadronEt")),
    useCorrectedEcal_(pset.getParameter<bool>("useCorrectedEcal"))
{
    assert(grid_->size() == ecal.raw().grid().size());
    assert(grid_->size() == hcal.raw().grid().size());
}

l1pf_calo::SimpleCaloLinker::~SimpleCaloLinker() 
{
}

void l1pf_calo::SimpleCaloLinker::run() 
{
    unsigned int i, ncells = grid_->size();

    const EtGrid & hraw = hcal_.raw();
    const ClusterGrid & ecals = ecal_.clusters();
    const ClusterGrid & hcals = hcal_.clusters();

    // for each ECal cluster, get the corresponding HCal cluster and the sum of the neighbour HCal clusters
    ecalToHCal_.clear();
    for (i = 0; i < ncells; ++i) {
        if (ecals[i].et > 0) {
            if (hcals[i].et > 0) {
                ecalToHCal_[i].ptLocalMax = hcals[i].et;
            } else {
                float tot = 0;
                for (int ineigh = 0; ineigh < 8; ++ineigh) {
                    tot += hcals.neigh(i,ineigh).et;
                }
                ecalToHCal_[i].ptOverNeighLocalMaxSum = tot ? (useCorrectedEcal_ ? ecals[i].et_corr : ecals[i].et)/tot : 0;
            }
        }
    }

    cluster_.clear();
    // promote HCal clusters to final clusters
    for (i = 0; i < ncells; ++i) {
        if (hcals[i].et > 0) {
            if (ecalToHCal_[i].ptLocalMax > 0) {
                // direct linking is easy
                cluster_[i].ecal_et = (useCorrectedEcal_ ? ecals[i].et_corr : ecals[i].et);
                cluster_[i].hcal_et = hcals[i].et;
                cluster_[i].et = cluster_[i].ecal_et + cluster_[i].hcal_et;
                float wecal = cluster_[i].ecal_et/cluster_[i].et, whcal = 1.0 - wecal;
                cluster_[i].eta = ecals[i].eta * wecal + hcals[i].eta * whcal;
                cluster_[i].phi = ecals[i].phi * wecal + hcals[i].phi * whcal;
                // wrap around phi
                if (cluster_[i].phi >  M_PI) cluster_[i].phi -= 2*M_PI;
                if (cluster_[i].phi < -M_PI) cluster_[i].phi += 2*M_PI;
            } else {
                // sidewas linking is more annonying
                float myet = hcals[i].et;
                float etot = 0;
                float avg_eta = 0;
                float avg_phi = 0;
                for (int ineigh = 0; ineigh < 8; ++ineigh) {
                    int ineighcell = grid_->neighbour(i, ineigh);
                    if (ineighcell == -1) continue; // skip dummy cells
                    float fracet = myet * ecalToHCal_.neigh(i,ineigh).ptOverNeighLocalMaxSum;
                    etot  += fracet;
                    avg_eta += fracet * (grid_->eta(ineighcell) - grid_->eta(i));
                    avg_phi += fracet * deltaPhi(grid_->phi(ineighcell), grid_->phi(i));
                }
                cluster_[i].hcal_et = hcals[i].et;
                cluster_[i].ecal_et = etot;
                cluster_[i].et  = myet + etot;
                cluster_[i].eta = hcals[i].eta + avg_eta / cluster_[i].et;
                cluster_[i].phi = hcals[i].phi + avg_phi / cluster_[i].et;
                // wrap around phi
                if (cluster_[i].phi >  M_PI) cluster_[i].phi -= 2*M_PI;
                if (cluster_[i].phi < -M_PI) cluster_[i].phi += 2*M_PI;
            }
        }
    }

    // promote Unlinked ECal clusters to final clusters
    for (i = 0; i < ncells; ++i) {
        if (ecals[i].et > 0 && ecalToHCal_[i].ptLocalMax == 0 && ecalToHCal_[i].ptOverNeighLocalMaxSum == 0) {
            // direct linking is easy
            cluster_[i].ecal_et = (useCorrectedEcal_ ? ecals[i].et_corr : ecals[i].et);
            cluster_[i].hcal_et = hraw[i];
            cluster_[i].et = cluster_[i].ecal_et + cluster_[i].hcal_et;
            cluster_[i].eta = ecals[i].eta;
            cluster_[i].phi = ecals[i].phi;
            // no need to wrap around phi
        }
    }
 
}



std::vector<l1tpf::Particle> l1pf_calo::SimpleCaloLinker::fetch(bool corrected) const {
    std::vector<l1tpf::Particle> ret;
    for (unsigned int i = 0, ncells = grid_->size(); i < ncells; ++i) {
        if (cluster_[i].et > 0) {
            bool photon = (cluster_[i].hcal_et < hoeCut_* cluster_[i].ecal_et);
            if (cluster_[i].et > (photon ? minPhotonEt_ : minHadronEt_)) {
                ret.emplace_back(corrected ? cluster_[i].et_corr : cluster_[i].et, cluster_[i].eta, cluster_[i].phi, 0.0, photon ? 3 : 2);  
                ret.back().setCaloEtaPhi(cluster_[i].eta, cluster_[i].phi);
            }
        }
    }
    return ret;
}
