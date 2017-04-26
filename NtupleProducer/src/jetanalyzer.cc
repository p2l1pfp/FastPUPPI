#include "../interface/jetanalyzer.hh"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/JetDefinition.hh>

jetanalyzer::jetanalyzer(std::string iFile, int debug) {
  fFile = new TFile(iFile.c_str(),"RECREATE");
  fTree = new TTree("jet","jet");
  //Reco
  std::vector<std::string> lVec;
  lVec.push_back("jpt_1");
  lVec.push_back("jeta_1");
  lVec.push_back("jphi_1");
  lVec.push_back("jm_1" );
  lVec.push_back("jdz_1");
  lVec.push_back("jpt_2");
  lVec.push_back("jeta_2");
  lVec.push_back("jphi_2");
  lVec.push_back("jm_2" );
  lVec.push_back("jdz_2");
  lVec.push_back("njets");
  lVec.push_back("jptgen_1");
  lVec.push_back("jdrgen_1");
  lVec.push_back("jptgen_2");
  lVec.push_back("jdfgen_2");
  int lNIndex = 7;
  fSize  = lVec.size();
  fBase  = 6;
  const int lNVars = fBase+lNIndex*fSize;
  int lIndex = fBase;
  fVar = new double[lNVars];
  for(int i0 = 0; i0 < lNVars; i0++) fVar[i0] = 0;
  fTree->Branch("m_z"   ,&fVar[0] ,"m_z/D");
  fTree->Branch("pt_z"  ,&fVar[1] ,"pt_z/D");
  fTree->Branch("eta_z" ,&fVar[2] ,"eta_z/D");
  fTree->Branch("phi_z" ,&fVar[3] ,"phi_z/D");
  fTree->Branch("dzmu"  ,&fVar[4] ,"dzmu/D");
  fTree->Branch("dz"    ,&fVar[5] ,"dzmu/D");
  for(unsigned int i0 = 0; i0 < lVec.size(); i0++) fTree->Branch(     lVec[i0] .c_str(),&fVar[i0+lIndex],(    lVec[i0]+"/D").c_str());
  lIndex += lVec.size();
  for(unsigned int i0 = 0; i0 < lVec.size(); i0++) fTree->Branch(("g"+lVec[i0]).c_str(),&fVar[i0+lIndex],("g"+lVec[i0]+"/D").c_str());
  lIndex += lVec.size();
  for(unsigned int i0 = 0; i0 < lVec.size(); i0++) fTree->Branch(("tk"+lVec[i0]).c_str(),&fVar[i0+lIndex],("tk"+lVec[i0]+"/D").c_str());
  lIndex += lVec.size();
  for(unsigned int i0 = 0; i0 < lVec.size(); i0++) fTree->Branch(("pvtk"+lVec[i0]).c_str(),&fVar[i0+lIndex],("pvtk"+lVec[i0]+"/D").c_str());
  lIndex += lVec.size();
  for(unsigned int i0 = 0; i0 < lVec.size(); i0++) fTree->Branch(("ucalo"+lVec[i0]).c_str(),&fVar[i0+lIndex],("ucalo"+lVec[i0]+"/D").c_str());
  lIndex += lVec.size();
  for(unsigned int i0 = 0; i0 < lVec.size(); i0++) fTree->Branch(("calo"+lVec[i0]).c_str(),&fVar[i0+lIndex],("calo"+lVec[i0]+"/D").c_str());
  lIndex += lVec.size();
  for(unsigned int i0 = 0; i0 < lVec.size(); i0++) fTree->Branch(("pup"+lVec[i0]).c_str(),&fVar[i0+lIndex],("pup"+lVec[i0]+"/D").c_str());
  // debug
  fDebug = debug;
}
void jetanalyzer::setZ(std::vector<combiner::Particle> &iParticle,double iDZ) {
  int nmuons = 0;
  for(unsigned   int i0 = 0; i0 < iParticle.size(); i0++) {
    if (iParticle[i0].pdgId() == 4) { nmuons++; }
    for(unsigned int i1 = 0; i1 < iParticle.size(); i1++) {
      if(iParticle[i0].pt() < 20 || iParticle[i1].pt() < 20) continue;
      if(iParticle[i0].pdgId() != 4 || iParticle[i1].pdgId() != 4) continue;
      auto pVec0 = iParticle[i0].p4(), pVec1 = iParticle[i1].p4();
      if((pVec0+pVec1).M() < 60. && (pVec0+pVec1).M() > 120.) continue;
      if(iParticle[i0].charge() * iParticle[i1].charge() == 1.) continue; //oppo charge
      pVec0+=pVec1;
      fVar[0]  = pVec0.M();
      fVar[1]  = pVec0.Pt();
      fVar[2]  = pVec0.Eta();
      fVar[3]  = pVec0.Phi();
      fVar[4] = (iParticle[i0].dz()+iParticle[i1].dz())/2.;
      fVar[5] = iDZ;
      break;
    }
  }
}
void jetanalyzer::setJets(std::vector<combiner::Particle> &iParticles,int iIndex) {
  int lIndex = fBase+iIndex*fSize;
  if (fDebug) std::cout << "---> " << iIndex << " -- " << fSize << " -- " << lIndex << std::endl;
  std::vector < fastjet::PseudoJet > particles;
  for(unsigned   int i0 = 0; i0 < iParticles.size(); i0++) {
    if(iParticles[i0].pt() < 0.1 || iParticles[i0].pdgId() == 4) continue;
    fastjet::PseudoJet curpar   = fastjet::PseudoJet(iParticles[i0].px(),iParticles[i0].py(),iParticles[i0].pz(),iParticles[i0].energy()); 
    curpar.set_user_index( i0 );
    particles.push_back(curpar);
  }
  std::vector < fastjet::PseudoJet > lJets = cluster(particles,0.4,30.);
  fVar[lIndex+10] = lJets.size();
  if(lJets.size() > 0) { 
    fVar[lIndex]   = lJets[0].pt();
    fVar[lIndex+1] = lJets[0].eta();
    fVar[lIndex+2] = lJets[0].phi();
    fVar[lIndex+3] = lJets[0].m();
    fVar[lIndex+4] = dz(lJets[0],iParticles);
    fVar[lIndex+11] = genmatch(1,lJets[0],fGenJets);
    fVar[lIndex+12] = genmatch(0,lJets[0],fGenJets);
  }
  if(lJets.size() > 1) { 
    fVar[lIndex+5] = lJets[1].pt();
    fVar[lIndex+6] = lJets[1].eta();
    fVar[lIndex+7] = lJets[1].phi();
    fVar[lIndex+8] = lJets[1].m();
    fVar[lIndex+9] = dz(lJets[1],iParticles);
    fVar[lIndex+13] = genmatch(1,lJets[1],fGenJets);
    fVar[lIndex+14] = genmatch(0,lJets[1],fGenJets);
  }
}
void jetanalyzer::setGenJets(const reco::GenParticleCollection &iGenParticles,int iIndex) { 
  int lIndex = fBase+iIndex*fSize;
  std::vector < fastjet::PseudoJet > particles;
  for (reco::GenParticleCollection::const_iterator itGenP = iGenParticles.begin(); itGenP!=iGenParticles.end(); ++itGenP) {
    if(itGenP->status() != 1 && abs(itGenP->pdgId()) != 13) continue;
    fastjet::PseudoJet curpar   = fastjet::PseudoJet(itGenP->px(),itGenP->py(),itGenP->pz(),itGenP->energy()); 
    particles.push_back(curpar);
  }
  if(particles.size() < 4) return;
  std::vector < fastjet::PseudoJet > lJets = cluster(particles,0.4,30.);
  fVar[lIndex+10] = lJets.size();
  if(lJets.size() > 0) { 
    fVar[lIndex]   = lJets[0].pt();
    fVar[lIndex+1] = lJets[0].eta();
    fVar[lIndex+2] = lJets[0].phi();
    fVar[lIndex+3] = lJets[0].m();
    fVar[lIndex+4] = -999;
  }
  if(lJets.size() > 1) { 
    fVar[lIndex+5] = lJets[1].pt();
    fVar[lIndex+6] = lJets[1].eta();
    fVar[lIndex+7] = lJets[1].phi();
    fVar[lIndex+8] = lJets[1].m();
    fVar[lIndex+9] = -999;
  }
  fGenJets.clear();
  for(unsigned int i0 = 0; i0 < lJets.size(); i0++) fGenJets.emplace_back(lJets[i0]);
}
std::vector<fastjet::PseudoJet> jetanalyzer::cluster(std::vector < fastjet::PseudoJet > &particles, double iRadius,double iPt){
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, iRadius); 
  int activeAreaRepeats = 1;
  double ghostArea = 0.01;
  double ghostEtaMax = 7.0;
  fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
  fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
  fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(particles, jetDef, fjAreaDefinition);
  std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering->inclusive_jets(iPt));
  //thisClustering->delete_self_when_unused();
  return out_jets;
}
double jetanalyzer::dz(fastjet::PseudoJet &iJet,std::vector<combiner::Particle> &iParticles) {
  double lDZ = 0;
  double lPtTot = 0;
  for(unsigned int i0 = 0; i0 < iJet.constituents().size(); i0++) { 
    combiner::Particle pParticle = iParticles[iJet.constituents()[i0].user_index()];
    lDZ = pParticle.dz() * pParticle.pt();
    lPtTot += pParticle.pt();
  }
  if(lPtTot > 0) lDZ/=lPtTot;
  return lDZ;
}
double jetanalyzer::genmatch(int iId, fastjet::PseudoJet &matchjet,std::vector < fastjet::PseudoJet > &genjets){
  double lDR = 1000;
  double lPt = -1;
  for(unsigned int i0 = 0; i0 < genjets.size(); i0++) { 
    double pDR = matchjet.delta_R(genjets[i0]);
    if(pDR > 0.25) continue;
    if(lDR < pDR)  continue;
    lDR = pDR;
    lPt = genjets[i0].pt();
  }
  if(iId == 0) return lDR;
  return lPt;
}
