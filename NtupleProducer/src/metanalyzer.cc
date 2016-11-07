#include "../interface/metanalyzer.hh"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"

metanalyzer::metanalyzer(std::string iFile) {
  fNVars = 30;
  for(int i0 = 0; i0 < fNVars; i0++) fVar[i0] = 0;
  fFile = new TFile(iFile.c_str(),"RECREATE");
  fTree = new TTree("met","met");
  fTree->Branch("m_z"   ,&fVar[0],"m_z/D");
  fTree->Branch("pt_z"  ,&fVar[1],"pt_z/D");
  fTree->Branch("eta_z" ,&fVar[2],"eta_z/D");
  fTree->Branch("phi_z" ,&fVar[3],"phi_z/D");
  fTree->Branch("phi_z" ,&fVar[3],"phi_z/D");
  fTree->Branch("met"   ,&fVar[4],"met/D");
  fTree->Branch("metphi",&fVar[5],"metphi/D");
  fTree->Branch("u1"         ,&fVar[6] ,"u1/D");
  fTree->Branch("u2"         ,&fVar[7] ,"u2/D");
  fTree->Branch("calomet"    ,&fVar[8] ,"calomet/D");
  fTree->Branch("calometphi" ,&fVar[9] ,"calometphi/D");
  fTree->Branch("calou1"     ,&fVar[10],"calou1/D");
  fTree->Branch("calou2"     ,&fVar[11],"calou2/D");
  fTree->Branch("tkmet"      ,&fVar[12],"tkmet/D");
  fTree->Branch("tkmetphi"   ,&fVar[13],"tkmetphi/D");
  fTree->Branch("tku1"       ,&fVar[14],"tku1/D");
  fTree->Branch("tku2"       ,&fVar[15],"tku2/D");
  fTree->Branch("ucalomet"   ,&fVar[16],"ucalomet/D");
  fTree->Branch("ucalometphi",&fVar[17],"ucalometphi/D");
  fTree->Branch("ucalou1"    ,&fVar[18],"ucalou1/D");
  fTree->Branch("ucalou2"    ,&fVar[19],"ucalou2/D");
  fTree->Branch("pvtkmet"    ,&fVar[20],"pvtkmet/D");
  fTree->Branch("pvtkmetphi" ,&fVar[21],"pvtkmetphi/D");
  fTree->Branch("pvtku1"     ,&fVar[22],"pvtku1/D");
  fTree->Branch("pvtku2"     ,&fVar[23],"pvtku2/D");
  fTree->Branch("pupmet"     ,&fVar[24],"pupmet/D");
  fTree->Branch("pupmetphi"  ,&fVar[25],"pupmetphi/D");
  fTree->Branch("pupu1"      ,&fVar[26],"pupu1/D");
  fTree->Branch("pupu2"      ,&fVar[27],"pupu2/D");
  fTree->Branch("dz"         ,&fVar[28],"dz/D");
  fTree->Branch("dzmu"       ,&fVar[29],"dzmu/D");
}
void metanalyzer::setZ(std::vector<combiner::Particle> &iParticle,double iDZ) {
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
      fVar[28] = iDZ;
      fVar[29] = (iParticle[i0].dz()+iParticle[i1].dz())/2.;
      break;
    }
  }
}
void metanalyzer::setMETRecoil(int iId,std::vector<combiner::Particle> &iParticle,bool iAdd) {
  int lId = 4+iId*4;
  TLorentzVector  lVec(0.,0.,0.,0.);
  TLorentzVector lZ(fVar[1],0.,fVar[3],0.);
  for(unsigned   int i0 = 0; i0 < iParticle.size(); i0++) {
    lVec-=iParticle[i0].tp4();
  }
  if(iAdd) lVec -= lZ;
  fVar[lId+0] = lVec.Pt();
  fVar[lId+1] = lVec.Phi();
  lVec += lZ;
  lVec.RotateZ(-lZ.Phi());
  fVar[lId+2] = lVec.Px();
  fVar[lId+3] = lVec.Py();
}
