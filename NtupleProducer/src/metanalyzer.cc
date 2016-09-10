#include "../interface/metanalyzer.hh"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"

metanalyzer::metanalyzer(std::string iFile) {
  fNVars = 20;
  for(int i0 = 0; i0 < fNVars; i0++) fVar[i0] = 0;
  fFile = new TFile(iFile.c_str(),"RECREATE");
  fTree = new TTree("met","met");
  fTree->Branch("m_z"   ,&fVar[0],"fMZ/D");
  fTree->Branch("pt_z"  ,&fVar[1],"fVar[1]/D");
  fTree->Branch("eta_z" ,&fVar[2],"fVar[2]/D");
  fTree->Branch("phi_z" ,&fVar[3],"fVar[3]/D");
  fTree->Branch("met"   ,&fVar[4],"fVar[4]/D");
  fTree->Branch("metphi",&fVar[5],"fVar[5]/D");
  fTree->Branch("u1"         ,&fVar[6] ,"fVar[6]/D");
  fTree->Branch("u2"         ,&fVar[7] ,"fVar[7]/D");
  fTree->Branch("calomet"    ,&fVar[8] ,"fVar[8]/D");
  fTree->Branch("calometphi" ,&fVar[9] ,"fVar[9]/D");
  fTree->Branch("calou1"     ,&fVar[10],"fVar[10]/D");
  fTree->Branch("calou2"     ,&fVar[11],"fVar[11]/D");
  fTree->Branch("tkmet"      ,&fVar[12],"fVar[12]/D");
  fTree->Branch("tkmetphi"   ,&fVar[13],"fVar[13]/D");
  fTree->Branch("tku1"       ,&fVar[14],"fVar[14]/D");
  fTree->Branch("tku2"       ,&fVar[15],"fVar[15]/D");
  fTree->Branch("ucalomet"   ,&fVar[16],"fVar[16]/D");
  fTree->Branch("ucalometphi",&fVar[17],"fVar[17]/D");
  fTree->Branch("ucalou1"    ,&fVar[18],"fVar[18]/D");
  fTree->Branch("ucalou2"    ,&fVar[19],"fVar[19]/D");
}
void metanalyzer::setZ(std::vector<combiner::Particle> &iParticle) {
  for(unsigned   int i0 = 0; i0 < iParticle.size(); i0++) {
    for(unsigned int i1 = 0; i1 < iParticle.size(); i1++) {
      if(iParticle[i0].Et < 20 || iParticle[i1].Et < 20) continue;
      TLorentzVector pVec0; pVec0.SetPtEtaPhiM(iParticle[i0].Et,iParticle[i0].Eta,iParticle[i0].Phi,iParticle[i0].M);
      TLorentzVector pVec1; pVec1.SetPtEtaPhiM(iParticle[i1].Et,iParticle[i1].Eta,iParticle[i1].Phi,iParticle[i1].M);
      if((pVec0+pVec1).M() < 60.) continue;
      pVec0+=pVec1;
      fVar[0] = pVec0.M();
      fVar[1] = pVec0.Pt();
      fVar[2] = pVec0.Eta();
      fVar[3] = pVec0.Phi();
      break;
    }
  }
}
void metanalyzer::setMETRecoil(int iId,std::vector<combiner::Particle> &iParticle,bool iAdd) {
  int lId = 4+iId*4;
  TLorentzVector lVec; lVec.SetPtEtaPhiM(0.,0.,0.,0.);
  TLorentzVector lZ;   lZ  .SetPtEtaPhiM(fVar[1],0.,fVar[3],0.);
  for(unsigned   int i0 = 0; i0 < iParticle.size(); i0++) {
    TLorentzVector pVec1; pVec1.SetPtEtaPhiM(iParticle[i0].Et,0,iParticle[i0].Phi,0);
    lVec-=pVec1;
  }
  //if(iAdd) lVec -= lZ;
  fVar[lId+0] = lVec.Pt();
  fVar[lId+1] = lVec.Phi();
  lVec += lZ;
  lVec.RotateZ(-lZ.Phi());
  fVar[lId+2] = lVec.Px();
  fVar[lId+3] = lVec.Py();
}
