#include "../interface/metanalyzer.hh"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"

metanalyzer::metanalyzer(std::string iFile) {
  fNVars = 20;
  for(int i0 = 0; i0 < fNVars; i0++) fVar[i0] = 0;
  fFile = new TFile(iFile.c_str(),"RECREATE");
  fTree = new TTree("met","met");
  fTree->Branch("m_z"   ,&fVar[0],"m_z/D");
  fTree->Branch("pt_z"  ,&fVar[1],"pt_z/D");
  fTree->Branch("eta_z" ,&fVar[2],"eta_z/D");
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
  fTree->Branch("puppimet"   ,&fVar[20],"puppimet/D");
  fTree->Branch("puppimetphi",&fVar[21],"puppimetphi/D");
  fTree->Branch("puppiu1"    ,&fVar[22],"puppiu1/D");
  fTree->Branch("puppiu2"    ,&fVar[23],"puppiu2/D");
  fTree->Branch("tkvtxmet"   ,&fVar[24],"tkvtxmet/D");
  fTree->Branch("tkvtxmetphi",&fVar[25],"tkvtxmetphi/D");
  fTree->Branch("tkvtxu1"    ,&fVar[26],"tkvtxu1/D");
  fTree->Branch("tkvtxu2"    ,&fVar[27],"tkvtxu2/D");
}
void metanalyzer::setZ(std::vector<combiner::Particle> &iParticle) {
  int nmuons = 0;
  for(unsigned   int i0 = 0; i0 < iParticle.size(); i0++) {
    if (iParticle[i0].id == 4) { nmuons++; }
    for(unsigned int i1 = 0; i1 < iParticle.size(); i1++) {
      if(iParticle[i0].Et < 20 || iParticle[i1].Et < 20) continue;
      if(iParticle[i0].id != 4 || iParticle[i1].id != 4) continue;
      TLorentzVector pVec0; pVec0.SetPtEtaPhiM(iParticle[i0].Et,iParticle[i0].Eta,iParticle[i0].Phi,iParticle[i0].M);
      TLorentzVector pVec1; pVec1.SetPtEtaPhiM(iParticle[i1].Et,iParticle[i1].Eta,iParticle[i1].Phi,iParticle[i1].M);
      if((pVec0+pVec1).M() < 60. && (pVec0+pVec1).M() > 120.) continue;
      if(iParticle[i0].charge * iParticle[i1].charge == 1.) continue; //oppo charge
      pVec0+=pVec1;
      fVar[0] = pVec0.M();
      fVar[1] = pVec0.Pt();
      fVar[2] = pVec0.Eta();
      fVar[3] = pVec0.Phi();
      break;
    }
  }
  std::cout << "zinfo: m = " << fVar[0] << ", pt = " << fVar[1] << std::endl;
}
void metanalyzer::setMETRecoil(int iId,std::vector<combiner::Particle> &iParticle,bool iAdd) {
  int lId = 4+iId*4;
  TLorentzVector lVec; lVec.SetPtEtaPhiM(0.,0.,0.,0.);
  TLorentzVector lZ;   lZ  .SetPtEtaPhiM(fVar[1],0.,fVar[3],0.);
  for(unsigned   int i0 = 0; i0 < iParticle.size(); i0++) {
    TLorentzVector pVec1; pVec1.SetPtEtaPhiM(iParticle[i0].Et,0,iParticle[i0].Phi,0);
    lVec-=pVec1;
  }
  if(iAdd) lVec -= lZ;
  fVar[lId+0] = lVec.Pt();
  fVar[lId+1] = lVec.Phi();
  lVec += lZ;
  lVec.RotateZ(-lZ.Phi());
  fVar[lId+2] = lVec.Px();
  fVar[lId+3] = lVec.Py();
}
