#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include <string>
#include <sstream>
#include <iostream>


TH1F *plot(std::string iName,std::string iVar,int iId,std::string iCut,int iColor,int iNBins,float iMin,float iMax,TTree *iTree) { 
  std::stringstream lSS; lSS << iName << iId;
  TH1F *lH = new TH1F(lSS.str().c_str(),lSS.str().c_str(),iNBins,iMin,iMax);//150.0);
  iTree->Draw((iVar+">>"+lSS.str()).c_str(),("pt > 1 " + iCut).c_str());
  //lH[iId]->SetFillStyle(3001);
  //lH[iId]->SetFillColor(iColor);
  lH->SetLineColor(iColor);
  lH->SetLineWidth(2);
  lH->SetMarkerSize(0);
  lH->SetMarkerColor(iColor);
  //lH->GetXaxis()->SetTitle("Iso");
  lH->GetYaxis()->SetTitle("N");
  return lH;
}

TGraphAsymmErrors *plotEff(std::string iName,int iId,std::string iBCut,std::string iCut,int iColor,TTree *iTree,float iScale=1.0) { 
  std::stringstream lSS; lSS << iName << iId;   std::stringstream lS0; lS0 << iName << iId << "F";
  double lX[] = {0,10,20,30,40,50,70,90,150,200};
  TH1F *lH0 = new TH1F(lSS.str().c_str(),lSS.str().c_str(),9,lX);//10,0.0,200.0); 
  TH1F *lH1 = new TH1F(lS0.str().c_str(),lS0.str().c_str(),9,lX);//10,0.0,200.0);
  //iTree->Draw(("chargedIso>>"+lSS.str()).c_str(),("pt > 30" + iCut).c_str());
  std::stringstream lScale; lScale << "*" << iScale;
  iTree->Draw(("genpt>>"+lSS.str()).c_str(),("genpt > 1. && abs(geneta) < 2.5"+iBCut).c_str());
  iTree->Draw(("genpt>>"+lS0.str()).c_str(),("gendr > 0. && abs(geneta) < 2.5 && genpt > 1 && pt > 10"+lScale.str()+iBCut + iCut).c_str());
  //lH[iId]->SetFillStyle(3001);
  //lH[iId]->SetFillColor(iColor);
  TGraphAsymmErrors* lHD = new TGraphAsymmErrors(lH1,lH0);
  //lHD->Sumw2();
  lHD->SetLineColor(iColor);
  lHD->GetYaxis()->SetTitle("Eff");
  lHD->GetXaxis()->SetTitle("gen p_{T} (GeV)");
  lHD->GetYaxis()->SetRangeUser(0.,1.0);
  lHD->SetMarkerStyle(kFullCircle); 
  lHD->SetMarkerColor(iColor);
  return lHD;
}

TGraphAsymmErrors *plotEffEta(std::string iName,int iId,std::string iBCut,std::string iCut,int iColor,TTree *iTree,float iScale=1.0) { 
  std::stringstream lSS; lSS << iName << iId;   std::stringstream lS0; lS0 << iName << iId << "F";
  double lX[] = {-2.5,-2.3,-2.0,-1.8,-1.6,-1.4,-1.0,-0.6,-0.2,0.2,0.6,1.0,1.4,1.6,1.8,2.0,2.3,2.5};
  TH1F *lH0 = new TH1F(lSS.str().c_str(),lSS.str().c_str(),17,lX);//10,0.0,200.0); 
  TH1F *lH1 = new TH1F(lS0.str().c_str(),lS0.str().c_str(),17,lX);//10,0.0,200.0);
  //iTree->Draw(("chargedIso>>"+lSS.str()).c_str(),("pt > 30" + iCut).c_str());
  std::stringstream lScale; lScale << "*" << iScale;
  iTree->Draw(("geneta>>"+lSS.str()).c_str(),("genpt > 40. && abs(geneta) < 2.5"+iBCut).c_str());
  iTree->Draw(("geneta>>"+lS0.str()).c_str(),("gendr > 0. && abs(geneta) < 2.5 && genpt > 40 && pt > 10"+lScale.str()+iBCut + iCut).c_str());
  //lH[iId]->SetFillStyle(3001);
  //lH[iId]->SetFillColor(iColor);
  TGraphAsymmErrors* lHD = new TGraphAsymmErrors(lH1,lH0);
  //lHD->Sumw2();
  lHD->SetLineColor(iColor);
  lHD->GetYaxis()->SetTitle("Eff");
  lHD->GetXaxis()->SetTitle("gen #eta");
  lHD->GetYaxis()->SetRangeUser(0.,1.0);
  lHD->SetMarkerStyle(kFullCircle); 
  lHD->SetMarkerColor(iColor);
  return lHD;
}

TH1F *plotRate(std::string iName,int iId,std::string iBCut,std::string iCut,int iColor,TTree *iTree,float iScale=1.0) { 
  std::stringstream lSS; lSS << iName << iId;   std::stringstream lS0; lS0 << iName << iId << "F";
  std::stringstream lScale; lScale << "*" << iScale;
  double lX[] = {0,10,20,30,40,50,70,90,150,200};
  TH1F *lH0 = new TH1F(lSS.str().c_str(),lSS.str().c_str(),9,lX);
  TH1F *lH1 = new TH1F(lS0.str().c_str(),lS0.str().c_str(),9,lX);
  //iTree->Draw(("chargedIso>>"+lSS.str()).c_str(),("pt > 30" + iCut).c_str());
  iTree->Draw(("pt"+lScale.str()+">>"+lSS.str()).c_str(),("pt"+lScale.str()+" > 10"+iBCut).c_str());
  iTree->Draw(("pt"+lScale.str()+">>"+lS0.str()).c_str(),("pt"+lScale.str()+" > 10"+iBCut + iCut).c_str());
  //lH[iId]->SetFillStyle(3001);
  //lH[iId]->SetFillColor(iColor);
  lH1->Sumw2();
  //lH1->Divide(lH0);
  lH1->SetLineColor(iColor);
  lH1->GetYaxis()->SetTitle("Rate");
  lH1->GetXaxis()->SetTitle("p_{T} scaled (GeV)");
  lH1->SetMarkerStyle(kFullCircle); 
  lH1->SetMarkerColor(iColor);
  lH1->Scale(400.);
  return lH1;
}
float unique(TTree* iTree) { 
  ULong_t lLast = -1; 
  int lTotal = 0; 
  int lCount = 0; 
  ULong64_t lId = 0; iTree->SetBranchAddress("event2",&lId); 
  float lPt = 0;     iTree->SetBranchAddress("pt",&lPt); 
  for(int i0 = 0; i0 < iTree->GetEntriesFast(); i0++) { 
    iTree->GetEntry(i0); 
    if(lPt < 2) continue;
    lTotal++;
    if(lId != lLast) lCount++; 
    lLast = lId; 
  }
  std::cout <<" total " << lTotal << " -- count -- " << lCount << std::endl;
  return float(lTotal)/float(lCount);
}
TH2F** scan2D(std::string iLabel,TTree *iBTree,TTree *iSTree,double iRate,int iNBins,double iMin,double iMax,int iPtBins,double iPtMin,double iPtMax,float iScale=1.0) {
  double pIterX = (iMax-iMin)/double(iNBins);
  double pIterY = (iPtMax-iPtMin)/double(iPtBins);
  TH2F* lH2D    = new TH2F(("A"+iLabel).c_str(),("A"+iLabel).c_str(),iNBins+1,iMin-pIterX*0.5,iMax+pIterX*0.5,iPtBins+1,iPtMin-pIterY*0.5,iPtMax+pIterY*0.5);
  TH2F* lH2DEff = new TH2F(("E"+iLabel).c_str(),("E"+iLabel).c_str(),iNBins+1,iMin-pIterX*0.5,iMax+pIterX*0.5,iPtBins+1,iPtMin-pIterY*0.5,iPtMax+pIterY*0.5);
  float lUniqueScale = unique(iBTree);
  int lTotal  = 100000.;//*lUniqueScale;//iBTree->GetEntries();
  int lSTotal = iSTree->GetEntries("genpt > 25 && abs(geneta) < 2.5");
  std::stringstream lScale; lScale << "*" << iScale;
  for(int i0 = 0; i0 < lH2D->GetNbinsX()+1; i0++) {
    std::stringstream pCut; 
    double pCVal = lH2D->GetXaxis()->GetBinCenter(i0);
    if(iMax > 1.1) pCut << "fullIso    < " << pCVal;
    if(iMax < 1.1) pCut << "chargedIso > " << pCVal;
    std::stringstream pName; pName << iLabel << "tmp" << i0;
    TH1F *pH = new TH1F(pName.str().c_str(),pName.str().c_str(),iPtBins+1,iPtMin-pIterY*0.5,iPtMax+pIterY*0.5);
    iBTree->Draw(("pt"+lScale.str()+">>"+pName.str()).c_str(),pCut.str().c_str());
    double pVal = -1;
    for(int i1 = 0; i1 < pH->GetNbinsX()+1; i1++) {
      int pInt = pH->Integral(i1,1000000);
      float frac = 40000000.*(float(pInt)/float(lTotal));
      pVal = pH->GetBinCenter(i1);
      //if(frac < iRate) break;
      lH2D->SetBinContent(i0,i1,frac/1000.);
      if(i1 == 1) std::cout <<" ==> " << i1 << " -- " << pCut.str() << " -- " << pH->GetBinCenter(i1) << " -- " << pInt << " -- " << lTotal << " -- " << float(pInt)/float(lTotal) << " -- " << frac << std::endl;
      std::stringstream pTotCut;
      pTotCut << pCut.str() << " && pt > " << pVal << " && genpt > 25 && abs(geneta) < 2.5";
      int lSPass = iSTree->GetEntries(pTotCut.str().c_str());
      //std::cout << "--> cut " << pVal << " -- eff " << double(lSPass) << " -- " << double(lSTotal) << " -- " << pTotCut.str() << std::endl;
      double lEff = double(lSPass)/double(lSTotal);
      lH2DEff->SetBinContent(i0,i1,lEff);
    }
    //delete pH;
    iRate=iRate;
  }
  TH2F** lH2Ds = new TH2F*[2];
  lH2Ds[0] = lH2D;
  lH2Ds[1] = lH2DEff;
  lH2Ds[0]->GetZaxis()->SetTitle("Rate (kHz)");
  lH2Ds[1]->GetZaxis()->SetTitle("Signal Eff (p_{T}^{Gen} > 25 |#eta^{Gen}| < 2.5)");
  for(int i0 = 0; i0 < 2; i0++) lH2Ds[i0]->GetXaxis()->SetTitle("Disc)");
  for(int i0 = 0; i0 < 2; i0++) lH2Ds[i0]->GetYaxis()->SetTitle("p_{T} Scaled (GeV)");
  return lH2Ds;
}
TH1F** scan1D(std::string iLabel,std::string iIso,float iScale,TTree *iBTree,TTree *iSTree,double iRate,int iNBins,double iMin,double iMax,int iPtBins,double iPtMin,double iPtMax,int iColor) {
  double pIterX = (iMax-iMin)/double(iNBins);
  double pIterY = (iPtMax-iPtMin)/double(iPtBins);
  TH1F* lH2D    = new TH1F(("A"+iLabel).c_str(),("A"+iLabel).c_str(),iPtBins+1,iPtMin-pIterY*0.5,iPtMax+pIterY*0.5);
  TH1F* lH2DEff = new TH1F(("E"+iLabel).c_str(),("E"+iLabel).c_str(),iPtBins+1,iPtMin-pIterY*0.5,iPtMax+pIterY*0.5);
  float lUniqueScale = unique(iBTree);
  int lTotal  = 100000.;//*lUniqueScale;//iBTree->GetEntries();
  int lSTotal = iSTree->GetEntries("genpt > 25 && abs(geneta) < 2.5");
  std::cout <<" ===> " << lTotal << " -- " << lUniqueScale << std::endl;
  for(int i0 = 0; i0 < iPtBins; i0++) {
    if(i0 % 10 == 0) std::cout << "===> " << i0 << std::endl;
    std::stringstream pCut; pCut << "abs(eta) < 2.5 && pt*" << iScale;
    double pCVal = (pIterY*double(i0));
    pCut << " > " << pCVal;
    std::stringstream pName; pName << iLabel << "tmp" << i0;
    TH1F *pH = new TH1F(pName.str().c_str(),pName.str().c_str(),iNBins+1,iMin-pIterX*0.5,iMax+pIterX*0.5);
    iBTree->Draw((iIso+">>"+pName.str()).c_str(),pCut.str().c_str());
    double pVal = -1;
    double frac =  0;
    for(int i1 = 0; i1 < pH->GetNbinsX()+2; i1++) {
      int pInt = 0;
      pVal = pH->GetBinCenter(iNBins-i1);
      if(iMax < 1.1) pVal = pH->GetBinCenter(i1);
      if(iMax < 1.1) pInt = pH->Integral(i1,1000000);
      if(iMax > 1.1) pInt = pH->Integral(-1,iNBins-i1); 
      frac = 40000000.*(float(pInt)/float(lTotal));
      if(frac < iRate) break;
    }
    delete pH;
    lH2D->SetBinContent(i0,pVal);
    lH2D->SetBinError(i0,pIterX);
    std::stringstream pTotCut;
    if(iMax < 1.1) pTotCut << pCut.str() << " && "+iIso+" > " << pVal << " && genpt > 25 && abs(geneta) < 2.5";
    if(iMax > 1.1) pTotCut << pCut.str() << " && "+iIso+" < " << pVal << " && genpt > 25 && abs(geneta) < 2.5";
    int lSPass = iSTree->GetEntries(pTotCut.str().c_str());
    double lEff = double(lSPass)/double(lSTotal);
    lH2DEff->SetBinContent(i0,lEff);
    if(lSPass > 0) lH2DEff->SetBinError(i0,lEff*(1-lEff)/sqrt(lSTotal));
  }
  TH1F** lH2Ds = new TH1F*[2];
  lH2Ds[0] = lH2D;
  lH2Ds[1] = lH2DEff;
  lH2Ds[0]->GetYaxis()->SetTitle("Iso Cut");
  lH2Ds[1]->GetYaxis()->SetTitle("Signal Eff(p^{Gen}_{T} > 25 GeV |#eta| < 2.5)");
  for(int i0 = 0; i0 < 2; i0++) lH2Ds[i0]->GetXaxis()->SetTitle("p_{T} Scaled (GeV)");
  for(int i0 = 0; i0 < 2; i0++) lH2Ds[i0]->SetLineColor(iColor);
  for(int i0 = 0; i0 < 2; i0++) lH2Ds[i0]->SetMarkerColor(iColor);
  for(int i0 = 0; i0 < 2; i0++) lH2Ds[i0]->SetMarkerStyle(kFullCircle);
  return lH2Ds;
}
int fColor[6] = {kBlue,kGreen+1,kRed,kOrange,kMagenta+2,kCyan+2};
void plotRateOpt(TTree **iTree,float *iScale) { 
  TH1F** lNN   =  scan1D("NN"  ,"chargedIso",iScale[2],iTree[2],iTree[6],50000.,500,0.,1.0,100,5,100,fColor[2]);
  TH1F** lNNP   = scan1D("NNP" ,"chargedIso",iScale[3],iTree[3],iTree[7],50000.,500,0.,1.0,100,5,100,fColor[3]);
  TH1F** lPF    = scan1D("PF"  ,"chargedIso",iScale[0],iTree[0],iTree[4],50000.,500,0.,150,100,5,100,fColor[0]);
  TH1F** lPup   = scan1D("PUP" ,"chargedIso",iScale[1],iTree[1],iTree[5],50000.,500,0.,150,100,5,100,fColor[1]);
  TH1F** lPFF   = scan1D("PFF" ,"fullIso"   ,iScale[0],iTree[0],iTree[4],50000.,500,0.,150,100,5,100,fColor[4]);
  TH1F** lPupF  = scan1D("PUPF","fullIso"   ,iScale[1],iTree[1],iTree[5],50000.,500,0.,150,100,5,100,fColor[5]);
  lPF[0]  ->Draw();
  lPup[0] ->Draw("sames");
  lNN[0]  ->Draw("sames");
  lNNP[0] ->Draw("sames");
  lPFF[0] ->Draw("sames");
  lPupF[0]->Draw("sames");
 
  TCanvas *lC1 = new TCanvas("B","B",800,600); 
  lPF[1]  ->Draw();
  lPup[1] ->Draw("sames");
  lPFF[1] ->Draw("sames");
  lPupF[1]->Draw("sames");
  lNN[1]  ->Draw("sames");
  lNNP[1] ->Draw("sames");

  TLegend *lL = new TLegend(0.45,0.7,0.65,0.9); lL->SetFillColor(0); 
  lL->SetBorderSize(0); 
  lL->AddEntry(lPF[1]  ,"PF","lp");
  lL->AddEntry(lPup[1] ,"Puppi","lp");
  lL->AddEntry(lPFF[1] ,"PF Full","lp");
  lL->AddEntry(lPupF[1],"Puppi Full","lp");
  lL->AddEntry(lNN[1]  ,"NN+PF","lp");
  lL->AddEntry(lNNP[1] ,"NN+Puppi","lp");
  lL->Draw();
  return;
}
void plotScale(int iN, TTree **iTree,float *iScale) { 
  TH1F **lH = new TH1F*[iN*3];
  for(int i0 = 4; i0 < iN; i0++) { 
    iScale[i0] = 1.;
    std::stringstream pVar; pVar << iScale[i0] << "*pt/genpt";
    lH[0+3*i0] = plot("Scale",pVar.str(),0+iN*i0,"&& pt > 20 && genpt > 10"                                    ,fColor[i0 % 4],40,0,3.0,iTree[i0]); 
    lH[1+3*i0] = plot("Scale",pVar.str(),1+iN*i0,"&& pt > 20 && genpt > 10 && abs(eta) < 1.5"                  ,fColor[i0 % 4],40,0,3.0,iTree[i0]); 
    lH[2+3*i0] = plot("Scale",pVar.str(),2+iN*i0,"&& pt > 20 && genpt > 10 && abs(eta) > 1.5 && abs(eta) < 2.5",fColor[i0 % 4],40,0,3.0,iTree[i0]); 
  }
  for(int i0 = 12; i0 < 24; i0++) lH[i0]->Scale(1./lH[i0]->Integral());
  lH[12+0]->GetXaxis()->SetTitle("p_{T}/p_{T}^{Gen}");
  lH[12+0]->GetYaxis()->SetTitle("Normalized");
  lH[12+0]->Draw();
  lH[12+3]->Draw("sames");
  lH[12+6]->Draw("sames");
  lH[12+9]->Draw("sames");
  std::cout << "==> " << lH[12]->GetMean() << " -- " << lH[15]->GetMean() << " -- " << lH[18]->GetMean() << " -- " << lH[21]->GetMean() << std::endl;

 TLegend *lL = new TLegend(0.45,0.7,0.65,0.9); lL->SetFillColor(0); 
 lL->SetBorderSize(0); 
 lL->AddEntry(lH[12+0] ,"PF"      ,"lp");
 lL->AddEntry(lH[12+3] ,"Puppi"   ,"lp");
 lL->AddEntry(lH[12+6] ,"NN"      ,"lp");
 lL->AddEntry(lH[12+9],"NN+Puppi","lp");
 lL->Draw();
}
void plotIso(int iN,TTree **iTree,float *iScale) { 
  TH1F **lH = new TH1F*[iN*3];
  for(int i0 = 0; i0 < iN; i0++) { 
    iScale[i0] = 1.;
    std::stringstream pVar; pVar << "chargedIso";
    float lMax = 30.;
    if(i0 % 4 > 1) lMax = 1.;
    std::string lSig = " && genpt > 1";
    std::string lCut = " && pt > 20";
    if(i0 > 3) lCut += lSig;
    lH[0+3*i0] = plot("Scale",pVar.str(),0+iN*i0,lCut                                       ,fColor[i0 % 4],20,0,lMax,iTree[i0]); 
    lH[1+3*i0] = plot("Scale",pVar.str(),1+iN*i0,lCut+" && abs(eta) < 1.5"                  ,fColor[i0 % 4],20,0,lMax,iTree[i0]); 
    lH[2+3*i0] = plot("Scale",pVar.str(),2+iN*i0,lCut+" && abs(eta) > 1.5 && abs(eta) < 2.5",fColor[i0 % 4],20,0,lMax,iTree[i0]); 
    if(i0 < 4) for(int i1 = 0; i1 < 3; i1++) lH[i1+i0*3]->SetLineStyle(kDashed);
    for(int i1 = 0; i1 < 3; i1++) lH[i1+i0*3]->Scale(1./lH[i1+i0*3]->Integral());
  }
  lH[0]->GetXaxis()->SetTitle("Iso (GeV)");
  lH[0]->Draw();
  lH[3]->Draw("sames");
  lH[12+0]->Draw("sames");
  lH[12+3]->Draw("sames");
  TLegend *lL = new TLegend(0.45,0.7,0.65,0.9); lL->SetFillColor(0); 
  lL->SetBorderSize(0); 
  lL->AddEntry(lH[12+0]    ,"PF Sig"      ,"lp");
  lL->AddEntry(lH[12+3]    ,"Puppi Sig"   ,"lp");
  lL->AddEntry(lH[0] ,"PF Bkg"      ,"lp");
  lL->AddEntry(lH[3] ,"Puppi Bkg"   ,"lp");
  lL->Draw();  

  TCanvas *lC1 = new TCanvas("B","B",800,600);
  lH[0+7]->GetXaxis()->SetTitle("NN");
  lH[0+7]  ->Draw();
  lH[0+10] ->Draw("sames");
  lH[12+7] ->Draw("sames");
  lH[12+10]->Draw("sames");

  TLegend *lL0 = new TLegend(0.45,0.7,0.65,0.9); lL0->SetFillColor(0); 
  lL0->SetBorderSize(0); 
  lL0->AddEntry(lH[12+6] ,"NN Sig"      ,"lp");
  lL0->AddEntry(lH[12+9] ,"NN+Puppi Sig","lp");
  lL0->AddEntry(lH[0+6],"NN Bkg"      ,"lp");
  lL0->AddEntry(lH[0+9],"NN+Puppi Bkg","lp");
  lL0->Draw();
}
void plotPt(int iN,TTree **iTree,float *iScale) { 
  TH1F **lH = new TH1F*[iN*3];
  for(int i0 = 0; i0 < iN; i0++) { 
    //iScale[i0] = 1.;
    std::stringstream pVar; pVar << "pt" << "*" << iScale[i0];
    float lMax = 200.;
    std::string lSig = " && genpt > 1";
    std::string lCut = " && pt > 0";
    if(i0 > 3) lCut += lSig;
    lH[0+3*i0] = plot("Scale",pVar.str(),0+iN*i0,lCut                                       ,fColor[i0 % 4],40,0,lMax,iTree[i0]); 
    lH[1+3*i0] = plot("Scale",pVar.str(),1+iN*i0,lCut+" && abs(eta) < 1.5"                  ,fColor[i0 % 4],40,0,lMax,iTree[i0]); 
    lH[2+3*i0] = plot("Scale",pVar.str(),2+iN*i0,lCut+" && abs(eta) > 1.5 && abs(eta) < 2.5",fColor[i0 % 4],40,0,lMax,iTree[i0]); 
    if(i0 < 4) for(int i1 = 0; i1 < 3; i1++) lH[i1+i0*3]->SetLineStyle(kDashed);
  }
  //for(int i0 = 0; i0 < 12; i0++) lH[i0]->Scale(lH[i0+12]->Integral()/lH[i0]->Integral());
  for(int i0 = 0; i0 < 12; i0++) lH[i0]->Scale(0.1);
  lH[0]->GetXaxis()->SetTitle("p_{T}");
  lH[0]->Draw();
  lH[3]->Draw("sames");
  lH[12+0]->Draw("sames");
  lH[12+3]->Draw("sames");
  TLegend *lL = new TLegend(0.45,0.7,0.65,0.9); lL->SetFillColor(0); 
  lL->SetBorderSize(0); 
  lL->AddEntry(lH[12+0]    ,"PF"      ,"lp");
  lL->AddEntry(lH[12+3]    ,"Puppi"   ,"lp");
  lH[6]    ->Draw("sames");
  lH[9]    ->Draw("sames");
  lH[12+6] ->Draw("sames");
  lH[12+9] ->Draw("sames");
  lL->AddEntry(lH[12+6] ,"NN "      ,"lp");
  lL->AddEntry(lH[12+9] ,"NN+Puppi ","lp");
  lL->Draw();
}
void plotEffFake(int iN,TTree **iTree,float *iScale,std::string *iCuts) { 
  TGraphAsymmErrors **lH = new TGraphAsymmErrors*[iN*3];
  for(int i0 = 4; i0 < iN; i0++) { 
    lH[0+3*i0] = plotEff("Eff",0+iN*i0,""                                      ,iCuts[i0],fColor[i0 % 4],iTree[i0],iScale[i0]); 
    lH[1+3*i0] = plotEff("Eff",1+iN*i0,"&& abs(geneta) < 1.5"                  ,iCuts[i0],fColor[i0 % 4],iTree[i0],iScale[i0]); 
    lH[2+3*i0] = plotEff("Eff",2+iN*i0,"&& abs(geneta) > 1.5 && abs(geneta) < 2.5",iCuts[i0],fColor[i0 % 4],iTree[i0],iScale[i0]); 
  }
  lH[12+0+0]->Draw("ap");
  lH[12+3+0]->Draw("p sames");
  lH[12+6+0]->Draw("p sames");
  lH[12+9+0]->Draw("p sames");
  TLegend *lL = new TLegend(0.45,0.7,0.65,0.9); lL->SetFillColor(0); 
  lL->SetBorderSize(0); 
  lL->AddEntry(lH[12+0] ,"PF"      ,"lp");
  lL->AddEntry(lH[12+3] ,"Puppi"   ,"lp");
  lL->AddEntry(lH[12+6] ,"NN"      ,"lp");
  lL->AddEntry(lH[12+9],"NN+Puppi","lp");
  lL->Draw();

  TCanvas *lC2 = new TCanvas("E","E",800,600);
  TGraphAsymmErrors **lHE = new TGraphAsymmErrors*[iN*3];
  for(int i0 = 4; i0 < iN; i0++) { 
    lHE[0+3*i0] = plotEffEta("EEff",0+iN*i0,""                                      ,iCuts[i0],fColor[i0 % 4],iTree[i0],iScale[i0]); 
    lHE[1+3*i0] = plotEffEta("EEff",1+iN*i0,"&& abs(geneta) < 1.5"                  ,iCuts[i0],fColor[i0 % 4],iTree[i0],iScale[i0]); 
    lHE[2+3*i0] = plotEffEta("EEff",2+iN*i0,"&& abs(geneta) > 1.5 && abs(geneta) < 2.5",iCuts[i0],fColor[i0 % 4],iTree[i0],iScale[i0]); 
  }
  lHE[12+0+0]->Draw("ap");
  lHE[12+3+0]->Draw("p sames");
  lHE[12+6+0]->Draw("p sames");
  lHE[12+9+0]->Draw("p sames");
  TLegend *lLE = new TLegend(0.45,0.7,0.65,0.9); lL->SetFillColor(0); 
  lLE->SetBorderSize(0); 
  lLE->AddEntry(lHE[12+0] ,"PF"      ,"lp");
  lLE->AddEntry(lHE[12+3] ,"Puppi"   ,"lp");
  lLE->AddEntry(lHE[12+6] ,"NN"      ,"lp");
  lLE->AddEntry(lHE[12+9],"NN+Puppi","lp");
  lLE->Draw();

  TH1F **lH2 = new TH1F*[iN*3];
  TCanvas *lC1 = new TCanvas("B","B",800,600);
  for(int i0 = 0; i0 < 4; i0++) { 
    lH2[0+3*i0] = plotRate("Rate",0+iN*i0,""                                   ,iCuts[i0],fColor[i0 % 4],iTree[i0],iScale[i0]); 
    lH2[1+3*i0] = plotRate("Rate",1+iN*i0,"&& abs(eta) < 1.5"                  ,iCuts[i0],fColor[i0 % 4],iTree[i0],iScale[i0]); 
    lH2[2+3*i0] = plotRate("Rate",2+iN*i0,"&& abs(eta) > 1.5 && abs(eta) < 2.5",iCuts[i0],fColor[i0 % 4],iTree[i0],iScale[i0]); 
  }
  lH2[0+0]->Draw();
  lH2[3+0]->Draw("sames");
  lH2[6+0]->Draw("sames");
  lH2[9+0]->Draw("sames");
  TLegend *lL0 = new TLegend(0.45,0.7,0.65,0.9); lL0->SetFillColor(0); 
  lL0->SetBorderSize(0); 
  lL0->AddEntry(lH2[0],"PF"      ,"lp");
  lL0->AddEntry(lH2[3],"Puppi"   ,"lp");
  lL0->AddEntry(lH2[6],"NN"      ,"lp");
  lL0->AddEntry(lH2[9],"NN+Puppi","lp");
  lL0->Draw();
  return;
}
void plot2D(TTree **iTree,float *iScale) { 
  //TH2F** lPF = scan2D("PF",iTree[0],iTree[4],10000.,100,0.,30.,25,10,75,iScale[0]);
  //lPF[0]->Draw("colz");
  //TCanvas *lC1 = new TCanvas("B","B",800,600); 
  //lPF[1]->Draw("colz");
  
  //TCanvas *lC2 = new TCanvas("C","C",800,600); 
  TH2F** lNN = scan2D("NN",iTree[3],iTree[7],10000.,100,0.,1.0,50,0,75,iScale[2]);
  lNN[0]->Draw("colz");
  TCanvas *lC3 = new TCanvas("D","D",800,600); 
  lNN[1]->Draw("colz");
}

void quickPlot() { 
  int lN = 8;
  TTree **lTree = new TTree*[lN];
  TFile *l0File  = new TFile("MB_flat_emcand_scale2_10_dZ10.root");//MB.root");
  lTree[0]      = (TTree*) l0File->Get("ntuplePF/tree");
  lTree[1]      = (TTree*) l0File->Get("ntuplePF2/tree");
  lTree[2]      = (TTree*) l0File->Get("ntupleNN/tree");
  //TFile *l2File  = new TFile("MB.root");
  lTree[3]      = (TTree*) l0File->Get("ntupleNN2/tree");

  TFile *l1File  = new TFile("Sig_flat_emcand_scale2_10_dZ10.root");//Sig_flat_emcand_scale_10.root");//Sig_flat_emcand_scale_10_dZ10.root ");//Sig_flat_emcand_scale_10.root");//_flat_emcand.root");//Sig_flat_Invisble100.root");//Sig_flat_emcand.root");
  lTree[4]      = (TTree*) l1File->Get("ntuplePF/tree");
  lTree[5]      = (TTree*) l1File->Get("ntuplePF2/tree");
  lTree[6]      = (TTree*) l1File->Get("ntupleNN/tree");
  //TFile *l3File  = new TFile("Sig.root");
  lTree[7]      = (TTree*) l1File->Get("ntupleNN2/tree");

  TCanvas *lC0 = new TCanvas("A","A",800,600); 
  std::string *lArr = new std::string[8]; 
  lArr[0] = " && fullIso < 2000";
  lArr[1] = " && fullIso < 2000";
  lArr[2] = " && chargedIso > -0.01";
  lArr[3] = " && chargedIso > -0.01";
  lArr[4] = " && fullIso < 2000";
  lArr[5] = " && fullIso < 2000";
  lArr[6] = " && chargedIso > -0.01";
  lArr[7] = " && chargedIso > -0.01";
  float *lScale = new float[8]; 
  lScale[0] = 1./1.25;//1.2;
  lScale[1] = 1./1.20;//1.2;
  lScale[2] = 1.;
  lScale[3] = 1.;
  lScale[4] = 1./1.25;//1.2;
  lScale[5] = 1./1.20;//1.2;
  lScale[6] = 1.;
  lScale[7] = 1.;
  //plotScale(lN,lTree,lScale);
  //plotIso(lN,lTree,lScale);
  //plotPt(lN,lTree,lScale);
  plotEffFake(lN,lTree,lScale,lArr);
  //plot2D(lTree,lScale);
  //plotRateOpt(lTree,lScale);
  return; 
}

