//#include "NYStyle.h"
#include "tdrstyle.C"
#include <sstream>

TH1F   *fH ; 
TGraphErrors *fG0;
TGraphErrors *fG1;
TF1    *fF1;
TF1    *fF2;
bool fCorrect = false;

void formatCanvas(TCanvas *iC0) { 
  TLatex *lat =  new TLatex(); 
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->SetTextSize(0.04);
  lat->DrawLatex(0.25,0.88,"#bf{CMS Simulation} ");
  lat->DrawLatex(0.25,0.83,"#it{Preliminary}");
  //lat->Draw();

  TLatex *lat2 =  new TLatex(); 
  lat2->SetNDC();
  lat2->SetTextFont(42);
  lat2->SetTextSize(0.04);
  lat2->DrawLatex(0.75,0.96,"13 TeV, PU = 150");

}
TH1F* Xratio(TH1F *iRatio,TH1F *iBase,bool iSame=false) { 
  TH1F* lRatio = (TH1F*) iRatio->Clone("");
  lRatio->Divide(iBase);
  lRatio->SetMarkerColor(lRatio->GetLineColor());
  //lRatio->SetLineColor(kBlack);
  lRatio->SetTitle("");
  lRatio->SetMarkerStyle(kFullCircle);
  lRatio->GetXaxis()->SetTitleSize(10);
  lRatio->GetYaxis()->SetTitleSize(0.08);
  lRatio->GetYaxis()->SetTitleOffset(0.6);
  lRatio->GetYaxis()->SetTitle("X/PF");
  lRatio->GetYaxis()->SetRangeUser(0.001,10);
  lRatio->GetYaxis()->SetLabelSize(0.1);
  if(!iSame) lRatio->Draw("ep");
  if(iSame)  lRatio->Draw("ep sames");
  TLine *lLine = new TLine(5.,1.,125.,1.);
  lLine->SetLineStyle(kDashed);
  lLine->Draw();
  return lRatio;
}
TF1* check(TGraphErrors *iG) { 
  double *lX = new double[iG->GetN()+1];
  double lXOld = 0;
  for(int i0 = 0; i0 < iG->GetN(); i0++) { 
    lX[i0] = (iG->GetX()[i0]+lXOld)/2.;
    lXOld  = iG->GetX()[i0];
  }
  lX[iG->GetN()] = (iG->GetX()[iG->GetN()-1]-lX[iG->GetN()-1])+iG->GetX()[iG->GetN()-1];
  TH1F *lH = new TH1F("A","A",iG->GetN(),lX);
  for(int i0 = 0; i0 < iG->GetN(); i0++) { 
    lH->Fill(iG->GetX()[i0],iG->GetY()[i0]);
  }
  lH->Fit("pol3","","R",1,200);
  return lH->GetFunction("pol3");
}
void makeMETPlotsFile(std::string iVar="",std::string iMet="",std::string iName="Pup_vDZ03MinR02_A0_v2.root",
		      int iColor=1) { 
  TFile *lFile = new TFile(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("met");
  const int lN = 9;
  double lEnd = 750; 
  //double lRange[lN+1]  = {0,10,20,30,40,50,60,80,100,125,150,200};
  double lRange[lN+1]  = {0,10,20,30,40,60,80,100,150,200};
  //double lRange[lN+1]  = {0,10,20,30,40,50,60,80,100,125,200};
  double lX    [lN]; 
  double lEX   [lN]; 
  double lMean [lN]; 
  double lRMS  [lN]; 
  double lEMean[lN]; 
  double lERMS [lN]; 
  std::string lCut  = "(pt_z > 5 && m_z > 60 && m_z < 120 && abs(dzmu-dz) < 10.2)";//*("+iMet+"met > 0)";
  TH1F *lXH = new TH1F(("Met"+iVar).c_str(),("Met"+iVar).c_str(),35,0,530); lXH->Sumw2();
  lXH->GetXaxis()->SetTitle("#slash{E_{T}} (GeV)");
  lXH->SetLineColor(iColor);
  std::string lHDraw = iMet+"met>>Met"+iVar;
  std::cout << "===> " << lHDraw << " -- " << lCut << std::endl;
  lTree->Draw(lHDraw.c_str(),lCut.c_str());
  fH = lXH;
  fH->GetYaxis()->SetTitle("Events/5 GeV");

  TCanvas *lC0 = new TCanvas("A","A",800,600); 
  for(int i0 = 0; i0 < lN; i0++) { 
    TH1F *lH  = new TH1F(("Res" +iVar).c_str(),( "Res"+iVar).c_str(),100,-1000,1000);
    TH1F *lH0 = new TH1F(("HRes"+iVar).c_str(),("HRes"+iVar).c_str(),100,-2,3);
    //TH1F *lH0 = new TH1F(("HRes"+iVar).c_str(),("HRes"+iVar).c_str(),100,-1000,1000);
    std::string lDraw  = "("+iVar+"-pt_z)>>Res"+iVar;
    std::string lRDraw = "("+iVar+"/pt_z)>>HRes"+iVar;
    if(iVar.find("u2") != std::string::npos) lDraw = "("+iVar+")>>Res"+iVar;
    std::stringstream pSS; pSS << lCut << "*(pt_z > " << lRange[i0] << " && pt_z < " << lRange[i0+1] << ")";
    lTree->Draw(lDraw.c_str(),pSS.str().c_str());
    if(iVar.find("u2") == std::string::npos) lTree->Draw(lRDraw.c_str(),pSS.str().c_str());
    lH ->Fit("gaus");
    if(iVar.find("u2") == std::string::npos) lH0->Fit("gaus");
    lX    [i0] = (lRange[i0] + lRange[i0+1])/2.;
    lEX   [i0] = fabs(lRange[i0] - lRange[i0+1])/2.;//sqrt(12.);
    if(iVar.find("u2") != std::string::npos) lMean [i0] = lH->GetFunction("gaus")->GetParameter(1);
    if(iVar.find("u2") != std::string::npos) lEMean[i0] = lH->GetFunction("gaus")->GetParError (1);
    if(iVar.find("u2") == std::string::npos) lMean [i0] = lH0->GetMean();
    if(iVar.find("u2") == std::string::npos) lEMean[i0] = lH0->GetRMS()/sqrt(lTree->GetEntries(pSS.str().c_str()));//GetFunction("gaus")->GetParError (1);
    lRMS  [i0] = lH->GetFunction("gaus")->GetParameter(2);
    lERMS [i0] = lH->GetFunction("gaus")->GetParError (2);
    delete lH;
  }
  /*
  for(int i0 = 0; i0 < lN; i0++) { 
    lMean[i0]  = (1.+ lMean[i0]/lX[i0]);
    lEMean[i0] = lEMean[i0]/lX[i0];
  }
  */
  TGraphErrors *lG0 = new TGraphErrors(lN,lX,lMean,lEX,lEMean);
  lG0->SetMarkerStyle(kFullCircle);
  //lG0->Fit("pol2","R","",30,lEnd);
  lG0->SetLineColor  (iColor); 
  lG0->SetMarkerColor(iColor); 
  lG0->GetYaxis()->SetRangeUser(0.,1.5);
  lG0->GetXaxis()->SetTitle("p_{T}^{Z} (GeV)");
  lG0->GetYaxis()->SetTitle("<u_{#parallel}>/p_{T}");
  fG0 = lG0;
  if(iVar.find("u1") != std::string::npos) fF2 = check(fG0);
  
  if(fCorrect) { for(int i0 = 0; i0 < lN; i0++) lRMS[i0] *= 1./fF2->Eval(lX[i0]);}
  //if(fCorrect) { for(int i0 = 0; i0 < lN; i0++) lRMS[i0] *= 1./fG0->Eval(lX[i0]);}
  TGraphErrors *lG1 = new TGraphErrors(lN,lX,lRMS,lEX,lERMS);
  lG1->SetMarkerStyle(kFullCircle);
  TF1 *lF1 = new TF1("func","[0]+[1]/sqrt(x)+[2]/x",30,250);
  //TF1 *lF1 = new TF1("func","pol1",30,300);
  lF1->SetLineColor(iColor);
  //lG1->Fit("func","","R",30,lEnd);
  lG1->SetLineColor  (iColor); 
  lG1->SetMarkerColor(iColor);   
  lG1->GetYaxis()->SetRangeUser(0,400);
  if(iVar.find("pup")  != std::string::npos) lG1->GetYaxis()->SetRangeUser(0,55);
  fG1 = lG1;
  lG1->GetXaxis()->SetTitle("p_{T}^{Z} (GeV)");
  //lG1->GetXaxis()->SetTitle("n_{PV}");
  if(iVar.find("u2") != std::string::npos) lG1->GetYaxis()->SetTitle("#sigma(u_{#perp} ) (GeV)");
  if(iVar.find("u1") != std::string::npos) lG1->GetYaxis()->SetTitle("#sigma(u_{#parallel} ) (GeV)");
  if(fCorrect) { 
    if(iVar.find("u2") != std::string::npos) lG1->GetYaxis()->SetTitle("#sigma(u_{#perp})/<u_{#parallel}>/p_{T}");
    if(iVar.find("u1") != std::string::npos) lG1->GetYaxis()->SetTitle("#sigma(u_{#parallel})/<u_{#parallel}>/p_{T}");
  }
  fF1 = lF1;
}
void makePlots_v3(
		  std::string iFile0="ZMM_def_140_v30.root",
		  std::string iFile1="ZMM_def_140_v8.root",
		  std::string iFile2="ZMM_def_140_v30.root",
		  std::string iFile3="ZMM_def_140_v8.root"
	       ) {  
  //Prep();
  setTDRStyle();
  std::string lLabel="";
  std::string lLabel2="";
  std::string lLabel_b="";
  std::string lLabel2_b="";
  std::string lLast  ="";
  makeMETPlotsFile(lLabel+"u1"+lLast,"",iFile0,kRed);
  TH1F * lH0    = (TH1F*)fH ->Clone("hist0"); 
  TGraphErrors  * lPupF0 = (TGraphErrors*) fG0->Clone("Rep0"); 
  TGraphErrors  * lPupG0 = (TGraphErrors*) fG1->Clone("Res0"); 

  makeMETPlotsFile(lLabel+"u1"+lLast,"",iFile1,kGreen+1);
  TH1F * lH1    = (TH1F*)fH ->Clone("hist0"); 
  TGraphErrors  * lPupF1 = (TGraphErrors*) fG0->Clone("Rep0"); 
  TGraphErrors  * lPupG1 = (TGraphErrors*) fG1->Clone("Res0"); 

  makeMETPlotsFile(lLabel+"pupu1"+lLast,"pup",iFile2,kOrange-9);
  TH1F * lH2    = (TH1F*)fH ->Clone("hist0"); 
  TGraphErrors  * lPupF2 = (TGraphErrors*) fG0->Clone("Rep0"); 
  TGraphErrors  * lPupG2 = (TGraphErrors*) fG1->Clone("Res0"); 

  makeMETPlotsFile(lLabel+"pupu1"+lLast,"pup",iFile3,kBlue-9);
  TH1F * lH3    = (TH1F*)fH ->Clone("hist0"); 
  TGraphErrors  * lPupF3 = (TGraphErrors*) fG0->Clone("Rep0"); 
  TGraphErrors  * lPupG3 = (TGraphErrors*) fG1->Clone("Res0"); 

  makeMETPlotsFile(lLabel+"u2"+lLast,"",iFile0,kRed);
  TGraphErrors  * lPupF20 = (TGraphErrors*) fG0->Clone("Rep0"); 
  TGraphErrors  * lPupG20 = (TGraphErrors*) fG1->Clone("Res0"); 

  makeMETPlotsFile(lLabel+"u2"+lLast,"",iFile1,kGreen+1);
  TGraphErrors  * lPupF21 = (TGraphErrors*) fG0->Clone("Rep0"); 
  TGraphErrors  * lPupG21 = (TGraphErrors*) fG1->Clone("Res0"); 

  makeMETPlotsFile(lLabel+"pupu2"+lLast,"pup",iFile2,kOrange-9);
  TGraphErrors  * lPupF22 = (TGraphErrors*) fG0->Clone("Rep0"); 
  TGraphErrors  * lPupG22 = (TGraphErrors*) fG1->Clone("Res0"); 

  makeMETPlotsFile(lLabel+"pupu2"+lLast,"pup",iFile3,kBlue-9);
  TGraphErrors  * lPupF23 = (TGraphErrors*) fG0->Clone("Rep0"); 
  TGraphErrors  * lPupG23 = (TGraphErrors*) fG1->Clone("Res0"); 

  std::string lLeg1 = "Calo";
  std::string lLeg2 = "PF";//PF #Delta #eta";
  std::string lLeg3 = "TK";//PF #Delta #eta 075";
  std::string lLeg4 = "Pup";//PF #Delta #eta 25";

  //std::string lLeg1 = "4 GeV";
  //std::string lLeg2 = "3 GeV";
  //std::string lLeg3 = "2 GeV";
  //std::string lLeg4 = "1 GeV";
  TLegend *lL = new TLegend(0.65,0.6,0.6,0.88); lL->SetBorderSize(0); lL->SetFillColor(0); 
  lL->AddEntry(lH0,lLeg1.c_str(),"l");  
  lL->AddEntry(lH1,lLeg2.c_str(),"l");
  lL->AddEntry(lH2,lLeg3.c_str(),"l");  
  lL->AddEntry(lH3,lLeg4.c_str(),"l");

  TCanvas *lCX = new TCanvas("B","B",800,600);  
  lCX->Divide(1,2); lCX->cd();  lCX->cd(1)->SetPad(0,0.3,1.0,1.0); gPad->SetLeftMargin(0.2) ;
  lCX->cd(1)->SetLogy();
  lH1->Scale(lH0->Integral()/lH1->Integral());
  lH2->Scale(lH0->Integral()/lH2->Integral());
  lH3->Scale(lH0->Integral()/lH3->Integral());
  lH0->GetYaxis()->SetRangeUser(0.1,lH0->GetMaximum()*15.);
  lH0->Draw("hist");
  lH1->Draw("hist sames");
  lH2->Draw("hist sames");
  lH3->Draw("hist sames");
  lL->Draw();
  formatCanvas(lCX);
  lCX->cd(2)->SetPad(0,0,1.0,0.3); gPad->SetLeftMargin(0.2) ;
  Xratio(lH2,lH1);
  Xratio(lH3,lH1,true);
  //lCX->SaveAs("Met200.pdf");
  //lCX->SaveAs((iFile1+"Met200.png").c_str());
  //lCX->SaveAs("Met200.C");
 
  TLegend *lL1 = new TLegend(0.6,0.15,0.9,0.4); lL1->SetBorderSize(0); lL1->SetFillColor(0); 
  lL1->AddEntry(lPupF0,lLeg1.c_str(),"lp");
  lL1->AddEntry(lPupF1,lLeg2.c_str(),"lp");
  lL1->AddEntry(lPupF2,lLeg3.c_str(),"lp");
  lL1->AddEntry(lPupF3,lLeg4.c_str(),"lp");
  
  TCanvas *lC0 = new TCanvas("C","C",800,600); 
  lPupF0 ->Draw("ape");
  lPupF1 ->Draw("pe sames");
  lPupF2 ->Draw("pe sames");
  lPupF3 ->Draw("pe sames");
  TLine *lLine = new TLine(10,1,200,1); lLine->SetLineStyle(kDashed); lLine->Draw();
  lL1->Draw();
  formatCanvas(lC0);
  //lC0->SaveAs("Response200.pdf");
  //lC0->SaveAs((iFile1+"Response.png").c_str());
  //lC0->SaveAs("Response200.C");

  TLegend *lL2 = new TLegend(0.6,0.15,0.9,0.4); lL2->SetBorderSize(0); lL2->SetFillColor(0); 
  lL2->AddEntry(lPupG0,lLeg1.c_str(),"lp");  
  lL2->AddEntry(lPupG1,lLeg2.c_str(),"lp");
  lL2->AddEntry(lPupG2,lLeg3.c_str(),"lp");
  lL2->AddEntry(lPupG3,lLeg4.c_str(),"lp");

  TCanvas *lC1 = new TCanvas("D","D",800,600); 
  lPupG0 ->Draw("ape");
  lPupG1 ->Draw("pe sames");
  lPupG2 ->Draw("pe sames");
  lPupG3 ->Draw("pe sames");
  lL2->Draw();
  formatCanvas(lC1);
  //lC1->SaveAs((iFile1+"U1Res.png").c_str());
  //lC1->SaveAs("U1Res200.pdf");
  //lC1->SaveAs("U1Res200.C");

  TLegend *lL3 = new TLegend(0.6,0.15,0.9,0.4); lL3->SetBorderSize(0); lL3->SetFillColor(0); 
  lL3->AddEntry(lPupG20,lLeg1.c_str(),"lp");  
  lL3->AddEntry(lPupG21,lLeg2.c_str(),"lp");
  lL3->AddEntry(lPupG22,lLeg3.c_str(),"lp");
  lL3->AddEntry(lPupG23,lLeg4.c_str(),"lp");

  TCanvas *lC3 = new TCanvas("E","E",800,600); formatCanvas(lC3);
  lPupG20 ->Draw("ape");
  lPupG21 ->Draw("pe sames");
  lPupG22 ->Draw("pe sames");
  lPupG23 ->Draw("pe sames");
  lL3->Draw();
  formatCanvas(lC3);
  //lC3->SaveAs((iFile1+"U2Res.png").c_str());
  //lC3->SaveAs("U2Res200.pdf");
  //lC3->SaveAs("U2Res200.C");
  return;
  TFile *lFile = new TFile("Output.root","RECREATE");
  lPupG20->SetName("NoTKUPerpRes");
  lPupG21->SetName ("TkUPerpRes");
  lPupG0->SetName("NoTKUParRes");
  lPupG1->SetName ("TkUParaRes");
  lPupF0->SetName("NoTKUParResponse");
  lPupF1->SetName ("TkUParaResponse");
  lH0   ->SetName("MetNoTk");
  lH1   ->SetName("MetTk");

  lPupG20->Write();
  lPupG21->Write();
  lPupG0->Write();
  lPupG1->Write();
  lPupF0->Write();
  lPupF1->Write();
  lH0   ->Write();
  lH1   ->Write();

  
  //TLegend *lL4 = new TLegend(0.2,0.7,0.45,0.9); lL4->SetBorderSize(0); lL4->SetFillColor(0); 
  //lL4->AddEntry(lPupF20,"Track |#eta| < 2.5","lp");  
  //lL4->AddEntry(lPupF21,"Track |#eta| < 4.0","lp");

  //TCanvas *lC4 = new TCanvas("F","F",800,600); 
  //lPupF20 ->Draw("alpe");
  //lPupF21 ->Draw("lpe sames");
  //lL4->Draw();
  
  //check(lPupF1);
}
