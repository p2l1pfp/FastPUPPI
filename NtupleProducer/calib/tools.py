#! /usr/bin/env python
import ROOT as r,sys,math,array,os
from optparse import OptionParser
from ROOT import std,RooDataHist

import tdrstyle
tdrstyle.setTDRStyle()

#parser = OptionParser()
#parser.add_option('--input'  ,action='store',type='string',dest='input'    ,default='mlfit.root',help='input file')
#parser.add_option('--data'   ,action='store',type='string',dest='data'     ,default='combine_Zprime150_0.root',help='data file')
#parser.add_option('--cats'    ,action='store',type='string',dest='cats'    ,default='1',help='categories 0,1,2,...')
#parser.add_option('--passfail',action='store_true',         dest='passfail',default=False,help='fit pass and failng') 
#parser.add_option('--divide'  ,action='store_true',         dest='divide'  ,default=False,help='ratio') 
#parser.add_option('--mass'    ,action='store',              dest='mass'    ,default=100  ,help='mass') 

#(options,args) = parser.parse_args()

def end():
    rep = ''
    while not rep in [ 'q', 'Q','a',' ' ]:
        rep = raw_input( 'enter "q" to quit: ' )
        if 1 < len(rep):
            rep = rep[0]
            
def divide(iData,iHists):
    lPass = -1
    lFail = -1
    for i0 in range(0,len(iHists)):
        iHists[i0].Sumw2()
        if iHists[i0].GetName().find("") > -1 and lPass == -1:
            lPass = i0
        if iHists[i0].GetName().find("fail") > -1 and lFail == -1:
            lFail = i0
    if lPass != -1:
        iData[0].Divide(iHists[lPass])
        for i0 in range(0,iData[0].GetNbinsX()+1):
            if iHists[lPass].GetBinContent(i0) > 0:
                iData[0].SetBinError(i0,iData[0].GetBinError(i0)/iHists[lPass].GetBinContent(i0))
        for i0 in range(0,len(iHists)):
            if i0 != lPass and iHists[i0].GetName().find("") > -1:
                iHists[i0].Divide(iHists[lPass])
        iHists[lPass].Divide(iHists[lPass])

    if lFail != -1:
        iData[1].Divide(iHists[lFail])
        for i0 in range(0,len(iHists)):
            if i0 != lFail and iHists[i0].GetName().find("fail") > -1:
                iHists[i0].Divide(iHists[lFail])
        iHists[lFail].Divide(iHists[lFail])

def draw(iData,iHists,iName="A",iPF=True,iDivide=False):
    if iDivide:
        divide(iData,iHists)

    lC0 = r.TCanvas("Can"+iName,"Can"+iName,800,600);
    if iPF:
        lC0.Divide(2)
    lC0.cd(1)
    lHist=fix(iHists[0],iData[0])
    lHist.Draw("ep")
    pDraw=False
    for pHist in iHists:
        pHist.Draw("hist sames") if pDraw else pHist.Draw("e2 sames")
        pDraw=True
    lHist.Draw("ep sames")
    lLegend = r.TLegend(0.63,0.63,0.88,0.88)
    lLegend.SetFillColor(0)
    lLegend.SetBorderSize(0)
    lLegend.AddEntry(iData [0],"data","lp")
    lLegend.AddEntry(iHists[0],"QCD","lf")
    lLegend.AddEntry(iHists[1],"W#rightarrow qq","l")
    lLegend.AddEntry(iHists[2],"Z#rightarrow qq","l")
    #lLegend.AddEntry(iHists[3],"Signal","lf")
    lLegend.Draw()
    if not iPF:
        lC0.SaveAs(iName+".png")
        lC0.SaveAs(iName+".pdf")
        end()
        return

    lC0.cd(2)
    iData[1].Draw("ep")
    pDraw=False
    for pHist in iHists:
        if pHist.GetName().find("fail") > -1:
            pHist.Draw("hist sames") if pDraw else pHist.Draw("e2 sames")
            pDraw=True
    lC0.Modified()
    lC0.Update()
    iData[1].Draw("ep sames")
    lC0.SaveAs(iName+".png")
    lC0.SaveAs(iName+".pdf")
    end()

def norm(iFile,iH,iName) :
    lNorms = iName.split("/")[0].replace("shapes","norm")
    print "---",lNorms
    lArg = iFile.Get(lNorms)
    lName = iName.replace(iName.split("/")[0]+"/","")
    lVal = r.RooRealVar(lArg.find(lName)).getVal();
    print "!!!",lName,lNorms,iH.Integral(),lVal
    iH.Scale(lVal/iH.Integral())

def load(iFile,iName):
    lHist = iFile.Get(iName)
    norm(iFile,lHist,iName)
    lHist.SetName(iName.replace("/","_"))
    lHist.SetTitle(iName.replace("/","_"))
    return lHist

def fix(iH0,iH1):
    lTHist = r.TH1F(iH1.GetName(),iH1.GetName(),iH0.GetNbinsX(),iH0.GetXaxis().GetXmin(),iH0.GetXaxis().GetXmax())
    for i0 in range(0,iH0.GetNbinsX()+1):
        lTHist.SetBinContent(i0,iH1.GetBinContent(i0))
        lTHist.SetBinError(i0,iH1.GetBinError(i0))
    lTHist.SetFillStyle(iH1.GetFillStyle())
    lTHist.SetFillColor(iH1.GetFillColor())
    lTHist.SetLineColor(iH1.GetLineColor())
    lTHist.SetLineStyle(iH1.GetLineStyle())
    return lTHist
    
def loadHist(iFile,iCat,iMass,iEWK=True,iS=True):
    lHists = []
    lFit = "shapes_fit_s/"+iCat+"/" if iS else "shapes_fit_b/"+iCat+"/"
    #lFit = "shapes_prefit/"+iCat+"/" if iS else "shapes_fit_b/"+iCat+"/"
    lHists.append(load(iFile,lFit+"qcd"))
    if iEWK:
        lId = len(lHists)
        lHists.append(load(iFile,lFit+"Winc"))
        lHists.append(load(iFile,lFit+"Zinc"))
        if iS:
            lHists.append(load(iFile,lFit+"sig"))
        lHists[lId].SetLineColor(46)
        lHists[lId+1].SetLineColor(9)
        if iS:
            lHists[lId+2].SetLineColor(8)
        lHists[lId].SetLineWidth(2)
        lHists[lId+1].SetLineWidth(2)
        if iS :
            lHists[lId+2].SetLineWidth(2)
        lHists[lId].Add(lHists[0])
        lHists[lId+1].Add(lHists[1])
        if iS :
            lHists[lId+2].Add(lHists[2])
    lHists[0].SetFillColor(16)
    lHists[0].SetFillStyle(3001)
    lHists[0].SetLineStyle(2)
    return lHists

def loadData(iDataFile,iCat):
    lData = iDataFile.Get("data_obs")
    lData.GetXaxis().SetTitle("m_{J} (GeV)")
    lData.SetMarkerStyle(20)
    lData.Sumw2()
    return [lData]

if __name__ == "__main__":
    lHFile = r.TFile(options.input)
    lDFile = r.TFile(options.data)
    lDSum=[]
    lSum=[]
    for cat in options.cats.split(','):
        lData  = loadData(lDFile,"pass_cat"+cat)
        lHists = loadHist(lHFile,"ch0",options.mass)
        if options.passfail:
            lData .extend(loadData(lDFile,"fail_cat"+cat))
            lHists.extend(loadHist(lHFile,"ch"+cat+"_fail_cat"+cat,options.mass,True,True))
        if len(lSum) == 0:
            lDSum = lData
            lSum  = lHists
        else:
            for i0 in range(0,len(lDSum)):
                lDSum[i0].Add(lData[i0])
            for i0 in range(0,len(lSum)):
                lSum[i0].Add(lHists[i0])
        draw(lData,lHists,cat,options.passfail,options.divide)
    draw(lDSum,lSum,"sum",options.passfail,options.divide)

