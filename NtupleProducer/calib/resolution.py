#! /usr/bin/env python
import ROOT as r,sys,math,array,os,argparse,numpy
from array import array
from tools import end
import tdrstyle
tdrstyle.setTDRStyle()

parser = argparse.ArgumentParser(description='Process benchmarks.')
parser.add_argument('-Tree'   ,'--Tree'      ,action='store',dest='Tree',default='HcalInfo',help='Tree Name')
parser.add_argument('-input'  ,'--input'     ,action='store',dest='input'  ,default='Pion3.root',help='input file')
parser.add_argument('-var'    ,'--var'       ,action='store',dest='var'    ,default='hcal_corr_et'         ,help='mass') 
parser.add_argument('-genvar' ,'--genvar'    ,action='store',dest='genvar' ,default='genPt'         ,help='mass') 
parser.add_argument('-eta'    ,'--eta'       ,dest='etarange',nargs='+',type=float,default=[-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0],help='eta range')
parser.add_argument('-pt'     ,'--pt'        ,dest='ptrange' ,nargs='+',type=float,default=range(2,60,2),help='eta range')
#parser.add_argument('-frac'   ,'--frac'      ,dest='fracrange',nargs='+',type=float,default=[x * 0.1 for x in range(0, 11)],help='em fraction range')
parser.add_argument('-frac'   ,'--frac'      ,dest='fracrange',nargs='+',type=float,default=[0],help='em fraction range')

args = parser.parse_args()

fColor=[r.kRed,r.kBlue,r.kGreen+1,r.kViolet,r.kRed-2,r.kBlue-2,r.kSpring,r.kAzure+1,r.kOrange]

def draw(iGraphs,iEta,iName="A",iGaus=False):
    lC0 = r.TCanvas("Can"+iName,"Can"+iName,800,600);
    lLeg = r.TLegend(0.2,0.65,0.45,0.85)
    lLeg.SetBorderSize(0)
    lLeg.SetFillColor(0)
    #iGraphs[0].GetYaxis().SetRangeUser(-0.1,0.1)
    iGraphs[0].Draw("alep")
    lLeg.AddEntry(iGraphs[0],str(iEta[0])+"< #eta < "+str(iEta[0+1]),"lpe")
    iGraphs[0].SetName (str(iEta[0])+"< #eta < "+str(iEta[0+1])+iName)
    iGraphs[0].SetTitle(str(iEta[0])+"< #eta < "+str(iEta[0+1])+iName)
    pF1 = r.TF1("f1","[0]+[1]/sqrt(x)+[2]/x",0.1,200);
    if iGaus:
        iGraphs[0].Fit("f1")
    for i0 in range(1,len(iGraphs)):
        print i0,iGraphs,fColor[i0 % len(fColor)]
        iGraphs[i0].SetLineColor(fColor[i0 % len(fColor)])
        iGraphs[i0].SetMarkerColor(fColor[i0  % len(fColor)])
        iGraphs[i0].Draw("lpe")
        lLeg.AddEntry(iGraphs[i0],str(iEta[i0])+"< #eta < "+str(iEta[i0+1]),"lpe")
        iGraphs[i0].SetName(str(iEta[i0])+"< #eta < "+str(iEta[i0+1])+iName)
        iGraphs[i0].SetTitle(str(iEta[i0])+"< #eta < "+str(iEta[i0+1])+iName)
        pF1 = r.TF1("f1","[0]+[1]/sqrt(x)+[2]/x",0.1,200);
        #if iGaus:
            #iGraphs[i0].Fit("f1")
    lLeg.Draw()
    lC0.Modified()
    lC0.Update()
    lC0.SaveAs(iName+".png")
    lC0.SaveAs(iName+".pdf")
    end()

def profileRes(iTree,iName,iVar,iCut,iPt):
    lPt = array('d', [])
    for pPt in iPt:
        lPt.append(pPt)
    lProf = r.TProfile(lName,lName,len(iPt)-1,lPt)
    lProf.BuildOptions(-1,1,"s")
    iTree.Draw(lVar+">>"+lName,lCut)   
    lX = array('d', [])
    lY = array('d', [])
    lE = array('d', [])
    for i0 in range(1,lProf.GetNbinsX()+1):
        print lProf.GetXaxis().GetBinCenter(i0),lProf.GetBinContent(i0),lProf.GetBinError(i0)
        lX.append(lProf.GetXaxis().GetBinCenter(i0))
        lY.append(lProf.GetBinContent(i0))
        lE.append(lProf.GetBinError(i0))
        #lE.append(lProf.GetBinError(i0)/(1.+lProf.GetBinContent(i0)))
    return (lX,lY,lE)

def gausRes(iTree,iName,iVar,iCut,iPt):
    lX   = array('d', [])
    lY   = array('d', [])
    lE   = array('d', [])
    lEX  = array('d', [])
    lEY  = array('d', [])
    lEE  = array('d', [])
    for pPt in range(0,len(iPt)-1):
        pH = r.TH1F("ptmp","ptmp",50,-50.2,50.2)
        #pH = r.TH1F("ptmp","ptmp",20,-2.2,2.2)
        pCut = iCut+"*(genPt > "+str(iPt[pPt])+" && genPt < "+str(iPt[pPt+1])+")"
        iTree.Draw(iVar+">>ptmp",pCut)
        pH.Fit("gaus")
        lX.append((iPt[pPt]+iPt[pPt+1])/2.)
        lY.append(pH.GetFunction("gaus").GetParameter(1)/lX[len(lX)-1])
        lE.append(pH.GetFunction("gaus").GetParameter(2)/lX[len(lX)-1])
        lEX.append(abs(iPt[pPt]-iPt[pPt+1])/math.sqrt(12.))
        lEY.append(pH.GetFunction("gaus").GetParError(1)/lX[len(lX)-1])
        lEE.append(pH.GetFunction("gaus").GetParError(2)/lX[len(lX)-1])
    return (lX,lY,lE,lEX,lEY,lEE)

def resolution(iTree,iVar,iGVar,iEtaMin,iEtaMax,iPt):
    #lCut="("+str(iEtaMin)+" < genEta  && genEta < "+str(iEtaMax)+")*("+iVar+"> -0.1)*(trkPt > 1)*(abs(trkPt-"+iGVar+") < 10)"
    lCut="("+str(iEtaMin)+" < genEta  && genEta < "+str(iEtaMax)+")*("+iVar+"> -0.1)*(hcal_corr_et > 1)*(abs(hcal_corr_et-"+iGVar+") < 0.5*hcal_corr_et+10)"
    #lVar="("+iVar+"-"+iGVar+")/genPt"
    lVar="("+iVar+"-"+iGVar+")"
    print lVar
    lName=iVar+iGVar+str(iEtaMin)+str(iEtaMax)
    lName=lName.replace("-","")
    lName=lName.replace(".","")
    #lX,lY,lE = profileRes(iTree,lName,lVar+":genPt",lCut,iPt)
    lX,lY,lE,lEX,lEY,lEE = gausRes(iTree,lName,lVar,lCut,iPt)
    lMGraph = r.TGraphErrors(len(iPt)-1,lX,lY,lEX,lEY)
    lEGraph = r.TGraphErrors(len(iPt)-1,lX,lE,lEX,lEE)
    lMGraph.GetXaxis().SetTitle("p_{T}^{Gen} (GeV)")
    lEGraph.GetXaxis().SetTitle("p_{T}^{Gen} (GeV)")
    lMGraph.GetYaxis().SetTitle("(p_{T}^{TK}-p_{T}^{Gen})/p_{T}^{Gen}")
    lEGraph.GetYaxis().SetTitle("#sigma(p_{T}^{TK}-p_{T}^{Gen})/p_{T}^{Gen}")
    return (lMGraph,lEGraph)

if __name__ == "__main__":
    print args.Tree,args.input
    lFile = r.TFile(args.input)
    lTree = lFile.Get(args.Tree)
    lMu   = []
    lRes  = []
    for i0 in range(0,len(args.etarange)-1):
        pMu,pRes=resolution(lTree,args.var,args.genvar,args.etarange[i0],args.etarange[i0+1],args.ptrange)
        lMu.append(pMu)
        lRes.append(pRes)
    draw(lMu,args.etarange,"a")
    draw(lRes,args.etarange,"b",True)
    #exit()
    lOFile = r.TFile("ResOutput.root","RECREATE")
    lMu.extend(lRes)
    for pGraph in lMu:
        pGraph.Write()
    
