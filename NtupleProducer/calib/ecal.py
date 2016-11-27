#! /usr/bin/env python
import ROOT as r,sys,math,array,os,argparse,numpy
from array import array
from tools import end
import tdrstyle
tdrstyle.setTDRStyle()

parser = argparse.ArgumentParser(description='Process benchmarks.')
parser.add_argument('-Tree'   ,'--Tree'      ,action='store',dest='Tree'   ,default='EcalInfo'      ,help='Tree Name')
parser.add_argument('-input'  ,'--input'     ,action='store',dest='input'  ,default='PhotonGun.root',help='input file')
parser.add_argument('-var'    ,'--var'       ,action='store',dest='var'    ,default='ecal_clust_et' ,help='var') 
parser.add_argument('-genvar' ,'--genvar'    ,action='store',dest='genvar' ,default='genPt'         ,help='mass') 
parser.add_argument('-build'  ,'--build'     ,action='store_true',dest='build'  ,default=False      ,help='build')
#parser.add_argument('-eta'    ,'--eta'       ,dest='etarange',nargs='+',type=float,default=[-2.5,-2.25,-2.0,-1.75,-1.5,-1.0,-0.5,0.,0.5,1.0,1.5,1.75,2.0,2.25,2.5],help='eta range')
parser.add_argument('-eta'    ,'--eta'       ,dest='etarange',nargs='+',type=float,default=[-3.0,-2.0,-1.5,0.,1.5,2.0,3.0],help='eta range')
parser.add_argument('-pt'     ,'--pt'        ,dest='ptrange' ,nargs='+',type=float,default=range(0,50,4),help='pt range')
parser.add_argument('-frac'   ,'--frac'      ,dest='fracrange',nargs='+',type=float,default=[0],help='em fraction range')

args = parser.parse_args()

fColor=[r.kRed,r.kBlue,r.kGreen+1,r.kViolet,r.kRed-2,r.kBlue-2,r.kSpring,r.kAzure+1,r.kOrange,r.kRed,r.kBlue,r.kGreen+1,r.kViolet,r.kRed-2,r.kBlue-2,r.kSpring,r.kAzure+1,r.kOrange]

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
        print i0,iGraphs,fColor[i0]
        iGraphs[i0].SetLineColor(fColor[i0])
        iGraphs[i0].SetMarkerColor(fColor[i0])
        iGraphs[i0].Draw("lpe")
        lLeg.AddEntry(iGraphs[i0],str(iEta[i0])+"< #eta < "+str(iEta[i0+1]),"lpe")
        iGraphs[i0].SetName(str(iEta[i0])+"< #eta < "+str(iEta[i0+1])+iName)
        iGraphs[i0].SetTitle(str(iEta[i0])+"< #eta < "+str(iEta[i0+1])+iName)
        pF1 = r.TF1("f1","[0]+[1]/sqrt(x)+[2]/x",0.1,200);
        #if iGaus:
        #    iGraphs[i0].Fit("f1")
    lLeg.Draw()
    lC0.Modified()
    lC0.Update()
    #lC0.SaveAs(iName+".png")
    #lC0.SaveAs(iName+".pdf")
    end()

def profileRes(iTree,iName,iVar,iCut,iPt):
    lPt = array('d', [])
    for pPt in iPt:
        lPt.append(pPt)
    lProf = r.TProfile(iName,iName,len(iPt)-1,lPt)
    lProf.BuildOptions(-1,1,"s")
    iTree.Draw(iVar+">>"+iName,iCut)   
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
        pH = r.TH1F("ptmp","ptmp",100,-20.2,20.2)
        pCut = iCut+"*(genPt > "+str(iPt[pPt])+" && genPt < "+str(iPt[pPt+1])+")"
        iTree.Draw(iVar+">>ptmp",pCut)
        pH.Fit("gaus","L")
        lX.append((iPt[pPt]+iPt[pPt+1])/2.)
        lY.append(pH.GetFunction("gaus").GetParameter(1)/lX[len(lX)-1])
        lE.append(pH.GetFunction("gaus").GetParameter(2)/lX[len(lX)-1])
        lEX.append(abs(iPt[pPt]-iPt[pPt+1])/math.sqrt(12.))
        lEY.append(pH.GetFunction("gaus").GetParError(1)/lX[len(lX)-1])
        lEE.append(pH.GetFunction("gaus").GetParError(2)/lX[len(lX)-1])
    return (lX,lY,lE,lEX,lEY,lEE)

def resolution(iTree,iVar,iGVar,iEtaMin,iEtaMax,iPt):
    #lCut="("+str(iEtaMin)+" < genEta  && genEta < "+str(iEtaMax)+")*("+iVar+"> -0.1)*(hcal_corr_et > 1)"#-"+iGVar+") < 50)"
    lCut="("+str(iEtaMin)+" < genEta  && genEta < "+str(iEtaMax)+")*("+iVar+" > 1)*(("+iVar+"/"+iGVar+") > 0.8)"
    #lVar="("+iVar+"-"+iGVar+")/genPt"
    lVar="("+iVar+"-"+iGVar+")"
    lName=iVar+iGVar+str(iEtaMin)+str(iEtaMax)
    lName=lName.replace("-","")
    lName=lName.replace(".","")
    lX,lY,lE,lEX,lEY,lEE = gausRes(iTree,lName,lVar,lCut,iPt)
    #lX,lY,lE = profileRes(iTree,lName,lVar+":genPt",lCut,iPt)
    lMGraph = r.TGraphErrors(len(iPt)-1,lX,lY,lEX,lEY)
    lEGraph = r.TGraphErrors(len(iPt)-1,lX,lE,lEX,lEE)
    lMGraph.GetXaxis().SetTitle("p_{T}^{Gen} (GeV)")
    lEGraph.GetXaxis().SetTitle("p_{T}^{Gen} (GeV)")
    lMGraph.GetYaxis().SetTitle("(p_{T}^{Clust}-p_{T}^{Gen})/p_{T}^{Gen}")
    lEGraph.GetYaxis().SetTitle("#sigma(p_{T}^{Clust}-p_{T}^{Gen})/p_{T}^{Gen}")
    return (lMGraph,lEGraph)

def loadCorrector(iFile,iEta):
    lFile = r.TFile(iFile)
    lPar0=[]
    lPar1=[]
    for i0 in range(0,len(iEta)-1):
        lG0=lFile.Get(str(iEta[i0])+"< #eta < "+str(iEta[i0+1])+"a")
        lG1=lFile.Get(str(iEta[i0])+"< #eta < "+str(iEta[i0+1])+"b")
        lPar0.append(lG0)
        lPar1.append(lG1)
    return lPar0,lPar1

def clean(iX,iY):
    lX   = array('d', [])
    lY   = array('d', [])
    for i0 in range(0,len(iX)-1):
        if iX[i0] > iX[i0+1]-3:
                continue
        lX.append(iX[i0])
        lY.append(iY[i0])
    return lX,lY

def translateIEta(eta):
    towerEtas = [0,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696,0.783,0.870,0.957,1.044,1.131,1.218,1.305,1.392,1.479,1.566,1.653,1.740,1.830,1.930,2.043,2.172,2.322,2.5,2.650,2.853,3.139,3.314,3.489,3.664,3.839,4.013,4.191,4.363,4.538,4.716,4.889,5.191]
    eta=eta-41
    outeta = towerEtas[abs(eta)-1]*0.5+towerEtas[abs(eta)]*0.5
    if eta < 0:
        outeta = outeta * -1
    return outeta

def buildTables(iFrac,iPar0,iPar1,iHPt,iHEta):
    lCorr = []
    for i0 in range(1,82):
        pEta0=translateIEta(i0)
        pEta = -1
        for pEta1 in range(0,len(iHEta)-1):
            if pEta0 < iHEta[pEta1] or pEta0 > iHEta[pEta1+1]:
                continue
            pEta=pEta1
        if pEta == -1 and pEta0 < iHEta[0]:
            pEta = 0
        if pEta == -1 and pEta0 > iHEta[len(iHEta)-1]:
            pEta = len(iHEta)-2
        print pEta,pEta0,i0,len(iHEta)-1,iHEta[len(iHEta)-1]
        lCorrFrac=[]
        for i1 in range(0,len(iFrac)):
            lX   = array('d', [])
            lY   = array('d', [])
            for i2 in range(0,len(iHPt)-1):
                lA = (iPar0[pEta].Eval(iHPt[i2])+iPar0[pEta].Eval(iHPt[i2]))/2.
                lB = (iPar1[pEta].Eval(iHPt[i2])+iPar1[pEta].Eval(iHPt[i2]))/2.
                pPt = iHPt[i2]*lA + iHPt[i2] + lB*iFrac[i1]
                lY.append(iHPt[i2])
                lX.append(pPt)
            lX,lY=clean(lX,lY)
            pGraph = r.TGraph(len(lX),lX,lY)
            pGraph.SetTitle("eta_"+str(i0)+"_frac_"+str(i1))
            pGraph.SetName ("eta_"+str(i0)+"_frac_"+str(i1))
            lCorrFrac.append(pGraph)
        lCorr.append(lCorrFrac)
    return lCorr

def writeCorr(iCorr):
    lFile = r.TFile("Output.root","RECREATE")
    for pCorrFrac in lCorr:
        for pGraph in pCorrFrac:
            pGraph.Write()

if __name__ == "__main__":
    print args.Tree,args.input
    if args.build:
        lPar0,lPar1=loadCorrector(args.input,args.etarange)
        lCorr=buildTables(args.fracrange,lPar0,lPar1,args.ptrange,args.etarange)
        writeCorr(lCorr)
        exit()
    lFile = r.TFile(args.input)
    lTree = lFile.Get(args.Tree)
    lMu   = []
    lRes  = []
    for i0 in range(0,len(args.etarange)-1):
        pMu,pRes=resolution(lTree,args.var,args.genvar,args.etarange[i0],args.etarange[i0+1],args.ptrange)
        lMu.append(pMu)
        lRes.append(pRes)
        #lProf.append(pProf)
    draw(lMu,args.etarange,"a")
    draw(lRes,args.etarange,"b",True)
    lOFile = r.TFile("Output.root","RECREATE")
    lMu.extend(lRes)
    for pGraph in lMu:
        pGraph.Write()
    
