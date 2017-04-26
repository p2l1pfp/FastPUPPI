#! /usr/bin/env python
import ROOT as r,sys,math,array,os,argparse,numpy
from array import array
from tools import end
import tdrstyle
tdrstyle.setTDRStyle()

parser = argparse.ArgumentParser(description='Process benchmarks.')
parser.add_argument('-Tree'   ,'--Tree'      ,action='store',dest='Tree',default='HcalInfo',help='Tree Name')
parser.add_argument('-input'  ,'--input'     ,action='store',dest='input'  ,default='root://eoscms//eos/cms/store/cmst3/user/pharris/FastPUPPI/ntuple_pion_81X_v3.root',help='input file')
parser.add_argument('-var'    ,'--var'       ,action='store',dest='var'    ,default='min(ecal_et,0)+hcal_clust_et'         ,help='mass') 
parser.add_argument('-genvar' ,'--genvar'    ,action='store',dest='genvar' ,default='genPt'         ,help='mass') 
parser.add_argument('-pt'     ,'--pt'        ,dest='ptrange' ,nargs='+',type=float,default=range(0,70,5),help='eta range')
parser.add_argument('-frac'   ,'--frac'      ,dest='fracrange',nargs='+',type=float,default=[x * 0.1 for x in range(0, 11)],help='em fraction range')
parser.add_argument('-ieta'   ,'--ieta'      ,action='store',dest='ieta',default=-30,help='ieta')
parser.add_argument('-iphi'   ,'--iphi'      ,action='store',dest='iphi',default=1,help='iphi')
args = parser.parse_args()

fColor=[r.kRed,r.kBlue,r.kGreen+1,r.kViolet,r.kRed-2,r.kBlue-2,r.kSpring,r.kAzure+1,r.kOrange]

def draw(iGraphs,iEta,iName="A"):
    lC0 = r.TCanvas("Can"+iName,"Can"+iName,800,600);
    lLeg = r.TLegend(0.2,0.65,0.45,0.85)
    lLeg.SetBorderSize(0)
    lLeg.SetFillColor(0)
    #iGraphs[0].GetYaxis().SetRangeUser(-0.1,0.1)
    iGraphs[0].Draw("alep")
    lLeg.AddEntry(iGraphs[0],str(iEta[0])+"< #eta < "+str(iEta[0+1]),"lpe")
    iGraphs[0].SetName (str(iEta[0])+"< #eta < "+str(iEta[0+1])+iName)
    iGraphs[0].SetTitle(str(iEta[0])+"< #eta < "+str(iEta[0+1])+iName)
    for i0 in range(1,len(iGraphs)):
        print i0,iGraphs,fColor[i0]
        iGraphs[i0].SetLineColor(fColor[i0])
        iGraphs[i0].SetMarkerColor(fColor[i0])
        iGraphs[i0].Draw("lpe")
        lLeg.AddEntry(iGraphs[i0],str(iEta[i0])+"< #eta < "+str(iEta[i0+1]),"lpe")
        iGraphs[i0].SetName(str(iEta[i0])+"< #eta < "+str(iEta[i0+1])+iName)
        iGraphs[i0].SetTitle(str(iEta[i0])+"< #eta < "+str(iEta[i0+1])+iName)
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
        pH = r.TH1F("ptmp","ptmp",100,-100.2,200.2)
        pCut = iCut+"*(genPt > "+str(iPt[pPt])+" && genPt < "+str(iPt[pPt+1])+")"
        iTree.Draw(iVar+">>ptmp",pCut)
        pH.Fit("gaus")
        lX.append((iPt[pPt]+iPt[pPt+1])/2.)
        lY.append(pH.GetFunction("gaus").GetParameter(1))
        lE.append(pH.GetFunction("gaus").GetParameter(2))
        lEX.append(abs(iPt[pPt]-iPt[pPt+1])/math.sqrt(12.))
        lEY.append(pH.GetFunction("gaus").GetParError(1))
        lEE.append(pH.GetFunction("gaus").GetParError(2))
    return (lX,lY,lE,lEX,lEY,lEE)

def resolution(iTree,iVar,iGVar,iEtaMin,iEtaMax,iPt):
    lCut="("+str(iEtaMin)+" < genEta  && genEta < "+str(iEtaMax)+")*("+iVar+"> -0.1)*(hcal_clust_et > 1)"#-"+iGVar+") < 50)"
    #lVar="("+iVar+"-"+iGVar+")/genPt"
    lVar="("+iVar+"-"+iGVar+")"
    lName=iVar+iGVar+str(iEtaMin)+str(iEtaMax)
    lName=lName.replace("-","")
    lName=lName.replace(".","")
    #lX,lY,lE = profileRes(iTree,lName,lVar+":genPt",lCut,iPt)
    lX,lY,lE,lEX,lEY,lEE = gausRes(iTree,lName,lVar,lCut,iPt)
    lMGraph = r.TGraphErrors(len(iPt)-1,lX,lY,lEX,lEY)
    lEGraph = r.TGraphErrors(len(iPt)-1,lX,lE,lEX,lEE)
    return (lMGraph,lEGraph)

def resp2D(iTree,iPt=20):
    lH = r.TH2F("TD","TD",40,-5,5,40,-3.15,3.15)
    lN = r.TH2F("TN","TN",40,-5,5,40,-3.15,3.15)
    lLabel="genPhi:genEta"
    lWeight="(hcal_et*2.+ecal_et)/genPt"
    lCut   ="(genPt > "+str(iPt)+" &&  hcal_et+ecal_et > 5)"
    iTree.Draw(lLabel+">>"+lH.GetName(),lCut+"*"+lWeight,"colz")
    iTree.Draw(lLabel+">>"+lN.GetName(),lCut            ,"colz")
    lH.Divide(lN)
    lH.GetXaxis().SetTitle("#eta")
    lH.GetYaxis().SetTitle("#phi")
    lH.Draw("colz")
    end()

#def surfaceFit(iTree):
#    lData = r.RooDataSet()
#    for i0 in range(0,iTree.GetEntriesFast()):
#        pGen  = iTree.genPt
#        pEcal = iTree.ecal_et
#        pHcal = iTree.hcal_clust_et
#        lData.add(r.RooArgSet(pGen,pEcal,pHcal))
#    lFit  = r.RooGenericPdf("A","A"

def binnedSurfaceFit(iTree,iPtMin,iPtMax,iEtaMin,iEtaMax,iPhiMin,iPhiMax):
    lEFrac = array('d', [])
    lResp  = array('d', [])
    n=0
    lOld=-999
    lPass=True
    for i0 in range(0,iTree.GetEntriesFast()):
        iTree.GetEntry(i0)
        if i0 % 1000000 == 0:
            print i0,"---",iPtMin,"-",iPtMax,"--",iEta,iPhi
        pPt   = iTree.genPt
        pEta  = iTree.hcal_ieta #genEta
        pPhi  = iTree.hcal_iphi #genPhi
        if pPt < iPtMin or pPt > iPtMax or iTree.hcal_clust_et ==0:
            continue
        #if pEta < iEtaMin or pEta > iEtaMax:
        #    continue
        pDPhi = abs(pPhi-iPhiMin)
        if 72-pDPhi < pDPhi :
            pDPhi = 72 - pDPhi
        if pEta != iEtaMin or pDPhi > 20:
            continue
        pEcal = iTree.ecal_et
        pHcal = iTree.hcal_clust_et
        pEcal = max(pEcal,0)
        lEFrac.append(pEcal/(pHcal+pEcal))
        lResp .append((pEcal+pHcal)-pPt)
        #print pEcal,pHcal,pPt,lResp[len(lRes)-1],lEFrac[len(lEFrac)-1]
        if lOld != -999 and lOld != lEFrac[len(lEFrac)-1]:
            lPass = False
        lOld = lEFrac[len(lEFrac)-1] 
        n=n+1
    if n == 0:
        return (0,0,0,0)
    lGraph = r.TGraph(n,lEFrac,lResp)
    lGraph.Fit("pol1")
    lGraph.Draw("ap")
    lMean   = lGraph.GetFunction("pol1").GetParameter(0)
    lSlope  = lGraph.GetFunction("pol1").GetParameter(1)
    lEMean  = lGraph.GetFunction("pol1").GetParError(0)
    lESlope = lGraph.GetFunction("pol1").GetParError(1)
    if lPass:
        lMean=0
        lEMean=0
        lSlope=0
        lESlope=0.1
        for pVar in lResp:
            lMean = pVar + lMean
        lMean=lMean/len(lResp)
        lEMean=lMean/math.sqrt(len(lResp))
    return (lMean,lSlope,lEMean,lESlope)

def ptFit(iTree,iEtaMin,iEtaMax,iPhiMin,iPhiMax):
    lPt    = array('d', [])
    lResp  = array('d', [])
    n=0
    lOld=-999
    lPass=True
    for i0 in range(0,iTree.GetEntriesFast()):
        iTree.GetEntry(i0)
        if i0 % 1000000 == 0:
            print i0,"---",iEta,iPhi
        pPt   = iTree.genPt
        pEta  = iTree.hcal_ieta #genEta
        pPhi  = iTree.hcal_iphi #genPhi
        if iTree.hcal_clust_et == 0:
            continue
        pDPhi = abs(pPhi-iPhiMin)
        if 72-pDPhi < pDPhi :
            pDPhi = 72 - pDPhi
        if pEta != iEtaMin or pDPhi > 20:
            continue
        pEcal = iTree.ecal_et
        pHcal = iTree.hcal_clust_et
        pEcal = max(pEcal,0)
        if pEcal != 0:
            continue
        lResp .append((pEcal+pHcal)-pPt)
        lPt   .append(pPt)
        n=n+1
    if n == 0:
        return (0,0,0,0)
    lGraph = r.TGraph(n,lPt,lResp)
    lGraph.Fit("pol1")
    lMean   = lGraph.GetFunction("pol1").GetParameter(0)
    lSlope  = lGraph.GetFunction("pol1").GetParameter(1)
    lEMean  = lGraph.GetFunction("pol1").GetParError(0)
    lESlope = lGraph.GetFunction("pol1").GetParError(1)
    return (lMean,lSlope,lEMean,lESlope)

def surfacePtFit(iTree,iPt,iEtaMin,iEtaMax,iPhiMin,iPhiMax):
    lX   = array('d', [])
    lEX  = array('d', [])
    lY0  = array('d', [])
    lEY0 = array('d', [])
    lMean,lSlope,lEMean,lESlope = ptFit(iTree,iEtaMin,iEtaMax,iPhiMin,iPhiMax)
    for pPt in range(0,len(iPt)-1):
        pX = (iPt[pPt]+iPt[pPt+1])/2.
        lX.append(pX)
        lEX.append(abs(iPt[pPt]-iPt[pPt+1])/math.sqrt(12.))
        lY0 .append(lMean + lSlope * pX)
        lEY0.append(math.sqrt(lEMean*lEMean + lESlope*lESlope*pX*pX))
    lMGraph = r.TGraphErrors(len(iPt)-1,lX,lY0,lEX,lEY0)
    lMGraph.SetTitle("Mean0_" +str(iEtaMin)+"_"+str(iPhiMin))
    lMGraph.SetName ("Mean0_" +str(iEtaMin)+"_"+str(iPhiMin))
    return lMGraph

def surfaceFit(iTree,iPt,iEtaMin,iEtaMax,iPhiMin,iPhiMax):
    lX   = array('d', [])
    lEX  = array('d', [])
    lY0  = array('d', [])
    lEY0 = array('d', [])
    lY1  = array('d', [])
    lEY1 = array('d', []) 
    for pPt in range(0,len(iPt)-1):
        pY0,pY1,pEY0,pEY1 = binnedSurfaceFit(iTree,iPt[pPt],iPt[pPt+1],iEtaMin,iEtaMax,iPhiMin,iPhiMax)
        lX.append((iPt[pPt]+iPt[pPt+1])/2.)
        lEX.append(abs(iPt[pPt]-iPt[pPt+1])/math.sqrt(12.))
        lY0 .append(pY0)
        lEY0.append(pEY0)
        lY1 .append(pY1)
        lEY1.append(pEY1)
    lMGraph = r.TGraphErrors(len(iPt)-1,lX,lY0,lEX,lEY0)
    lEGraph = r.TGraphErrors(len(iPt)-1,lX,lY1,lEX,lEY1)
    lMGraph.SetTitle("Mean_" +str(iEtaMin)+"_"+str(iPhiMin))
    lEGraph.SetTitle("Slope_"+str(iEtaMin)+"_"+str(iPhiMin))
    lMGraph.SetName("Mean_" +str(iEtaMin)+"_"+str(iPhiMin))
    lEGraph.SetName("Slope_"+str(iEtaMin)+"_"+str(iPhiMin))
    return lMGraph,lEGraph

def clean(iX,iY):
    lX   = array('d', [])
    lY   = array('d', [])
    for i0 in range(0,len(iX)-1):
        if iX[i0] > iX[i0+1]-3:
                continue
        lX.append(iX[i0])
        lY.append(iY[i0])
    return lX,lY

def buildTables(iFrac,iPar0,iPar1,iHPt,iHEta):
    lCorr = []
    lHEta = [-2.25,-2.0,-1.25,-0.75,-0.0,0.75,1.25,2.0,2.25]
    for i0 in range(0,61):
        pEta0=(i0-30.)/10.
        pEta = -1
        pEtaFrac = 0.
        for pEta1 in range(0,len(lHEta)-1):
            if pEta0 < lHEta[pEta1] or pEta0 > lHEta[pEta1+1]:
                continue
            pEta=pEta1
            pEtaFrac = (pEta0-lHEta[pEta1])/(lHEta[pEta1+1]-lHEta[pEta1])
        if pEta == -1 and pEta0 < lHEta[0]:
            pEta = 0
            pEtaFrac = 0.
        if (pEta == -1 and pEta0 > lHEta[len(lHEta)-1]) or pEta == 8:
            pEta     = 7
            pEtaFrac = 1.
        lCorrFrac=[]
        print "Test",pEta0,pEtaFrac,pEta,len(iHEta)
        for i1 in range(0,len(iFrac)):
            lX   = array('d', [])
            lY   = array('d', [])
            for i2 in range(0,len(iHPt)-1):
                lA = (iPar0[pEta].Eval(iHPt[i2])*(1.-pEtaFrac)+iPar0[pEta+1].Eval(iHPt[i2])*(pEtaFrac))
                lB = (iPar1[pEta].Eval(iHPt[i2])*(1.-pEtaFrac)+iPar1[pEta+1].Eval(iHPt[i2])*(pEtaFrac))
                #pPt = iHPt[i2] + lA + lB*iFrac[i1]
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
    lFile = r.TFile.Open(args.input)
    lTree = lFile.Get(args.Tree)
    lMu   = []
    lMu2  = []
    lRes  = []
    for iEta in range(-41,41+1):
        for iPhi in range(1,72+1):
            if iEta != int(args.ieta) or iPhi != int(args.iphi):
                continue
            print "fitting",iEta,iPhi
            pMu,pRes=surfaceFit(lTree,args.ptrange,iEta,iEta,iPhi,iPhi)
            pMu2    =surfacePtFit(lTree,args.ptrange,iEta,iEta,iPhi,iPhi)
            lMu2.append(pMu2)
            lMu.append(pMu)
            lRes.append(pRes)
    lOFile = r.TFile("Output.root","RECREATE")
    lMu.extend(lRes)
    lMu.extend(lMu2)
    for pGraph in lMu:
        pGraph.Write()
    
