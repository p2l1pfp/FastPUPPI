#! /usr/bin/env python
import ROOT as r,sys,math,array,os,argparse,numpy
from array import array
from tools import end
import tdrstyle
tdrstyle.setTDRStyle()

parser = argparse.ArgumentParser(description='Process benchmarks.')
parser.add_argument('-Tree'   ,'--Tree'      ,action='store',dest='Tree',default='HcalInfo',help='Tree Name')
parser.add_argument('-input'  ,'--input'     ,action='store',dest='input'  ,default='Pion2.root',help='input file')
parser.add_argument('-var'    ,'--var'       ,action='store',dest='var'    ,default='hcal_corr_et'         ,help='mass') 
parser.add_argument('-genvar' ,'--genvar'    ,action='store',dest='genvar' ,default='genPt'         ,help='mass') 
parser.add_argument('-eta'    ,'--eta'       ,dest='etarange',nargs='+',type=float,default=[-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0],help='eta range')
parser.add_argument('-pt'     ,'--pt'        ,dest='ptrange' ,nargs='+',type=float,default=range(0,60,3),help='eta range')
parser.add_argument('-frac'   ,'--frac'      ,dest='fracrange',nargs='+',type=float,default=[0],help='em fraction range')
args = parser.parse_args()

fColor=[r.kRed,r.kBlue,r.kGreen+1,r.kViolet,r.kRed-2,r.kBlue-2,r.kSpring,r.kAzure+1,r.kOrange]

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

def binnedSurfaceFitCmd(iTree,iPtMin,iPtMax,iEtaMin,iEtaMax,iPhiMin,iPhiMax):
    lCut = "(genPt > "+str(iPtMin)+" && genPt < "+str(iPtMax)+" && hcal_ieta == "+str(iEtaMin)+" && hcal_corr_et > 0)" 
    lVar = "hcal_et/(hcal_et+ecal_et):hcal_corr_et-genPt"
    print lVar,lCut
    print lCut
    lN   = iTree.Draw(lVar,lCut,"goff")
    if lN == 0:
        return (0,0,0,0)
    lGraph = r.TGraph(lN,iTree.GetV1(),iTree.GetV2())
    lPass=False
    for i0 in range(0,lN-1):
        if iTree.GetV2()[i0] != iTree.GetV2()[i0+1]:
            lPass=True
    l1D = False if float(iTree.GetEntries("(ecal_et > 0)*"+lCut))/float(iTree.GetEntries(lCut)) > 0.05 else True
    label = "pol0" if l1D else "pol1"
    lGraph.Fit (label)
    lMean   = lGraph.GetFunction(label).GetParameter(0)
    lSlope  = lGraph.GetFunction(label).GetParameter(1) if not l1D else 0
    lEMean  = lGraph.GetFunction(label).GetParError(0)
    lESlope = lGraph.GetFunction(label).GetParError(1)  if not l1D else 0
    if not lPass:
        lMean=0
        lEMean=0
        lSlope=0
        lESlope=0.1
        for i0 in range(0,lN):
            lMean = iTree.GetV1()[i0] + lMean
        lMean=lMean/lN
        lEMean=lMean/math.sqrt(lN)
    return (lMean,lSlope,lEMean,lESlope)
    
def binnedSurfaceFit(iTree,iPtMin,iPtMax,iEtaMin,iEtaMax,iPhiMin,iPhiMax):
    lEFrac = array('d', [])
    lResp  = array('d', [])
    n=0
    lOld=-999
    lPass=True
    for i0 in range(0,iTree.GetEntriesFast()):
        iTree.GetEntry(i0)
        if i0 % 1000000 == 0:
            print i0
        pPt   = iTree.genPt
        pEta  = iTree.hcal_ieta #genEta
        pPhi  = iTree.hcal_iphi #genPhi
        if pPt < iPtMin or pPt > iPtMax or iTree.hcal_clust_et == 0:
            continue
        #if pEta < iEtaMin or pEta > iEtaMax:
        #    continue
        if pEta != iEtaMin:
            continue
        pEcal = iTree.ecal_et
        pHcal = iTree.hcal_corr_et
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

def surfaceFit(iTree,iPt,iEtaMin,iEtaMax,iPhiMin,iPhiMax):
    lX   = array('d', [])
    lEX  = array('d', [])
    lY0  = array('d', [])
    lEY0 = array('d', [])
    lY1  = array('d', [])
    lEY1 = array('d', []) 
    for pPt in range(0,len(iPt)-1):
        pY0,pY1,pEY0,pEY1 = binnedSurfaceFitCmd(iTree,iPt[pPt],iPt[pPt+1],iEtaMin,iEtaMax,iPhiMin,iPhiMax)
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

if __name__ == "__main__":
    print args.Tree,args.input
    lFile = r.TFile.Open(args.input)
    lTree = lFile.Get(args.Tree)
    #if True:
        #surfaceFit(lTree,args.ptrange,-1.0,1.0)
        #resp2D(lTree)
        #exit()
    lMu   = []
    lRes  = []
    for iEta in range(-40,41):
        print "Eta",iEta
        pMu,pRes=surfaceFit(lTree,args.ptrange,iEta,iEta,-1000,1000)
        lMu.append(pMu)
        lRes.append(pRes)
    lOFile = r.TFile("Output.root","RECREATE")
    lMu.extend(lRes)
    for pGraph in lMu:
        pGraph.Write()
    
