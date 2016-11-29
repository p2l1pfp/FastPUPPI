#! /usr/bin/env python
import ROOT as r,sys,math,array,os,argparse,numpy
from array import array
from tools import end
import tdrstyle
tdrstyle.setTDRStyle()

parser = argparse.ArgumentParser(description='Process benchmarks.')
parser.add_argument('-input'  ,'--input'     ,action='store',dest='input'  ,default='ResOutput.root'    ,help='input file')
parser.add_argument('-pt'     ,'--pt'        ,dest='ptrange' ,nargs='+',type=float,default=range(2,60,2),help='eta range')
parser.add_argument('-frac'   ,'--frac'      ,dest='fracrange',nargs='+',type=float,default=[0],help='em fraction range')
#parser.add_argument('-eta'    ,'--eta'       ,dest='etarange',nargs='+',type=float,default=[-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0],help='eta range')
parser.add_argument('-eta'    ,'--eta'       ,dest='etarange',nargs='+',type=float,default=[-2.5,-2.25,-2.0,-1.75,-1.5,-1.0,-0.5,0.,0.5,1.0,1.5,1.75,2.0,2.25,2.5],help='eta range')
#parser.add_argument('-frac'   ,'--frac'      ,dest='fracrange',nargs='+',type=float,default=[x * 0.1 for x in range(0, 11)],help='em fraction range')
#parser.add_argument('-eta'    ,'--eta'       ,dest='etarange',nargs='+',type=float,default=[-3.0,-2.0,-1.5,0.,1.5,2.0,3.0],help='eta range')

args = parser.parse_args()

fColor=[r.kRed,r.kBlue,r.kGreen+1,r.kViolet,r.kRed-2,r.kBlue-2,r.kSpring,r.kAzure+1,r.kOrange]

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

def loadCorrectorFits(iFile,iEta):
    lFile = r.TFile(iFile)
    lPar0=[]
    for i0 in range(0,len(iEta)-1):
        lG1=lFile.Get(str(iEta[i0])+"< #eta < "+str(iEta[i0+1])+"b")
        lPar0.append(lG1.GetFunction("f1"))
    return lPar0d

def translateIEta(eta):
    towerEtas = [0,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696,0.783,0.870,0.957,1.044,1.131,1.218,1.305,1.392,1.479,1.566,1.653,1.740,1.830,1.930,2.043,2.172,2.322,2.5,2.650,2.853,3.139,3.314,3.489,3.664,3.839,4.013,4.191,4.363,4.538,4.716,4.889,5.191]
    eta=eta-41
    outeta = towerEtas[abs(eta)-1]*0.5+towerEtas[abs(eta)]*0.5
    if eta < 0:
        outeta = outeta * -1
    return outeta

def buildFitTables(iPar0,iEtaMin,iEtaMax,iHPt,iHEta):
    lCorr = []
    for i0 in range(iEtaMin,iEtaMax):
        pEta0=translateIEta(i0)
        pEta=-1
        pEta2=-1
        pEtaFrac = 0
        for pEta1 in range(1,len(iHEta)-1):
            if pEta0 < iHEta[pEta1-1] or pEta0 > iHEta[pEta1]:
                continue
            pEta=pEta1-1
            pEta2=pEta1
            pEtaFrac = (pEta0-iHEta[pEta1])/(iHEta[pEta1+1]-iHEta[pEta1])
        if pEta == -1 and pEta0 < iHEta[0]:
            pEta  = 0
            pEta2 = 0
        if pEta == -1 and pEta0 > iHEta[len(iHEta)-1]:
            pEta  = len(iHEta)-2
            pEta2 = len(iHEta)-2
        print pEtaFrac,pEta0,pEta,i0
        lCorrFrac=[]
        lX   = array('d', [])
        lY   = array('d', [])
        for i2 in range(0,len(iHPt)-1):
            lA = (iPar0[pEta].Eval(iHPt[i2])*(1.-pEtaFrac)+iPar0[pEta2].Eval(iHPt[i2])*(pEtaFrac))
            pRes = iHPt[i2]*lA
            lX.append(iHPt[i2])
            lY.append(pRes)
        pGraph = r.TGraph(len(lX),lX,lY)
        pGraph.SetTitle("eta_res_"+str(i0))
        pGraph.SetName ("eta_res_"+str(i0))
        lCorrFrac.append(pGraph)
        lCorr.append(lCorrFrac)
    return lCorr
    
def writeCorr(iCorr):
    lFile = r.TFile("Output.root","RECREATE")
    for pCorrFrac in lCorr:
        for pGraph in pCorrFrac:
            pGraph.Write()

if __name__ == "__main__":
    print args.input
    lPar0,lPar1=loadCorrector(args.input,args.etarange)
    lCorr=buildFitTables(lPar1,1,82,args.ptrange,args.etarange)
    writeCorr(lCorr)
    
