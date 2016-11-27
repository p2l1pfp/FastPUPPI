#! /usr/bin/env python
import ROOT as r,sys,math,array,os,argparse,numpy
from array import array
from tools import end
import tdrstyle
tdrstyle.setTDRStyle()

parser = argparse.ArgumentParser(description='Process benchmarks.')
parser.add_argument('-input'  ,'--input'     ,action='store'  ,dest='input'  ,default='calib.root',help='input file')
parser.add_argument('-var'    ,'--var'       ,action='store'  ,dest='var'    ,default='max(ecal_et,0)+hcal_clust_et'         ,help='mass') 
parser.add_argument('-genvar' ,'--genvar'    ,action='store'  ,dest='genvar' ,default='genPt'         ,help='mass') 
parser.add_argument('-pt'     ,'--pt'        ,dest='ptrange'  ,nargs='+',type=float,default=range(0,100,2),help='eta range')
parser.add_argument('-frac'   ,'--frac'      ,dest='fracrange',nargs='+',type=float,default=[x * 0.1 for x in range(0, 11)],help='em fraction range')
args = parser.parse_args()

fColor=[r.kRed,r.kBlue,r.kGreen+1,r.kViolet,r.kRed-2,r.kBlue-2,r.kSpring,r.kAzure+1,r.kOrange]
fPar0=[]
fPar1=[]
fPar2=[]

def cleanG(iG0,iG1,iEta,iPhi):
    lX    = array('d', [])
    lY0   = array('d', [])
    lY1   = array('d', [])
    lEX   = array('d', [])
    lEY0  = array('d', [])
    lEY1  = array('d', [])
    notZero  = False
    for  i0 in range(0,iG1.GetN()):
        if iG1.GetY()[i0] != 0:
            notZero = True
    for i0 in range(0,iG0.GetN()):
        if iG1.GetY()[i0] == 0 and notZero: 
            continue
        if i0 < iG0.GetN()-1 and abs(iG0.GetErrorY(i0)) > 2.5*abs(iG0.GetErrorY(i0+1)):
            continue
        if i0 < iG1.GetN()-1 and abs(iG1.GetErrorY(i0)) > 2.5*abs(iG1.GetErrorY(i0+1)):
            continue
        if i0 > 0 and abs(iG0.GetErrorY(i0)) > 2.5*abs(iG0.GetErrorY(i0-1)):
            continue
        if i0 > 0 and abs(iG1.GetErrorY(i0)) > 2.5*abs(iG1.GetErrorY(i0-1)):
            continue
        lX  .append(iG0.GetX()[i0])
        lEX .append(iG0.GetErrorX(i0))
        lY0 .append(iG0.GetY()[i0])
        lEY0.append(iG0.GetErrorY(i0))
        lY1 .append(iG1.GetY()[i0])
        lEY1.append(iG1.GetErrorY(i0))
    lMGraph = r.TGraphErrors(len(lX),lX,lY0,lEX,lEY0)
    lEGraph = r.TGraphErrors(len(lX),lX,lY1,lEX,lEY1)
    lMGraph.SetTitle("Mean_" +str(iEta)+"_"+str(iPhi))
    lEGraph.SetTitle("Slope_"+str(iEta)+"_"+str(iPhi))
    lMGraph.SetName("Mean_" +str(iEta)+"_"+str(iPhi))
    lEGraph.SetName("Slope_"+str(iEta)+"_"+str(iPhi))
    return lMGraph,lEGraph

def loadCorrector(iFile,iEtaMin,iEtaMax,iPhiMin,iPhiMax):
    lFile = r.TFile(iFile)
    for iEta in range(iEtaMin,iEtaMax+1):
        for iPhi in range(iPhiMin,iPhiMax+1):
            print "Mean_" +str(iEta)+"_"+str(iPhi)
            lG0=lFile.Get("Mean_" +str(iEta)+"_"+str(iPhi))
            lG1=lFile.Get("Slope_"+str(iEta)+"_"+str(iPhi))
            #lG2=lFile.Get("Mean0_"+str(iEta)+"_"+str(iPhi))
            lG0,lG1=cleanG(lG0,lG1,iEta,iPhi)
            fPar0.append(lG0)
            fPar1.append(lG1)
            #fPar2.append(lG2)
    return fPar0,fPar1,fPar2

def writeX(iG):
    f = open('outgraph', 'a')
    f.write(iG.GetName()+' ')
    for i0 in range(0,iG.GetN()):
        f.write(str(iG.GetX()[i0])+' ')
        f.write(str(iG.GetY()[i0])+' ')
    f.write('\n')

def loadFile(iFile='outgraph'):
    f = open('outgraph', 'r')
    lPar0=[]
    lPar1=[]
    lX   = array('d', [])
    lY   = array('d', [])
    iC   = True
    for line in f:
        iC
        pName=''
        i0=-1
        for pVar in line.split(' '):
            if len(pVar) == 1:
                continue
            i0 = i0 + 1
            if i0 == 0:
                pName=pVar
                continue
            if i0 % 2 == 1:
                lX.append(float(pVar))
            if i0 % 2 == 1:
                lY.append(float(pVar))
        if pName.find("Mean") < 0:
            iC = False
        pGraph = r.TGraph(len(lX),lX,lY)
        pGraph.SetName(pName)    
        pGraph.SetTitle(pName)    
        lPar0.append(pGraph) if iC else lPar1.append(pGraph)
    return lPar0,lPar1

def write(iPar1,iPar2):
    for pG in iPar1:
        writeX(pG)
    for pG in iPar2:
        writeX(pG)

def clean(iX,iY):
    lX   = array('d', [])
    lY   = array('d', [])
    lSlope = 0
    lN=0
    for i0 in range(0,len(iX)):
        if iX[i0] == 0:
            continue
        lSlope = lSlope+(iY[i0]/iX[i0])
        lN=lN+1
    lSlope = lSlope/lN
    for i0 in range(0,len(iX)-1):
        if iX[i0] > iX[i0+1]:
            continue
        if i0 > 0 and (iY[i0]/iX[i0] > 2.3*lSlope or iY[i0]/iX[i0] < 0.4*lSlope):
            continue
        lX.append(iX[i0])
        lY.append(iY[i0])
    lX.append(iX[len(iX)-1])
    lY.append(iY[len(iX)-1])
    return lX,lY

def buildTables(iFrac,iPar0,iPar1,iPar2,iEtaMin,iEtaMax,iPhiMin,iPhiMax,iHPt):
    lCorr = []
    pId = -1
    for iEta     in range(iEtaMin,iEtaMax+1):
        if abs(iEta) != 41:
            pId=pId+1 
        for iPhi in range(iPhiMin,iPhiMax+1):
            pEta = -1
            lCorrFrac=[]
            print "Eta",iEta,"Phi",iPhi
            for i1 in range(0,len(iFrac)):
                lX   = array('d', [])
                lY   = array('d', [])
                for i2 in range(0,len(iHPt)-1):
                    lA  =fPar0[pId].Eval(iHPt[i2])
                    lB = fPar1[pId].Eval(iHPt[i2])
                    #if i1 == 0:
                    #    lB=0
                    #    lA=iPar2[pId].Eval(iHPt[i2])
                    pPt = iHPt[i2] + lA + lB*iFrac[i1]
                    #pPt = iHPt[i2]*lA + iHPt[i2] + lB*iFrac[i1]
                    lY.append(iHPt[i2])
                    lX.append(pPt)
                lX,lY=clean(lX,lY)
                pGraph = r.TGraph(len(lX),lX,lY)
                pGraph.SetTitle("eta_"+str(iEta)+"phi_"+str(iPhi)+"_frac_"+str(i1))
                pGraph.SetName ("eta_"+str(iEta)+"phi_"+str(iPhi)+"_frac_"+str(i1))
                lCorrFrac.append(pGraph)
            lCorr.append(lCorrFrac)
    return lCorr

def writeCorr(iCorr):
    lFile = r.TFile("Output.root","RECREATE")
    id0=0
    for pCorrFrac in iCorr:
        id0=id0+1
        for pGraph in pCorrFrac:
            pGraph.Write()

if __name__ == "__main__":
    print args.input
    #fPar0,fPar1,fPar2=loadCorrector(args.input        ,-41,41,1,72)
    fPar0,fPar1,fPar2=loadCorrector(args.input        ,-40,40,-1000,-1000)
    lCorr=buildTables(args.fracrange,fPar0,fPar1,fPar2,-41,41,1,72,args.ptrange)
    writeCorr(lCorr)
