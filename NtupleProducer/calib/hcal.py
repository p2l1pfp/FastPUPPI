#! /usr/bin/env python
import ROOT as r,sys,math,array,os,argparse,numpy
from array import array
sys.path.insert(0, '/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_9_1_0_pre3/src/FastPUPPI/NtupleProducer/calib')
from tools import end
from ecal  import resolution,gausRes,proj,draw
import tdrstyle
tdrstyle.setTDRStyle()

parser = argparse.ArgumentParser(description='Process benchmarks.')
parser.add_argument('-Tree'   ,'--Tree'      ,action='store'  ,dest='Tree',default='HcalInfo',help='Tree Name')
parser.add_argument('-input'  ,'--input'     ,action='store'  ,dest='input'  ,default='/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/L1PF/ChargedPion_TP.root',help='input file')
parser.add_argument('-pt'     ,'--pt'        ,dest='ptrange'  ,nargs='+',type=float,default=range(2,20,2)+range(20,100,5),help='pt range')
parser.add_argument('-build'  ,'--build'     ,action='store_true',dest='build'  ,default=False      ,help='build')
parser.add_argument('-emfrac' ,'--emfrac'    ,dest='emfrange' ,nargs='+',type=float,default=numpy.arange(-0.1,1.1,0.2),help='emf range')
parser.add_argument('-frac'   ,'--frac'      ,dest='fracrange',nargs='+',type=float,default=[x * 0.1 for x in range(0, 11)],help='em fraction range')
parser.add_argument('-eta'    ,'--eta'       ,dest='ieta'     ,nargs='+',type=int  ,default=[10],help='ieta')
args = parser.parse_args()

fColor=[r.kRed,r.kBlue,r.kGreen+1,r.kViolet,r.kRed-2,r.kBlue-2,r.kSpring,r.kAzure+1,r.kOrange]

def surfaceFitOld(iTree,iPt,iEMF,iEta):
    lGraphs=[]
    lRGraphs=[]
    for pEMF in range(0,len(iEMF)-1):
        pCut="*(ecal_et/(hcal_clust_et) > "+str(iEMF[pEMF])+" && ecal_et/(hcal_clust_et) < "+str(iEMF[pEMF+1])+")"
        print pCut,"!!!"
        pMGraph,pRGraph = resolution(iTree,"hcal_clust_et","genPt",-5,5,iPt,"(abs(hcal_ieta-"+str(iEta)+") < 2. && hcal_clust_et > 1.)"+pCut)
        lGraphs.append(pMGraph)
        lRGraphs.append(pRGraph)
    lX   = array('d', [])
    lEX  = array('d', [])
    lY0  = array('d', [])
    lEY0 = array('d', [])
    lY1  = array('d', [])
    lEY1 = array('d', []) 
    lEMF = array('d', [])
    lEEMF= array('d', [])
    for i1 in range(0,len(iEMF)-1):
        lEMF .append(iEMF[i1]*0.5+iEMF[i1+1]*0.5)
        if i1+1 != len(iEMF):
            lEEMF.append(abs(iEMF[i1]-iEMF[i1+1])*0.5) 
        else:
            lEEMF.append(abs(iEMF[i1-1]-iEMF[i1])*0.5)
    for pPt in range(0,len(iPt)-1):
        lX  .append(0.5*   (iPt[pPt]+iPt[pPt+1]))
        lEX .append(0.5*abs(iPt[pPt]-iPt[pPt+1]))
        pY0  = array('d', [])
        pEY0 = array('d', [])
        for pEMF in range(0,len(lGraphs)):
            pY0 .append(lGraphs[pEMF].GetY() [pPt])
            pEY0.append(lGraphs[pEMF].GetEY()[pPt])
        #pGraph = r.TGraphErrors(len(lEMF),lEMF,pY0,lEEMF,pEY0)
        if len(lGraphs) > 1:
            pGraph = r.TGraph(len(lEMF),lEMF,pY0)
            pGraph.Fit("pol1")
            lY0 .append(pGraph.GetFunction("pol1").GetParameter(0))
            lEY0.append(pGraph.GetFunction("pol1").GetParError (0))
            lY1 .append(pGraph.GetFunction("pol1").GetParameter(1))
            lEY1.append(pGraph.GetFunction("pol1").GetParError (1))
        else:
            lY0 .append(pY0[0])
            lEY0.append(pEY0[0])
            lY1 .append(0.)
            lEY1.append(0.)
    lMGraph = r.TGraphErrors(len(iPt)-1,lX,lY0,lEX,lEY0)
    lEGraph = r.TGraphErrors(len(iPt)-1,lX,lY1,lEX,lEY1)
    lMGraph.SetTitle("Mean_" +str(iEta))
    lEGraph.SetTitle("Slope_"+str(iEta))
    lMGraph.SetName("Mean_" +str(iEta))
    lEGraph.SetName("Slope_"+str(iEta))
    return lMGraph,lEGraph

def surfaceFit(iTree,iPt,iEMF,iEta):
    lGraphs=[]
    lRGraphs=[]
    for pEMF in range(0,len(iEMF)-1):
        pCut="*(ecal_et/(hcal_clust_et) > "+str(iEMF[pEMF])+" && ecal_et/(hcal_clust_et) < "+str(iEMF[pEMF+1])+")"
        print pCut,"!!!"
        pMGraph,pRGraph = resolution(iTree,"hcal_clust_et","genPt",-5,5,iPt,"(abs(hcal_ieta-"+str(iEta)+") < 2. && hcal_clust_et > 1.)"+pCut)
        lGraphs.append(pMGraph)
        lRGraphs.append(pRGraph)
        pMGraph.SetTitle("Mean_"+str(iEta)+"_frac"+str(pEMF))
        pRGraph.SetTitle("Res_" +str(iEta)+"_frac"+str(pEMF))
        pMGraph.SetName ("Mean_"+str(iEta)+"_frac"+str(pEMF))
        pRGraph.SetName ("Res_" +str(iEta)+"_frac"+str(pEMF))
    return lGraphs,lRGraphs

fPar0=[]
fPar1=[]
fPar2=[]

def cleanG(iG0,iG1,iEta):
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
    lMGraph.SetTitle(iG0.GetTitle())
    lEGraph.SetTitle(iG1.GetTitle())
    lMGraph.SetName (iG0.GetName())
    lEGraph.SetName (iG1.GetName())
    return lMGraph,lEGraph

def loadCorrector(iFile,iEtaMin,iEtaMax,iEMF,iBasicEMF):
    lFile = r.TFile(iFile)
    for iEta in range(iEtaMin,iEtaMax+1):
        print "Mean_" +str(iEta)
        lFrac=iEMF
        lG0=[]
        lG1=[]
        if abs(iEta) > 29:
            lFrac=iBasicEMF
        print lFrac
        for i0 in range(len(lFrac)-1):
            pG0=lFile.Get("Mean_"+str(iEta)+"_frac"+str(i0))
            pG1=lFile.Get("Res_" +str(iEta)+"_frac"+str(i0))
            #print "Mean_"+str(iEta)+"_frac"+str(i0),pG0,pG1,lFrac
            pG0,pG1=cleanG(pG0,pG1,iEta)
            lG0.append(pG0)
            lG1.append(pG1)
        fPar0.append(lG0)
        fPar1.append(lG1)
    return fPar0,fPar1,fPar2

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

def buildTables(iFrac,iFitFrac,iPar0,iPar1,iPar2,iEtaMin,iEtaMax,iPhiMin,iPhiMax,iHPt):
    lCorr = []
    pId = -1
    for iEta     in range(iEtaMin,iEtaMax+1):
        if abs(iEta) != 41:
            pId=pId+1 
        for iPhi in range(iPhiMin,iPhiMax+1):
            pEta = -1
            lCorrFrac=[]
            #print "Eta",iEta,"Phi",iPhi
            for i1 in range(0,len(iFrac)):
                pFrac=0
                if len(fPar0[pId]) > 1:
                    for i2 in range(0,len(iFitFrac)):
                        if iFitFrac[i2] > iFrac[i1]:
                            pFrac=i2-1
                            break
                lX   = array('d', [])
                lY   = array('d', [])
                for i2 in range(0,len(iHPt)-1):
                    lA  = fPar0[pId][pFrac].Eval(iHPt[i2])
                    pPt = iHPt[i2]*lA + iHPt[i2]
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
    print args.Tree,args.input
    lBasicEMF = [-1.,1.1]
    if args.build:
        fPar0,fPar1,fPar2=loadCorrector(args.input        ,-40,40,args.emfrange,lBasicEMF)
        lCorr=buildTables(args.fracrange,args.emfrange,fPar0,fPar1,fPar2,-41,41,1,72,args.ptrange)
        writeCorr(lCorr)
        exit()
    lFile = r.TFile.Open(args.input)
    lTree = lFile.Get(args.Tree)
    lMu   = []
    lRes  = []
    for iEta in args.ieta:
        print "Eta",iEta
        if abs(iEta) > 29:
            pMu,pRes=surfaceFit(lTree,args.ptrange,lBasicEMF,iEta)
            lMu.extend(pMu)
            lRes.extend(pRes)
        else:
            pMu,pRes=surfaceFit(lTree,args.ptrange,args.emfrange,iEta)
            lMu.extend(pMu)
            lRes.extend(pRes)
    lOFile = r.TFile("Output.root","RECREATE")
    lMu.extend(lRes)
    for pGraph in lMu:
        pGraph.Write()
    lOFile.Close()
    os.system('mv Output.root /afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_9_1_0_pre3/src/FastPUPPI/NtupleProducer/calib/hcal/hcal_ieta_%d.root' % args.ieta[0])
