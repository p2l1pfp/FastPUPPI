import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

from FastPUPPI.NtupleProducer.scripts.respPlots import doRespPt

##### EXAMPLE USAGE:
# 1) ECAL Corrections from Pi0
#       python scripts/respCorrSimple.py respTupleNew_SinglePion0_NoPU.root plots/corr/em -p pizero -w TPEcal_pt02 -e TPEcal_pt02 -E 3.0
#
#    Check closure, after re-running Pi0 with Ecal calibration applied
#       python scripts/respCorrSimple.py respTupleNew_SinglePion0_NoPU_C.root plots/corr/em/closure -p pizero -w L1Ecal_pt02 -e L1Ecal_pt02 -E 3.0
#       
#       
# 2) Hadron Corrections, after re-running the pions with the Ecal corrections applied
#       python scripts/respCorrSimple.py respTupleNew_SinglePion_NoPU_C.root plots/corr/had -p pion -w L1RawCalo_pt02 -e L1RawCalo_pt02 --emf-slices L1Ecal_pt02/L1RawCalo_pt02 0.125,0.5,0.875
#
#    Check closure, after re-running Pi with Ecal + Had calibration applied 
#       python scripts/respCorrSimple.py respTupleNew_SinglePion_NoPU_C.root plots/corr/had/closure -p pion -w L1Calo_pt02 -e L1Calo_pt02 --emf-slices L1Ecal_pt02/L1RawCalo_pt02 0.125,0.5,0.875
#   
#    
# 3) Resolutions, after all the corrections
#    Hadrons, calo
#        python scripts/respCorrSimple.py respTupleNew_SinglePion_NoPU_C.root plots/corr/resolution -p pion -w L1Calo_pt02 -e L1Calo_pt02  -r
#    EM, calo
#        python scripts/respCorrSimple.py respTupleNew_SinglePion0_NoPU_C.root plots/corr/resolution -p pizero -w L1Calo_pt02 -e L1Calo_pt02  -r
#    Hadrons, track
#        python scripts/respCorrSimple.py respTupleNew_SinglePion_NoPU_C.root plots/corr/resolution -p pion -w TPTK_pthighest -e TPTK_pthighest   -r

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser("%(prog) infile [ src [ dst ] ]")
    parser.add_option("-w", dest="what",  default=None)
    parser.add_option("-e", dest="expr",  default=None)
    parser.add_option("-p", dest="particle", default=None, help="Choose particle (electron, ...)")
    parser.add_option("-E", "--etaMax", dest="etaMax",  default=5.0, type=float)
    parser.add_option("--mc", "--mc", dest="mcpt",  default="mc_pt")
    parser.add_option("--eta", dest="eta", nargs=2, default=None, type="float")
    parser.add_option("--ptmin", dest="ptmin", nargs=1, default=0, type="float")
    parser.add_option("--emf-slices", dest="emfSlices", nargs=2, default=None)
    parser.add_option("-r","--resolution", dest="resolution", default=False, action="store_true")
    options, args = parser.parse_args()
    selparticles = options.particle.split(",") if options.particle else []
    sels = []; 
    offsets = []
    scales  = []
    etas    = []
    emfbins = []
    emfs    = []
    if options.emfSlices:
        definition, strbounds = options.emfSlices
        emfbins = map(float, strbounds.split(","))
        for i,emfmax in enumerate(emfbins):
            emfmin = emfbins[i-1] if i > 0 else 0
            emfs.append((emfmax, "emf_%04.0f_%04.0f" % (emfmin*1000,emfmax*1000), "%g <= %s && %s <= %g" % (emfmin, definition, definition, emfmax)))
    for (particle, pdgIdCut, minPt, maxEta) in [ 
            ("pion", "abs(mc_id) == 211", 2, 5),
            ("pizero", "abs(mc_id) == 111", 2, 5),
            ("pimix", "(abs(mc_id) == 211 || (abs(mc_id) == 111 && (event % 2) == 1))", 2, 5),
            ("photon", "abs(mc_id) == 22", 10, 5),
            ("electron", "abs(mc_id) == 11", 10, 5),
            ("muon", "abs(mc_id) == 13", 10, 5),
            ("tau", "(abs(mc_id) == 15 || abs(mc_id) == 211)", 20, 5),
            ("jet", "abs(mc_id) == 0", 20, 5),
        ]:
        if options.particle != particle: continue
        if options.eta:
            etas = [ options.eta ]
        elif options.resolution:
            if "TK" in options.expr:
                etas = [ 0.8, 1.2, 1.5, 2.0, 2.5 ]
            elif ("PF" in options.expr or "Puppi" in options.expr):
                etas = [ 1.3, 1.7, 2.5, 3.0, 4.0, 5.0 ]
            else:
                etas = [ 1.3, 1.7, 2.8, 3.2, 4.0, 5.0 ]
            etas = [ (etas[i-1] if i else 0, etas[i]) for  i in xrange(len(etas)) ]
        else: 
            for ieta in range(0,50,5):
                if ieta*0.1 >= options.etaMax: break
                etas.append((0.1*(ieta),0.1*(ieta+5)))
        for etamin,etamax in etas:
            sels.append(("%s_eta_%02d_%02d"  % (particle,int(etamin*10),int(etamax*10)), "abs(mc_eta) > %.1f && abs(mc_eta) < %.1f && mc_pt > %g &&  %s" % (etamin,etamax,max(options.ptmin,minPt),pdgIdCut)))
    if options.emfSlices:
        oldsels, oldetas = sels[:], etas[:]
        sels = []; etas = []; emfbins = []
        for sel,(etamin,etamax) in zip(oldsels,oldetas):
            if abs(etamin) < 3.0:
                for emfmax,emfname,emfcut in emfs:
                    sels.append((sel[0]+"_"+emfname, "(%s) && (%s)" % (sel[1],emfcut)))
                    etas.append((etamin,etamax))
                    emfbins.append(emfmax)
            else:
                sels.append(sel)
                etas.append((etamin,etamax))
                emfbins.append(1.1)
    tree = ROOT.TChain("ntuple/tree")
    for a in args[:]:
        if a.endswith(".root"): 
            tree.Add(a)
            args.remove(a)
    print args
    odir = args[0] # "plots/910pre2/test"
    os.system("mkdir -p "+odir)
    os.system("cp %s/src/FastPUPPI/NtupleProducer/python/display/index.php %s/" % (os.environ['CMSSW_BASE'], odir));
    ROOT.gROOT.ProcessLine(".x %s/src/FastPUPPI/NtupleProducer/python/display/tdrstyle.cc" % os.environ['CMSSW_BASE']);
    ROOT.gStyle.SetOptStat(False)
    ROOT.gStyle.SetOptFit(False)
    ROOT.gStyle.SetErrorX(0.5)
    c1 = ROOT.TCanvas("c1","c1")
    for oname,cut in sels:
        print "Plotting ",oname
        if "electron" in oname or "muon" in oname or "pi" in oname:
            cut += " && abs(mc_iso04) < 0.05" # isolated
        prof, pres = doRespPt(oname,tree,options.what,options.expr,cut,mcpt=options.mcpt,fitopt="WQ0C EX0")
        if options.resolution:
            if not pres or not getattr(pres,'fit',None):
                offsets.append(8.0)
                scales.append(0.2)
                continue
            else:
                offsets.append(pres.fit.GetParameter(0))
                scales.append(pres.fit.GetParameter(1))
        else:
            if not prof or not getattr(prof,'fit',None):
                offsets.append(0.)
                scales.append(1.)
                continue
            else:
                offsets.append(prof.fit.GetParameter(0))
                scales.append(prof.fit.GetParameter(1))
        for p in prof, pres:
            if not p: continue
            if getattr(p,'fit',None):
                p.fit.SetLineWidth(2); p.fit.SetLineColor(1)
            p.SetLineWidth(3); p.SetLineColor(1);  p.SetMarkerColor(1)
        for plot,ptype,pfix in (prof,"response","_resp"),(pres,"resolution","_resol"):
                if not plot: continue
                frame = ROOT.TH1F("stk","stk",100,0.0,250.0 if "jet" in oname else 100.0)
                frame.GetXaxis().SetTitle("p_{T} (GeV)")
                frame.GetYaxis().SetDecimals(True)
                if "resolution" in ptype:
                    frame.GetYaxis().SetTitle("#sigma_{eff}(p_{T}^{corr})/p_{T}^{corr}")
                    frame.GetYaxis().SetRangeUser(0.008,10.0)
                    c1.SetLogy(True)
                else:
                    frame.GetYaxis().SetTitle("median p_{T}^{rec}/p_{T}^{gen}")
                    frame.GetYaxis().SetRangeUser(0,1.6)
                    c1.SetLogy(False)
                frame.Draw()
                if "resolution" not in ptype:
                    line = ROOT.TLine()
                    line.SetLineStyle(7)
                    line.DrawLine(0.0,1,frame.GetXaxis().GetXmax(),1)
                plot.Draw("P SAME" if "TGraph" in plot.ClassName() else "SAME")
                if hasattr(plot,'fit'): plot.fit.Draw("SAME")
                out = odir+'/'+oname+pfix+"-"+options.what+".png"
                c1.Print(out)
                del frame
    if options.emfSlices:
        print "simpleCorrHad = cms.PSet( "
        print "\t\t\tetaBins = cms.vdouble(%s),"  % (", ".join("% 6.3f" % s[1] for s in etas))
        print "\t\t\temfBins = cms.vdouble(%s),"  % (", ".join("% 6.3f" % s for s in emfbins))
        print "\t\t\toffset  = cms.vdouble(%s),"  % (", ".join("% 6.3f" % s for s in offsets))
        print "\t\t\tscale   = cms.vdouble(%s),"  % (", ".join("% 6.3f" % s for s in scales))
        print ")"
    else:
        print "simpleCorrEM = cms.PSet( "
        print "\t\t\tetaBins = cms.vdouble(%s),"  % (", ".join("% 6.3f" % s[1] for s in etas))
        print "\t\t\toffset  = cms.vdouble(%s),"  % (", ".join("% 6.3f" % s for s in offsets))
        print "\t\t\tscale   = cms.vdouble(%s),"  % (", ".join("% 6.3f" % s for s in scales))
        if options.resolution:
            kind = "track" if "TK" in options.expr and "ele" not in oname else "calo"
            print "\t\t\tkind    = cms.string('%s'),"  % (kind)
        print ")"
