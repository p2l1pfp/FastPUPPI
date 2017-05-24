import os
import ROOT
from array import array
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
    parser.add_option("--xpt", "--xpt", dest="xpt",  default=("mc_pt","p_{T}^{gen}"), nargs=2)
    parser.add_option("--coarse-eta", dest="coarseEta", default=False, action="store_true")
    parser.add_option("--root", dest="rootfile", default=None)
    parser.add_option("--eta", dest="eta", nargs=2, default=None, type="float")
    parser.add_option("--ptmin", dest="ptmin", nargs=1, default=0, type="float")
    parser.add_option("--ptbins", dest="ptbins", nargs=1, default=None, type="string")
    parser.add_option("--fitrange", dest="fitrange", nargs=2, default=(0,999), type="float")
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
    if options.rootfile:
        if options.ptbins or options.eta or "," in options.particle: raise RuntimeError
        tfout = ROOT.TFile.Open(options.rootfile, "RECREATE");
    if options.emfSlices:
        definition, strbounds = options.emfSlices
        emfbins = map(float, strbounds.split(","))
        for i,emfmax in enumerate(emfbins):
            emfmin = emfbins[i-1] if i > 0 else 0
            emfs.append((emfmax, "emf_%04.0f_%04.0f" % (emfmin*1000,emfmax*1000), "%g <= %s && min(%s,.99999) <= %g" % (emfmin, definition, definition, emfmax)))
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
        elif options.coarseEta:
            etas = [ 1.5, 3.0, 5.0 ]
            etas = [ (etas[i-1] if i else 0, etas[i]) for  i in xrange(len(etas)) if etas[i] <= options.etaMax ]
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
        if options.rootfile:
            etaedges = [e[0] for e in oldetas] + [ oldetas[-1][1] ]
            emfedges = [0] + [ e[0] for e in emfs ];
            tfout.index = ROOT.TH2F("INDEX","",len(etaedges)-1,array('f',etaedges),len(emfedges)-1,array('f',emfedges));
            tfout.WriteTObject(tfout.index)
    else:
        if options.rootfile:
            etaedges = [e[0] for e in etas] + [ etas[-1][1] ]
            tfout.index = ROOT.TH1F("INDEX","",len(etaedges)-1,array('f',etaedges));
            tfout.WriteTObject(tfout.index)
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
    if options.ptbins:
        ptedges = map(float,options.ptbins.split(","))
        ptbins = [ (ptedges[i], ptedges[i+1]) for i in xrange(len(ptedges)-1) ]
    else:
        ptbins = [ options.fitrange ]
    for iptbin,(ptmin,ptmax) in enumerate(ptbins):
      for ipoint,(oname,cut) in enumerate(sels):
        print "Plotting ",oname
        if "electron" in oname or "muon" in oname or "pi" in oname:
            cut += " && abs(mc_iso04) < 0.05" # isolated
        prof, pres = doRespPt(oname,tree,options.what,options.expr,cut,mcpt=options.mcpt,xpt=options.xpt[0],fitopt="WQ0C EX0",fitrange=(ptmin,ptmax))
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
                if scales[-1] == 0: scales[-1] = 1
        for p in prof, pres:
            if not p: continue
            if getattr(p,'fit',None):
                p.fit.SetLineWidth(2); p.fit.SetLineColor(2)
                if options.ptbins:
                    p.fit.SetRange(ptmin,ptmax)
            p.SetLineWidth(3); p.SetLineColor(1);  p.SetMarkerColor(1)
        for plot,ptype,pfix in (prof,"response","_resp"),(pres,"resolution","_resol"):
                if not plot: continue
                frame = ROOT.TH1F("stk","stk",100,0.0,250.0 if "jet" in oname else 100.0)
                frame.GetXaxis().SetTitle(options.xpt[1]+" (GeV)")
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
                    line.DrawLine(0.0,.25,frame.GetXaxis().GetXmax(),0.25)
                    if ptmin > 0:   line.DrawLine(ptmin,0.0,ptmin,1.6)
                    if ptmax < 999: line.DrawLine(ptmax,0.0,ptmax,1.6)
                plot.Draw("P SAME" if "TGraph" in plot.ClassName() else "SAME")
                if hasattr(plot,'fit'): plot.fit.Draw("SAME")
                if options.ptbins: pfix += "_ptbin%d" % (iptbin+1)
                out = odir+'/'+oname+pfix+"-"+options.what+".png"
                c1.Print(out)
                del frame
                if options.rootfile and ptype == "response":
                    etaval = 0.999*etas[ipoint][1]
                    etabin = tfout.index.GetXaxis().FindBin(etaval)
                    if options.emfSlices:
                        emfval = 0.999*emfbins[ipoint] if etaval < 3.0 else 0
                        emfbin = tfout.index.GetYaxis().FindBin(emfval)
                        print "%s is in eta bin %d (%.2f, %.2f) emf bin %d (%.2f, %.2f)" % (out, etabin, tfout.index.GetXaxis().GetBinLowEdge(etabin), tfout.index.GetXaxis().GetBinUpEdge(etabin), emfbin, tfout.index.GetYaxis().GetBinLowEdge(emfbin),  tfout.index.GetYaxis().GetBinUpEdge(emfbin))
                        pclone = plot.Clone("eta_bin%d_emf_bin%d" % (etabin,emfbin))
                    else:
                        print "%s is in eta bin %d (%.2f, %.2f)" % (out, etabin, tfout.index.GetXaxis().GetBinLowEdge(etabin), tfout.index.GetXaxis().GetBinUpEdge(etabin))
                        pclone = plot.Clone("eta_bin%d" % etabin)
                    ## simple inversion, computed with reco pt, doesn't close well
                    #for i in xrange(pclone.GetN()):
                    #    ratio = pclone.GetY()[i]
                    #    corr = 1 if ratio == 0 else 1/ratio
                    #    if   corr > 4:   corr = 4
                    #    elif corr < 0.5: corr = 0.5
                    #    pclone.SetPoint(i, pclone.GetX()[i], corr);
                    ## more fancy inversion
                    # make a clean reco vs gen plot
                    for i in xrange(pclone.GetN()):
                        pclone.SetPoint(i, pclone.GetX()[i], pclone.GetX()[i]*pclone.GetY()[i]);
                    # invert the plot
                    for i in xrange(pclone.GetN()):
                        pclone.SetPoint(i, pclone.GetY()[i], pclone.GetX()[i]);
                    pclone.Sort()
                    #for i in xrange(pclone.GetN()):
                    #    corr = pclone.GetY()[i]/pclone.GetX()[i] if pclone.GetX()[i] else 1
                    #    pclone.SetPoint(i, pclone.GetX()[i], );
                    tfout.WriteTObject(pclone)
    if options.ptbins:
        # convert to string
        sptMins = ["% 6.2f" % b[0] for b in ptbins]
        sptMaxs = ["% 6.2f" % b[1] for b in ptbins]
        # extend the last bin
        sptMaxs[-1] = "%6.0f" % (999999.,)
        # then make copies
        sptMins = [b for b in sptMins for s in sels]
        sptMaxs = [b for b in sptMaxs for s in sels]
        # and make copies of etas
        etas    = [e for b in ptbins for e in etas]
        if options.emfSlices:
            emfbins = [e for b in ptbins for e in emfbins]
    if options.rootfile:
        print "Calibrations stored in %s" % options.rootfile
        tfout.ls()
        tfout.Close()
    elif options.emfSlices:
        print "simpleCorrHad = cms.PSet( "
        print "\t\t\tetaBins = cms.vdouble(%s),"  % (", ".join("% 6.3f" % s[1] for s in etas))
        print "\t\t\temfBins = cms.vdouble(%s),"  % (", ".join("% 6.3f" % s for s in emfbins))
        print "\t\t\toffset  = cms.vdouble(%s),"  % (", ".join("% 6.3f" % s for s in offsets))
        print "\t\t\tscale   = cms.vdouble(%s),"  % (", ".join("% 6.3f" % s for s in scales))
        if options.ptbins:
            print "\t\t\tptMin   = cms.vdouble(%s),"  % (", ".join(sptMins))
            print "\t\t\tptMax   = cms.vdouble(%s),"  % (", ".join(sptMaxs))
        print ")"
    else:
        print "simpleCorrEm = cms.PSet( "
        print "\t\t\tetaBins = cms.vdouble(%s),"  % (", ".join("% 6.3f" % s[1] for s in etas))
        print "\t\t\toffset  = cms.vdouble(%s),"  % (", ".join("% 6.3f" % s for s in offsets))
        print "\t\t\tscale   = cms.vdouble(%s),"  % (", ".join("% 6.3f" % s for s in scales))
        if options.ptbins:
            print "\t\t\tptMin   = cms.vdouble(%s),"  % (", ".join(sptMins))
            print "\t\t\tptMax   = cms.vdouble(%s),"  % (", ".join(sptMaxs))
        if options.resolution:
            kind = "track" if "TK" in options.expr and "ele" not in oname else "calo"
            print "\t\t\tkind    = cms.string('%s'),"  % (kind)
        print ")"
