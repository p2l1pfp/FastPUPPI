import os, sys
import ROOT
from array import array
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

from FastPUPPI.NtupleProducer.scripts.respPlots import doRespPt

##### EXAMPLE USAGE:
# 1) ECAL Corrections from Pi0
#       python scripts/respCorrSimple.py respTupleNew_HadronGun_PU0.root plots/em-barrel -p pizero -w L1RawEcal_pt02 -e L1RawEcal_pt02 --fitrange 20 50 -E 3.0 --ptmax 50 --root ecorr.root
#
#    Check closure, after re-running Pi0 with Ecal calibration applied
#       python scripts/respCorrSimple.py respTupleNew_SinglePion0_NoPU_C.root plots/corr/em/closure -p pizero -w L1Ecal_pt02 -e L1Ecal_pt02 -E 3.0
#       
#       
# 2) Hadron Corrections, after re-running the pions with the Ecal corrections applied
#       python scripts/respCorrSimple.py respTupleNew_HadronGun_PU0.root plots/corr/had -p piswitch -w L1RawCalo_pt02 -e L1RawCalo_pt02 --fitrange 15 50 --root hadcorr.root --emf-slices L1RawCaloEM_pt02/L1RawCalo_pt02 0.125,0.50,0.875,1.0  --ptmax 50
#
#    Check closure, after re-running Pi with Ecal + Had calibration applied 
#  
#    python scripts/respPlots.py respTupleNew_HadronGun_PU0.root plots/corr/ -w l1pf -p pizero,pimix,pion --no-resol
#    
# 3) Resolutions, after all the corrections
#    Hadrons, calo
#        python scripts/respCorrSimple.py respTupleNew_HadronGun_PU0.root plots/corr/resol -p pion -w L1Calo_pthighest -e L1Calo_pthighest  -r --ptmin 5 --ptmax 50
#    EM, calo
#        python scripts/respCorrSimple.py respTupleNew_HadronGun_PU0.root plots/corr/resol -p pizero -w L1Ecal_pt02 -e L1Ecal_pt02  -r --ptmin 5 --ptmax 50
#    Hadrons, track
#        python scripts/respCorrSimple.py respTupleNew_HadronGun_PU0.root plots/corr/resol -p pion -w L1TK_pthighest -e L1TK_pthighest   -r --ptmin 5 --ptmax 50


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser("%(prog) infile [ src [ dst ] ]")
    parser.add_option("-w", dest="what",  default=None)
    parser.add_option("-e", dest="expr",  default=None)
    parser.add_option("-p", dest="particle", default=None, help="Choose particle (electron, ...)")
    parser.add_option("--cut", dest="extracut", default=None, help="Extra cut", nargs=2)
    parser.add_option("-E", "--etaMax", dest="etaMax",  default=5.0, type=float)
    parser.add_option("--mc", "--mc", dest="mcpt",  default="mc_pt")
    parser.add_option("--xpt", "--xpt", dest="xpt",  default=("mc_pt","p_{T}^{gen}"), nargs=2)
    parser.add_option("--coarse-eta", dest="coarseEta", default=False, action="store_true")
    parser.add_option("--semicoarse-eta", dest="semiCoarseEta", default=False, action="store_true")
    parser.add_option("--root", dest="rootfile", default=None)
    parser.add_option("--eta", dest="eta", nargs=2, default=None, type="float")
    parser.add_option("--ptmin", dest="ptmin", nargs=1, default=0, type="float")
    parser.add_option("--ptmax", dest="ptmax", nargs=1, default=0, type="float")
    parser.add_option("--ptbins", dest="ptbins", nargs=1, default=None, type="string")
    parser.add_option("--fitrange", dest="fitrange", nargs=2, default=(0,999), type="float")
    parser.add_option("--emf-slices", dest="emfSlices", nargs=2, default=None)
    parser.add_option("-r","--resolution", dest="resolution", default=False, action="store_true")
    parser.add_option("--no-resol", dest="noResol", default=False, action="store_true", help="skip resolution plots")
    parser.add_option("-g","--gauss", dest="gauss", default=False, action="store_true", help="make also gaussian estimates")
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
            ("pion", "abs(mc_id) == 211", 2, 2),
            ("pizero", "abs(mc_id) == 111", 2, 5),
            ("klong", "mc_id == 130", 2, 5),
            ("pimix", "(abs(mc_id) == 211 || (abs(mc_id) == 111 && (event % 2) == 1))", 2, 5),
            ("piswitch", "(abs(mc_id) == 211 || (abs(mc_id) == 111 && (event % 2) == 1 && abs(mc_eta) > 2.5))", 2, 5),
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
        elif options.semiCoarseEta:
            #etas = [ 1.5, 3.0, 3.5, 3.8, 4.2, 4.4, 4.7, 5.0 ]
	    etas = [1.5, 3.0, 3.2, 3.5, 4.0, 4.5, 5.0]
            etas = [ (etas[i-1] if i else 0, etas[i]) for  i in xrange(len(etas)) if etas[i] <= options.etaMax ]
        elif options.resolution:
            if "TK" in options.expr:
                etas = [ 0.8, 1.2, 1.5, 2.0, 2.5 ]
            elif ("PF" in options.expr or "Puppi" in options.expr):
                etas = [ 1.3, 1.7, 2.5, 3.0, 4.0, 5.0 ]
            else:
                etas = [ 1.3, 1.7, 2.8, 3.2, 4.0, 5.0 ]
            etas = [ (etas[i-1] if i else 0, etas[i]) for  i in xrange(len(etas)) if etas[i] <= options.etaMax ]
        else: 
            for ieta in range(0,50,5):
                if ieta*0.1 >= options.etaMax: break
                etas.append((0.1*(ieta),0.1*(ieta+5)))

        if options.extracut:
            particle += "_"+options.extracut[0]
            pdgIdCut += " && ("+options.extracut[1]+")"
        
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
        prof, pres = doRespPt(oname,tree,options.what,options.expr,cut,mcpt=options.mcpt,xpt=options.xpt[0],fitopt="WNQ0C EX0",fitrange=(ptmin,ptmax))
        profg = getattr(prof,'gaus',None) if prof and options.gauss else None
        presg = getattr(pres,'gaus',None) if pres and options.gauss else None
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
        for p in prof, pres, profg, presg:
            if not p: continue
            if getattr(p,'fit',None):
                p.fit.SetLineWidth(2); p.fit.SetLineColor(2)
                if options.ptbins:
                    p.fit.SetRange(ptmin,ptmax)
            p.SetLineWidth(3); p.SetLineColor(1);  p.SetMarkerColor(1)
        for plot,ptype,pfix in (prof,"response","_resp"),(pres,"resolution","_resol"),(profg,"response","_resp_gaus"),(presg,"resolution","_resol_gaus"):
                if not plot: continue
                frame = ROOT.TH1F("stk","stk",100,0.0,options.ptmax if options.ptmax else (250.0 if "jet" in oname else 100.0))
                frame.GetXaxis().SetTitle(options.xpt[1]+" (GeV)")
                frame.GetYaxis().SetDecimals(True)
                if "resolution" in ptype:
                    if options.noResol: continue
                    frame.GetYaxis().SetTitle("#sigma_{eff}(p_{T}^{corr})/p_{T}^{corr}")
                    frame.SetMaximum(1.0)
                    #frame.GetYaxis().SetRangeUser(0.008,2.0)
                    #c1.SetLogy(True)
                else:
                    frame.GetYaxis().SetTitle("median p_{T}^{rec}/p_{T}^{gen}")
                    #frame.GetYaxis().SetRangeUser(0,1.6)
                    frame.SetMaximum(1.6)
                    #c1.SetLogy(False)
                frame.Draw()
                line = ROOT.TLine()
                line.SetLineStyle(7)
                if "resolution" not in ptype:
                    line.DrawLine(0.0,1,frame.GetXaxis().GetXmax(),1)
                    line.DrawLine(0.0,.25,frame.GetXaxis().GetXmax(),0.25)
                if ptmin > 0:   line.DrawLine(ptmin,0.0,ptmin,frame.GetMaximum())
                if ptmax < 999: line.DrawLine(ptmax,0.0,ptmax,frame.GetMaximum())
                plot.Draw("P SAME" if "TGraph" in plot.ClassName() else "SAME")
                if hasattr(plot,'fit'): plot.fit.Draw("SAME")
                if options.ptbins: pfix += "_ptbin%d" % (iptbin+1)
                out = odir+'/'+oname+pfix+"-"+options.what+".png"
                c1.Print(out)
                if options.rootfile and ptype == "response":
                    tfout.WriteTObject(plot.Clone(oname),oname)
                    etaval = 0.999*etas[ipoint][1]
                    etabin = tfout.index.GetXaxis().FindBin(etaval)
                    if options.emfSlices:
                        emfval = 0.999*emfbins[ipoint] if etaval < 3.0 else 0
                        emfbin = tfout.index.GetYaxis().FindBin(emfval)
                        pname = "eta_bin%d_emf_bin%d" % (etabin,emfbin)
                    else:
                        pname = "eta_bin%d" % etabin
                   ### simple inversion
                   # pclone = plot.Clone(pname)
                   # # make a reco vs gen plot
                   # for i in xrange(pclone.GetN()):
                   #     pclone.SetPoint(i, pclone.GetX()[i], pclone.GetX()[i]*pclone.GetY()[i]);
                   # # invert the plot
                   # for i in xrange(pclone.GetN()):
                   #     pclone.SetPoint(i, pclone.GetY()[i], pclone.GetX()[i]);
                   ### fancy inversion
                    pclone = ROOT.TGraph(); pclone.SetName(pname)
                    points = []
                    for i in xrange(150,0,-1):
                        genpt = i+ 1.5
                        if genpt > options.fitrange[0]: corr = plot.fit.Eval(genpt)
                        elif genpt > plot.GetX()[0]:    corr = plot.Eval(genpt)
                        else:                           corr = plot.GetY()[0]
                        recopt = genpt * max(0.25, corr)
                        if len(points) and recopt > points[-1][0]: 
                            continue # avoid going backwards
                        points.append((recopt,genpt))    
                    pclone = ROOT.TGraph(len(points)); pclone.SetName(pname)
                    for i,(recopt,genpt) in enumerate(points):
                        pclone.SetPoint(i, recopt, genpt)
                    pclone.Sort()
                    tfout.WriteTObject(pclone)
                    presp = pclone.Clone()
                    for i in xrange(presp.GetN()):
                        presp.SetPoint(i, presp.GetX()[i], min(presp.GetY()[i]/presp.GetX()[i] if presp.GetX()[i] > 0 else 1.0, 4));
                    frame.GetXaxis().SetTitle("p_{T}^{rec} (GeV)")
                    frame.GetYaxis().SetTitle("p_{T}^{corr}/p_{T}^{rec}")
                    frame.GetYaxis().SetRangeUser(0,4.0)
                    frame.Draw()
                    presp.SetMarkerStyle(7); presp.SetLineWidth(2);
                    presp.Draw("PLX SAME")
                    line.DrawLine(0.0,1,frame.GetXaxis().GetXmax(),1)
                    out = odir+'/'+oname+pfix+"-"+options.what+"_corr.png"
                    c1.Print(out)
                    frame.GetYaxis().SetTitle("p_{T}^{corr} (GeV)")
                    frame.GetYaxis().SetRangeUser(0.0, 150.0)
                    frame.Draw() 
                    line.DrawLine(0.0,1,frame.GetXaxis().GetXmax(),150.0)
                    pextra = ROOT.TGraph(199)
                    for i in xrange(200):
                        pextra.SetPoint(i, 0.5*(i+1), pclone.Eval(0.5*(i+1)))
                    pextra.SetLineColor(ROOT.kGray+1); pextra.SetMarkerColor(ROOT.kGray+1) 
                    pextra.SetLineWidth(2); pextra.SetMarkerStyle(6); pextra.Draw("PXL SAME")
                    pclone.SetMarkerStyle(8); pclone.Draw("PX SAME")
                    tf1 = ROOT.TF1("_p1","pol1", 0,  options.ptmax if options.ptmax else 100)
                    pclone.Fit(tf1, "WNQ0C EX0", "", 20,  options.ptmax if options.ptmax else 100)
                    tf1.SetRange(0, options.ptmax if options.ptmax else 100)
                    tf1.SetLineWidth(3)
                    tf1.SetLineColor(ROOT.kGreen+2)
                    tf1.Draw("SAME")
                    out = odir+'/'+oname+pfix+"-"+options.what+"_corr1.png"
                    c1.Print(out)
                    frame.GetYaxis().SetRangeUser(0.0, 15.0)
                    frame.GetXaxis().SetRangeUser(0.0, 7.0)
                    out = odir+'/'+oname+pfix+"-"+options.what+"_corr1z.png"
                    
                    c1.Print(out)
                del frame
 
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
