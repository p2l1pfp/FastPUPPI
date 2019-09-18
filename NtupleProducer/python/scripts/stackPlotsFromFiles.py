import os, re, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(False)
ROOT.gStyle.SetErrorX(0.5)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


from optparse import OptionParser
parser = OptionParser("%(prog) out what [name1=]file1[,color1]  [name1=]file2[,color2] ...  ")
parser.add_option("--logy", dest="logy", action="store_false", default=False, help="Log y axis")
parser.add_option("--liny", dest="logy", action="store_false", default=False, help="Linear y axis")
parser.add_option("--legend", dest="legend", default="BR", help="Legend placement: BR, TR")
#parser.add_option("-W", dest="what_reg",     default=None, help="Choose set (inputs, l1pf, ...)")
#parser.add_option("-P","--plots", dest="plots", default="rate,isorate,roc,effc", help="Choose plot or comma-separated list of plots")
#parser.add_option("-J","--jets", dest="jets", default="ak4", help="Choose jet flavour")
#parser.add_option("-j","--jecs", dest="jecs", default="jecs.root", help="Choose JEC file")
#parser.add_option("--jm","--jec-method", dest="jecMethod", default="", help="Choose JEC method")
#parser.add_option("-R","--raw", dest="rawJets", default=False, action="store_true", help="Don't appy JECs")
#parser.add_option("-s", dest="genht",  default=None, type="float", help="Choose gen ht")
#parser.add_option("-r", dest="rate",  default="10,20,50", type="string", help="Choose rate [kHz] (for isorate plots, can specify more than one)")
#parser.add_option("-l","--label", dest="label",  default=None, type="string", help="Extra label for plots")
#parser.add_option("-p", "--pt", dest="pt",  default=30, type="float", help="Choose pt cut")
#parser.add_option("-e", "--eta", dest="eta",  default=2.4, type="float", help="Choose eta")
#parser.add_option("-v", dest="var",  default="ht", help="Choose variable (ht, met, metCentral, mht, jet<N>, mjj, ptj-mjj<M>)")
#parser.add_option("--xlabel","--varlabel", dest="varlabel", default=None, help="X axis label for the variable")
#parser.add_option("--xmax", dest="xmax",  default=None, type=float, help="Choose variable")
#parser.add_option("--logxbins", dest="logxbins",  default=None, nargs=2, type=float, help="--logxbins N X will make N bins, the last being a factor X larger than the first")
options, args = parser.parse_args()

if len(args) < 4:
    parser.print_usage()
    sys.exit(1)

colors = [ ROOT.kBlack, ROOT.kRed, ROOT.kGreen+2, ROOT.kBlue, ROOT.kMagenta+1, ROOT.kOrange+7, ROOT.kCyan+1, ROOT.kGray+2, ROOT.kViolet+5, ROOT.kSpring+5, ROOT.kAzure+1, ROOT.kPink+7, ROOT.kOrange+3, ROOT.kBlue+3, ROOT.kMagenta+3, ROOT.kRed+2 ]; colors.reverse()
tfiles = []
tobjs = []
labels = []
frame = None
for i,fname in enumerate(args[2:]):
    if "," in fname:
        (fname,colorname) = fname.split(",")
        color = eval("ROOT.k"+colorname)
    else:
        color = colors.pop()
    if "=" in fname: 
        (name, fname) = fname.split("=")
    else: 
        name = fname.replace(".root","")
    tfile = ROOT.TFile.Open(fname)
    if not tfile: 
        print "ERROR opening %s" % tfile
        continue
    tobj = tfile.Get(args[1])
    if not tobj:
        print "ERROR fetching %r from %s" % (args[1], tfile)
        tfile.ls()
        continue
    if not frame:
        frame = tfile.Get("frame")
        if not frame:
            print "ERROR fetching %r from %s" % ("frame", tfile)
            tfile.ls()
    tfiles.append(tfile)
    tobjs.append(tobj)
    labels.append(name)
    tobj.SetLineColor(color)
    tobj.SetMarkerColor(color)
if not tobjs:
    print "Nothing to plot!"
    sys.exit(1)

odir = os.path.dirname(args[0])
if odir:
    if not os.path.isdir(odir):
        os.system("mkdir -p "+odir)
    os.system("cp %s/src/FastPUPPI/NtupleProducer/python/display/index.php %s/" % (os.environ['CMSSW_BASE'], odir));
ROOT.gROOT.ProcessLine(".x %s/src/FastPUPPI/NtupleProducer/python/display/tdrstyle.cc" % os.environ['CMSSW_BASE']);
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
c1 = ROOT.TCanvas("c1","c1")

c1.SetLogy(options.logy)
if options.legend == "TR":
    leg = ROOT.TLegend(0.6,0.99,0.95,0.99-0.06*len(labels))
elif options.legend == "BR":
    leg = ROOT.TLegend(0.6,0.19,0.95,0.19+0.06*len(labels))
else:
    raise RuntimeError("Unsupported legend option %r" % options.legend)

frame.Draw()
for n,p in zip(labels,tobjs): 
    p.Draw("PCX SAME" if ("TH1" not in p.ClassName()) else "C SAME")
    leg.AddEntry(p, n, "LP")
leg.Draw()

c1.Print(args[0])
print "Wrote %s" % args[0]
