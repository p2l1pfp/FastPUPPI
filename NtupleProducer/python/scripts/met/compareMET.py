import sys
import ROOT
import tdrstyle
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadRightMargin(0.10);
ROOT.gStyle.SetOptFit(0 );
ROOT.gStyle.SetPalette(55);
import math
from array import array
import os, sys
from METContainer import METContainer


##---------------------------------------
def main():

	idir = './'

	bkg =  ROOT.TFile(idir + "/jetmetTuple_DYJetsToLL_PU140.root")
	sig = ROOT.TFile(idir + "/jetmetTuple_TTbar_PU140_looseTkQuality.root")

	#metTypes = ["Calo"]
	metTypes = ["METCalo","METTK","METPuppi","METTKV","METPF"];
	#ttags = [];

	bkgMets = []
	bkgMetCuts_1percent = []
	sigMets = []
	for t in metTypes:
	 bkgMets.append( METContainer(bkg.Get("ntuple/tree"),t,'',400) );
	for b in bkgMets: 
	 bkgMetCuts_1percent.append( b._L1Met_1percent );
	 print b._tag, b._L1Met_1percent

	mctr = 0;
	for t in metTypes: 
	  sigMets.append( METContainer(sig.Get("ntuple/tree"),t,'',400,bkgMetCuts_1percent[mctr])  );
	  mctr += 1;

	rocs = [];
	roclegs = [];
	for m in metTypes:
	 curtype = m;
	 if m == "": curtype = "pf"
	 curmet_sig = getMetByTag(sigMets, curtype);
	 curmet_bkg = getMetByTag(bkgMets, curtype);
	 h = [];
	 h.append( curmet_sig.h_met_sig100 )
	 h.append( curmet_bkg.h_met )
	 leg = ['ttbar','DY']
	 makeCanvases(h,leg, "met_"+curtype);
	 print "making ROC: ", curtype
	 cur_tg = makeROCFromHisto(h);
	 cur_tg.SetName("roc"+curtype)
	 rocs.append( cur_tg );
	 roclegs.append( curtype );

	print "number of rocs = ", len(rocs);
	makeCanvasGraphsROCs(rocs,roclegs,"rocs",1e-3,1e-3,1,1,True);

	gr_turnons = [];
	gr_turnons_legs = [];
	for i,s in enumerate(sigMets):
		gr_turnons.append( s._MetTurnOn );
		curleg = convertleg(s._xtag[:-3]) + str(int(bkgMetCuts_1percent[i])) + " GeV)"

		gr_turnons_legs.append( curleg );
	makeCanvasGraphs( gr_turnons, gr_turnons_legs, "turnons",0,0,500,1,False)

### ------------ utils ------------
### ------------ utils ------------

def convertleg( s ): 
	if   s == "pf_pt2": return "PF (pT_{tk} = 2; "
	elif s == "pf_pt3": return "PF (pT_{tk} = 3; "
	elif s == "tk_pt2": return "TKs (pT_{tk} = 2; "
	elif s == "tk_pt3": return "TKs (pT_{tk} = 3; "
	elif s == "ucalo_pt2": return "Calo (pT_{tk} = 2; "
	elif s == "ucalo_pt3": return "Calo (pT_{tk} = 3; "
	elif s == "pup_pt2": return "PUPPI (pT_{tk} = 2; "
	elif s == "pup_pt3": return "PUPPI (pT_{tk} = 3; "
	elif s == "pvtk_pt2": return "TKs_{PV} (pT_{tk} = 2; "
	elif s == "pvtk_pt3": return "TKs_{PV} (pT_{tk} = 3; "
	else: return s;
### ------------ utils ------------

def getMetByTag(mets,tag):
	themet = None
	for m in mets: 
		if tag in m._xtag: 
			themet = m;
			break;
	return themet;

def makeROCFromHisto(hists,LtoR=True):

	hsig = hists[0];
	hbkg = hists[1];

	nbins = hsig.GetNbinsX();
	binsize = hsig.GetBinWidth(1);
	lowedge = hsig.GetBinLowEdge(1);

	#print "lowedge: ",lowedge

	hsigIntegral = hsig.Integral();
	hbkgIntegral = hbkg.Integral();

	xval = array('d', [])
	yval = array('d', [])
	ctr = 0;
	effBkgPrev = -9999;
	for i in range(1,nbins+1):

			effBkg = 0;
			effSig = 0;

			if LtoR: effBkg = hbkg.Integral( i, nbins )/hbkgIntegral;
			else: effBkg = hbkg.Integral( 1, i )/hbkgIntegral;

			if LtoR: effSig = hsig.Integral( i, nbins )/hsigIntegral;
			else: effSig = hsig.Integral( 1, i )/hsigIntegral;

			#if not effBkg == 0.: print "cut: ",(lowedge+(i-1)*binsize),"effBkg: ", effBkg, ", effSig: ", effSig;

			xval.append( effSig );
			yval.append( effBkg );

			if effBkg > 0.005 and effBkg < 0.015: print effBkg, hbkg.GetBinCenter(i)

			#effBkgPrev = effBkg;
			ctr = ctr + 1;

	#print nbins, "and ", ctr
	tg = ROOT.TGraph( nbins, xval, yval );
	tg.SetName( "tg"+hsig.GetName() );
	return tg;	


##---------------------------------------
def makeCanvasGraphsROCs(grs,legs,name,xlo,ylo,xhi,yhi,setlog=False):

	color  = [1,2,4,6,7,3,8,14,
			  1,2,4,6,7,3,8,14]
	styles = [1,1,1,1,1,1,1,1,
			  2,2,2,2,2,2,2,2]
	c = ROOT.TCanvas("c_"+name,"c_"+name,1200,800);
	
	leg = ROOT.TLegend(0.65,0.2,0.9,0.7);
	leg.SetBorderSize(0);
	leg.SetFillStyle(0);
	leg.SetTextSize(0.030);
	leg.SetTextFont(42)
	for i in range(len(legs)):
		leg.AddEntry(grs[i],legs[i],"l");

	hrl = c.DrawFrame(xlo,ylo,xhi,yhi);
	hrl.GetYaxis().SetTitle(grs[0].GetYaxis().GetTitle());
	hrl.GetYaxis().SetTitleOffset(0.85);	
	hrl.GetXaxis().SetTitle(grs[0].GetXaxis().GetTitle());
	
	hrl.GetXaxis().SetTitle('signal efficiency');
	hrl.GetYaxis().SetTitle('bkg efficiency'); 

	y2axis = ROOT.TGaxis(xhi,ylo,xhi,yhi,0.04,40,50510,"+L G")
	y2axis.SetTitle("Rate (MHz)")
	y2axis.SetLineColor(2)
	# y2axis.SetLabelSize(0.1)
	y2axis.SetLabelOffset(0.01);
	y2axis.SetLabelColor(2)
	y2axis.SetTitleColor(2)

	for i,gr in enumerate(grs):
		gr.SetLineColor(color[i]);
		gr.SetMarkerColor(color[i]);
		gr.SetLineStyle(styles[i]);
		gr.SetLineWidth(2);
		gr.Draw("l");
	leg.Draw();
	y2axis.Draw();

	if setlog: ROOT.gPad.SetLogy();
	c.SaveAs("plots/%s.pdf" % (name));
	c.SaveAs("plots/%s.png" % (name));
	c.SaveAs("plots/%s.root" % (name));
	
##---------------------------------------
def makeCanvasGraphs(grs,legs,name,xlo,ylo,xhi,yhi,setlog=False):

	color  = [1,2,4,6,7,3,8,14,
			  1,2,4,6,7,3,8,14]
	styles = [1,1,1,1,1,1,1,1,
			  2,2,2,2,2,2,2,2]
	c = ROOT.TCanvas("c_"+name,"c_"+name,1200,800);
	
	leg = ROOT.TLegend(0.55,0.15,0.9,0.5);
	leg.SetBorderSize(0);
	leg.SetFillStyle(0);
	leg.SetTextSize(0.030);
	leg.SetTextFont(42)
	for i in range(len(legs)):
		leg.AddEntry(grs[i],legs[i],"l");

	hrl = c.DrawFrame(xlo,ylo,xhi,yhi);
	hrl.GetYaxis().SetTitle(grs[0].GetYaxis().GetTitle());
	hrl.GetYaxis().SetTitleOffset(0.85);	
	hrl.GetXaxis().SetTitle(grs[0].GetXaxis().GetTitle());
	
	hrl.GetXaxis().SetTitle('GEN MET');
	hrl.GetYaxis().SetTitle('signal efficiency (100 kHz bkg rate)'); 

	# y2axis = ROOT.TGaxis(xhi,ylo,xhi,yhi,0.04,40,50510,"+L G")
	# y2axis.SetTitle("Rate (MHz)")
	# y2axis.SetLineColor(2)
	# # y2axis.SetLabelSize(0.1)
	# y2axis.SetLabelOffset(0.01);
	# y2axis.SetLabelColor(2)
	# y2axis.SetTitleColor(2)

	for i,gr in enumerate(grs):
		gr.SetLineColor(color[i]);
		gr.SetMarkerColor(color[i]);
		gr.SetLineStyle(styles[i]);
		gr.SetLineWidth(2);
		gr.Draw("l");
	leg.Draw();
	# y2axis.Draw();

	if setlog: ROOT.gPad.SetLogy();
	c.SaveAs("plots/%s.pdf" % (name));
	c.SaveAs("plots/%s.png" % (name));	
	c.SaveAs("plots/%s.root" % (name));	

##---------------------------------------
def makeCanvases(hists,legs,name,normalize=False):

	color = [1,2,4,6,7,3,5,8]
	options = ["hist",
			   "histsames",
			   "histsames",
			   "histsames",
			   "histsames",
			   "histsames",
			   "histsames",
			   "histsames",
			   "histsames"]

	leg = ROOT.TLegend(0.65,0.7,0.9,0.9);
	leg.SetBorderSize(0);
	leg.SetFillStyle(0);
	leg.SetTextSize(0.035);
	for i in range(len(hists)):
		leg.AddEntry(hists[i],legs[i],"l");

	c = ROOT.TCanvas("c"+hists[0].GetName(),"c"+hists[0].GetName(),1000,800);

	max = -999;

	for i in range(len(hists)):
		hists[i].SetLineColor(color[i])
		hists[i].SetMarkerColor(color[i])
		hists[i].SetMarkerStyle(20)

		if hists[i].GetMaximum() > max: 
			max = hists[i].GetMaximum();
			hists[0].SetMaximum(max*1.25);
		if normalize and hists[i].Integral() > 0: hists[i].Scale(1./hists[i].Integral())

		hists[i].Draw(options[i]);

	leg.Draw();
	c.SaveAs("plots/"+name+".pdf")
	c.SaveAs("plots/"+name+".png")
	c.SaveAs("plots/"+name+".root")
	ROOT.gPad.SetLogy();
	c.SaveAs("plots/"+name+"_log.pdf")
	c.SaveAs("plots/"+name+"_log.png")
	c.SaveAs("plots/"+name+"_log.root")

def makeAGraph(listx,listy,linecolor = 1, linestyle = 1):

        a_m = array('d', []);
        a_g = array('d', []);

        for i in range(len(listx)):
                a_m.append(listx[i]);
                a_g.append(listy[i]);

        gr = ROOT.TGraph(len(listx),a_m,a_g);

        gr.SetLineColor(linecolor)
        gr.SetLineStyle(linestyle)
        gr.SetLineWidth(3)

        return gr

# ##---------------------------------------
# def makeCanvas_2D(h,name,logy=False,tag=""):

# 	c = ROOT.TCanvas("c_"+h.GetName(),"c_"+h.GetName(),1000,800);
# 	h.SetLineWidth(2);
# 	h.Draw("colz");
# 	ROOT.gPad.SetLogz();
# 	if logy: ROOT.gPad.SetLogy();
# 	c.SaveAs("plots_%s/%s.pdf" % (tag,name));

##---------------------------------------
if __name__ == '__main__':
		main();
##---------------------------------------     
