import sys
import ROOT
import tdrstyle
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadRightMargin(0.08);
ROOT.gStyle.SetOptFit(0 );
ROOT.gStyle.SetPalette(55);
import math
from array import array
import os
from METContainer import METContainer


##---------------------------------------
def main():

	fin = ROOT.TFile("data_May22-offlinetracks/M_dy_140_3.0ZMM_140PU.root");
	tin = fin.Get("met");

	pfmet    = METContainer(tin,""     , "pf", 400);
	calomet  = METContainer(tin,"calo" , "calo", 400);
	tkmet    = METContainer(tin,"tk"   , "tk", 200);
	ucalomet = METContainer(tin,"ucalo", "ucalo", 400);
	puppimet = METContainer(tin,"pup"  , "pup", 100);
	tkvtxmet = METContainer(tin,"pvtk" , "pvtk", 100);
	ipfmet    = METContainer(tin,"ipf"  , "ipup", 300);
	ipuppimet = METContainer(tin,"ipup"  , "ipup", 300);

	leg = ["pf"]
	h = []; h.append( pfmet.h_mz );
	makeCanvases( h, leg, "mZ" );
	h = []; h.append( pfmet.h_ptz );
	makeCanvases( h, leg, "ptZ" );
	
	leg = ["pf","calo","tk","ucalo","puppi","tkvtx", "i-pf", "i-puppi"]
	mets = [pfmet,calomet,tkmet,ucalomet,puppimet,tkvtxmet, ipfmet, ipuppimet];
	# mets = [pfmet,calomet,tkmet,ucalomet,puppimet,tkvtxmet];
	hnames = ["h_met","h_upara","h_uparaMinusPt","h_uperp"];
	htags = ['met','upara','uparaMinusPt','uperp']
	for i,n in enumerate(hnames):
		curh = []
		for m in mets:
			curh.append( getattr(m,n) );
		makeCanvases(curh,leg,htags[i]); 

	gr = [];
	gr.append( pfmet.gr_upara_scale );
	gr.append( calomet.gr_upara_scale );
	gr.append( tkmet.gr_upara_scale );
	gr.append( ucalomet.gr_upara_scale );
	gr.append( puppimet.gr_upara_scale );
	gr.append( tkvtxmet.gr_upara_scale );
	gr.append( ipfmet.gr_upara_scale );
	gr.append( ipuppimet.gr_upara_scale );
	makeCanvasGraphs(gr,leg,"upara-scale",0,0,120,2.0);
	gr = [];
	gr.append( pfmet.gr_upara_res );
	gr.append( calomet.gr_upara_res );
	gr.append( tkmet.gr_upara_res );
	gr.append( ucalomet.gr_upara_res );
	gr.append( puppimet.gr_upara_res );	
	gr.append( tkvtxmet.gr_upara_res );	
	gr.append( ipfmet.gr_upara_res );	
	gr.append( ipuppimet.gr_upara_res );	
	makeCanvasGraphs(gr,leg,"upara-res",0,0,120,300);
	gr = [];
	gr.append( pfmet.gr_uperp_scale );
	gr.append( calomet.gr_uperp_scale );
	gr.append( tkmet.gr_uperp_scale );	
	gr.append( ucalomet.gr_uperp_scale );
	gr.append( puppimet.gr_uperp_scale );
	gr.append( tkvtxmet.gr_uperp_scale );	
	gr.append( ipfmet.gr_uperp_scale );	
	gr.append( ipuppimet.gr_uperp_scale );	
	makeCanvasGraphs(gr,leg,"uperp-scale",0,-10,120,10);
	gr = [];
	gr.append( pfmet.gr_uperp_res );
	gr.append( calomet.gr_uperp_res );
	gr.append( tkmet.gr_uperp_res );
	gr.append( ucalomet.gr_uperp_res );	
	gr.append( puppimet.gr_uperp_res );	
	gr.append( tkvtxmet.gr_uperp_res );		
	gr.append( ipfmet.gr_uperp_res );		
	gr.append( ipuppimet.gr_uperp_res );		
	makeCanvasGraphs(gr,leg,"uperp-res",0,0,120,300);


# ##---------------------------------------
# def makeMETData(tag="", t):

# 	for i in range(t.GetEntries()):
# 		t.GetEntry();

# 		if t.m_z < 60 and t.m_z > 120: continue;
# 		ptindex = getindex(t.pt_z,ptranges);
# 		uparaname = tag+"u1";
# 		uperpname = tag+"u2";
# 		upara_ptbin[ptindex].Fill( getattr(t,uparaname) - t.pt_z );
# 		uperp_ptbin[ptindex].Fill( getattr(t,uperpname) );

# def getindex(val,listr):

# 	theindex = -1;
# 	for i,l in range(len(listr)-1):
# 		if val > listr[i] and val < listr[i+1]: 
# 			theindex = i; break;
# 	if val > listr[len(listr)]: theindex = len(listr);
# 	if theindex == -1: theindex = 0;
# 	return theindex;

##---------------------------------------
def makeCanvasGraphs(grs,legs,name,xlo,ylo,xhi,yhi,setlog=False):

	color = [1,2,4,6,7,3,5,8]
	c = ROOT.TCanvas("c_"+name,"c_"+name,1000,800);
	
	leg = ROOT.TLegend(0.65,0.7,0.9,0.9);
	leg.SetBorderSize(0);
	leg.SetFillStyle(0);
	leg.SetTextSize(0.035);
	for i in range(len(legs)):
		leg.AddEntry(grs[i],legs[i],"l");

	hrl = c.DrawFrame(xlo,ylo,xhi,yhi);
	hrl.GetYaxis().SetTitle(grs[0].GetYaxis().GetTitle());
	hrl.GetYaxis().SetTitleOffset(0.85);	
	hrl.GetXaxis().SetTitle(grs[0].GetXaxis().GetTitle());
	
	for i,gr in enumerate(grs):
		gr.SetLineColor(color[i]);
		gr.SetMarkerColor(color[i]);
		gr.SetLineWidth(2);
		gr.Draw("pe");
	leg.Draw();

	if setlog: ROOT.gPad.SetLogy();
	c.SaveAs("plots/%s.pdf" % (name));
	c.SaveAs("plots/%s.png" % (name));

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

	leg = ROOT.TLegend(0.65,0.65,0.9,0.9);
	leg.SetBorderSize(0);
	leg.SetFillStyle(0);
	leg.SetTextSize(0.035);
	for i in range(len(hists)):
		leg.AddEntry(hists[i],legs[i],"l");

	c = ROOT.TCanvas("c"+hists[0].GetName(),"c"+hists[0].GetName(),1000,800);

	max = -999;

	for i in range(len(hists)):
		hists[i].SetLineColor(color[i])
		hists[i].SetLineWidth(2)
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
	ROOT.gPad.SetLogy();
	c.SaveAs("plots/"+name+"_log.pdf")
	c.SaveAs("plots/"+name+"_log.png")



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