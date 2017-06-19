import sys
import ROOT
from array import array
import math
import os

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

##---------------------------------------
class METContainer:

	def __init__(self, tree, tag, xtag, fitrange = 300, L1CutVal = -999):
		
		self._tt = tree;
		self._tag = tag; 
		
		if tag == "": self._xtag = 'pf_'+xtag; 
		else: self._xtag = tag + '_' + xtag;

		self._fitrange = fitrange;
		self._L1CutVal = L1CutVal;

		# self._nptbins = 11;
		# self._ptranges = [0,10,20,30,40,50,60,80,100,125,150,200]
		# self._nptbins = 5;
		# self._ptranges = [0,20,50,80,125,200]
		self._nptbins = 4;
		self._ptranges = [0,20,45,70,150]
		self.hl_uperp_ptbin = [];
		self.hl_upara_ptbin = [];
		self.hl_uparaMinusPt_ptbin = [];
		self.hl_pt_ptbin = [];
		for i in range(self._nptbins):
			self.hl_upara_ptbin.append( ROOT.TH1F("h_upara_ptbin_"+str(i)+self._xtag, ";U_{#parallel};N", 400, -400, 400) );
			self.hl_uparaMinusPt_ptbin.append( ROOT.TH1F("hl_uparaMinusPt_ptbin"+str(i)+self._xtag, ";U_{#parallel};N", 400, -400, 400) );
			self.hl_uperp_ptbin.append( ROOT.TH1F("h_uperp_ptbin_"+str(i)+self._xtag, ";U_{#perp};N", 400, -400, 400) );
			self.hl_pt_ptbin.append( ROOT.TH1F("hl_pt_ptbin"+str(i)+self._xtag, ";U_{#parallel};N", 100, 0, 200) );

		self.h_genMet       = ROOT.TH1F("h_genMet"+self._xtag,";MET (GeV);",100,0,1000);	
		self.h_met          = ROOT.TH1F("h_met"+self._xtag,";MET (GeV);",1000,0,1000);
		self.h_met_sig050   = ROOT.TH1F("h_sigmet_050"+self._xtag,";MET (GeV);",1000,0,1000);
		self.h_met_sig100   = ROOT.TH1F("h_sigmet_100"+self._xtag,";MET (GeV);",1000,0,1000);
		self.h_met_sig150   = ROOT.TH1F("h_sigmet_150"+self._xtag,";MET (GeV);",1000,0,1000);
		self.h_met_sig200   = ROOT.TH1F("h_sigmet_200"+self._xtag,";MET (GeV);",1000,0,1000);
		self.h_upara        = ROOT.TH1F("h_upara"+self._xtag,";U_{#parallel};",50,-200,200);
		self.h_uparaMinusPt = ROOT.TH1F("h_uparaMinusPt"+self._xtag,";U_{#parallel} - Z pt;",50,-200,200);
		self.h_uperp        = ROOT.TH1F("h_uperp"+self._xtag,";U_{#perp};",50,-200,200);
		self.h_mz           = ROOT.TH1F("h_mz"+self._xtag,";Z mass (GeV);",75,0,150 );
		self.h_ptz          = ROOT.TH1F("h_ptz"+self._xtag,";Z pt (GeV);",50,0,200);

		self.run();

	def run(self):

		ctag = self._tag;
		nent = self._tt.GetEntries();
		for i in range(self._tt.GetEntries()):
			
			self._tt.GetEntry(i);
			if(i % (1 * nent/100) == 0):
				sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done");
				sys.stdout.flush();

			curmz = self._tt.m_z;
			curpt = self._tt.pt_z;

			curgenmet = getattr( self._tt, "genmet" );
			self.h_genMet.Fill( curgenmet );
			if curgenmet >  50: self.h_met_sig050.Fill( getattr( self._tt, ctag+"met" ) );
			if curgenmet > 100: self.h_met_sig100.Fill( getattr( self._tt, ctag+"met" ) );
			if curgenmet > 150: self.h_met_sig150.Fill( getattr( self._tt, ctag+"met" ) );
			if curgenmet > 200: self.h_met_sig200.Fill( getattr( self._tt, ctag+"met" ) );
			self.h_met.Fill( getattr( self._tt, ctag+"met" ) );

			# stuff for recoil
			if curmz < 60 and curmz > 120: continue;
			if curpt == 0: continue;
			
			self.h_mz.Fill( curmz );
			self.h_ptz.Fill( curpt );
			self.h_uparaMinusPt.Fill ( getattr( self._tt, ctag+"u1" ) - curpt );
			self.h_upara.Fill( getattr( self._tt, ctag+"u1" ) );
			self.h_uperp.Fill( getattr( self._tt, ctag+"u2" ) );

			for i in range(self._nptbins):
				if curpt > self._ptranges[i] and curpt < self._ptranges[i+1]:
					
					self.hl_upara_ptbin[i].Fill( getattr( self._tt, ctag+"u1" ) );
					self.hl_uparaMinusPt_ptbin[i].Fill( getattr( self._tt, ctag+"u1" ) - curpt )
					self.hl_uperp_ptbin[i].Fill( getattr( self._tt, ctag+"u2" ) );
					self.hl_pt_ptbin[i].Fill( curpt );

		print "\n";
		
		##### do a met turn-on curve for a given rate! 
		self._L1Met_1percent = -1.; 
		for i in range(self.h_met.GetNbinsX()):
			cureff = self.h_met.Integral(i+1,self.h_met.GetNbinsX()) / self.h_met.Integral();
			# print cureff, self.h_met.GetBinCenter(i+1);
			if cureff < 0.0025: 
				self._L1Met_1percent = self.h_met.GetBinCenter(i+1);
				break;

		xvals = [5, 50,100,150,175,200,225,250,275,300,350,400,450];
		yvals = [];
		self._MetTurnOn = None;
		if self._L1CutVal != -999:
			for xi in range(len(xvals)):
				denom = 0;
				numer = 0;
				for i in range(self._tt.GetEntries()):
					self._tt.GetEntry(i);
					if self._tt.genmet > xvals[xi]: denom+=1;
					if self._tt.genmet > xvals[xi] and getattr( self._tt, ctag+"met" ) > self._L1CutVal: numer+=1;
				yvals.append( float(numer)/float(denom) );
			self._MetTurnOn = makeAGraph(xvals,yvals);

		if self.h_mz.Integral() == 0: return; # if there are no z's

		##### ------------------------------------------------
		##### ------ fitting for resolutions
		##### ------------------------------------------------

		lx = []; lxe = []; lyS = []; lyeS = []; lyR = []; lyeR = [];
		lySperp = []; lyeSperp = []; lyRperp = []; lyeRperp = [];
		for i in range(self._nptbins):

			# get pt bin center
			ptbincenter = self.hl_pt_ptbin[i].GetMean();

			# fit upara
			gausfunc_para = ROOT.TF1("gausfunc_para"+str(i),"gaus",-1.*self._fitrange+ptbincenter,self._fitrange+ptbincenter)
			gausfunc_para.SetParameter(0,20);
			gausfunc_para.SetParameter(1,0);
			gausfunc_para.SetParameter(2,60);
			self.hl_upara_ptbin[i].Fit(gausfunc_para,"L");
			lx.append( ptbincenter );
			lxe.append( self.hl_pt_ptbin[i].GetRMS() );
			# lyS.append( gausfunc_para.GetParameter(1) );
			# lyeS.append( gausfunc_para.GetParError(1) );
			lyS.append( self.hl_upara_ptbin[i].GetMean()/ptbincenter );
			lyeS.append( self.hl_upara_ptbin[i].GetMean()/ptbincenter );

			gausfunc_paraMinusPt = ROOT.TF1("gausfunc_paraMinusPt"+str(i),"gaus",-1.*self._fitrange+ptbincenter,self._fitrange+ptbincenter)
			gausfunc_paraMinusPt.SetParameter(0,20);
			gausfunc_paraMinusPt.SetParameter(1,0);
			gausfunc_paraMinusPt.SetParameter(2,60);
			self.hl_uparaMinusPt_ptbin[i].Fit(gausfunc_paraMinusPt,"L");
			lyR.append( gausfunc_paraMinusPt.GetParameter(2) );
			lyeR.append( gausfunc_paraMinusPt.GetParError(2) );

			gausfunc_perp = ROOT.TF1("gausfunc_perp"+str(i),"gaus",-1.*self._fitrange+ptbincenter,self._fitrange+ptbincenter)
			gausfunc_perp.SetParameter(0,20);
			gausfunc_perp.SetParameter(1,0);
			gausfunc_perp.SetParameter(2,60);
			self.hl_uperp_ptbin[i].Fit(gausfunc_perp,"L");
			lySperp.append( gausfunc_perp.GetParameter(1) );
			lyeSperp.append( gausfunc_perp.GetParError(1) );			
			lyRperp.append( gausfunc_perp.GetParameter(2) );
			lyeRperp.append( gausfunc_perp.GetParError(2) );

			self.makeCanvasTmp(self.hl_upara_ptbin[i])
			self.makeCanvasTmp(self.hl_uparaMinusPt_ptbin[i])
			self.makeCanvasTmp(self.hl_uperp_ptbin[i])

		for i in range(self._nptbins):
			print lyS[i], lyeS[i], lx[i], lxe[i]

		self.gr_upara_scale = self.makeAGraphErrors(lx,lyS,lxe,lyeS);
		self.gr_upara_res   = self.makeAGraphErrors(lx,lyR,lxe,lyeR);
		self.gr_uperp_scale = self.makeAGraphErrors(lx,lySperp,lxe,lyeSperp);
		self.gr_uperp_res   = self.makeAGraphErrors(lx,lyRperp,lxe,lyeRperp);

		self.gr_upara_scale.SetTitle(";qT;<U_{#parallel}>/<q_{T}>");
		self.gr_upara_res.SetTitle(";qT;#sigma (U_{#parallel} - q_{T})");
		self.gr_uperp_scale.SetTitle(";qT;<U_{#perp}>");
		self.gr_uperp_res.SetTitle(";qT;#sigma (U_{#perp})");

	##-----
	def makeCanvasTmp(self,h):
		ctmp = ROOT.TCanvas("ctmp","ctmp",1000,800);
		h.Draw()
		ctmp.SaveAs("fits/%s.pdf" % (h.GetName()))
		#ctmp.SaveAs("fits/%s.pdf" % (h.GetName()))


	def makeAGraph(self,listx,listy,linecolor = 1, linestyle = 1):
		gr = ROOT.TGraph(len(listx),array('d', listx),array('d', listy));
		gr.SetLineColor(linecolor)
		gr.SetLineStyle(linestyle)
		gr.SetLineWidth(2)
		return gr
	
	def makeAGraphErrors(self,listx,listy,listxe,listye,linecolor = 1, linestyle = 1):
		gr = ROOT.TGraphErrors(len(listx),array('d', listx),array('d', listy),array('d', listxe),array('d', listye));
		gr.SetLineColor(linecolor)
		gr.SetLineStyle(linestyle)
		gr.SetLineWidth(2)
		return gr
	##-----