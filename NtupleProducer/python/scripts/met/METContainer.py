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
		
		self._xtag = tag
		#if tag == "": self._xtag = 'pf_'+xtag; 
		#else: self._xtag = tag + '_' + xtag;

		self._fitrange = fitrange;
		self._L1CutVal = L1CutVal;

		self.h_genMet       = ROOT.TH1F("h_genMet"+self._xtag    ,";%sMET (GeV);"%self._xtag,100,0,1000);	
		self.h_met          = ROOT.TH1F("h_met"+self._xtag       ,";%sMET (GeV);"%self._xtag,1000,0,1000);
		self.h_met_sig050   = ROOT.TH1F("h_sigmet_050"+self._xtag,";%sMET (GeV);"%self._xtag,1000,0,1000);
		self.h_met_sig100   = ROOT.TH1F("h_sigmet_100"+self._xtag,";%sMET (GeV);"%self._xtag,1000,0,1000);
		self.h_met_sig150   = ROOT.TH1F("h_sigmet_150"+self._xtag,";%sMET (GeV);"%self._xtag,1000,0,1000);
		self.h_met_sig200   = ROOT.TH1F("h_sigmet_200"+self._xtag,";%sMET (GeV);"%self._xtag,1000,0,1000);

		self.run();

	def run(self):

		ctag = self._tag;
		nent = self._tt.GetEntries();
		for i in range(self._tt.GetEntries()):
			
			self._tt.GetEntry(i);
			if(i % (1 * nent/100) == 0):
				sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done");
				sys.stdout.flush();

			curgenmet = getattr( self._tt, "METGen" );
			self.h_genMet.Fill( curgenmet );
			if curgenmet >  50: self.h_met_sig050.Fill( getattr( self._tt, ctag ) );
			if curgenmet > 100: self.h_met_sig100.Fill( getattr( self._tt, ctag ) );
			if curgenmet > 150: self.h_met_sig150.Fill( getattr( self._tt, ctag ) );
			if curgenmet > 200: self.h_met_sig200.Fill( getattr( self._tt, ctag ) );
			self.h_met.Fill( getattr( self._tt, ctag ) );

		print "\n";
		
		##### do a met turn-on curve for a given rate! 
		self._L1Met_1percent = -1.; 
		for i in range(self.h_met.GetNbinsX()):
			cureff = self.h_met.Integral(i+1,self.h_met.GetNbinsX()) / self.h_met.Integral();
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
					if self._tt.METGen > xvals[xi]: denom+=1;
					if self._tt.METGen > xvals[xi] and getattr( self._tt, ctag ) > self._L1CutVal: numer+=1;
				yvals.append( float(numer)/float(denom) );
			self._MetTurnOn = makeAGraph(xvals,yvals);

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
