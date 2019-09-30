import ROOT 

def CreateCanvas(CanvasName = "myPlot", LogY = False, Grid = True):
  c = ROOT.TCanvas(CanvasName,CanvasName,800,800)
  c.SetLeftMargin(0.13)
  if Grid: c.SetGrid()
  if LogY: c.SetLogy()
  return c

def DrawPrelimLabel(c):
  c.cd()
  tex = ROOT.TLatex()
  tex.SetTextSize(0.03)
  tex.DrawLatexNDC(0.13,0.91,"#scale[1.5]{CMS}")
  tex.Draw("same")
  
def DrawLumiLabel(c, Lumi = "35.9"):
  c.cd()
  tex = ROOT.TLatex()
  tex.SetTextSize(0.035)
  tex.SetTextFont(42)
  if Lumi!="":
    toDisplay = Lumi + " fb^{-1} (13 TeV)"
    tex.DrawLatexNDC(0.66,0.91,toDisplay)
  else:
    toDisplay = "(13 TeV)"
    tex.DrawLatexNDC(0.77,0.91,toDisplay)
  tex.Draw("same")

def SaveCanvas(c, PlotName = "myPlotName"):
  c.cd()
  c.SaveAs(PlotName + ".pdf")
  c.SaveAs(PlotName + ".png")
  c.SaveAs(PlotName + ".root")

