import sys
import ROOT
from Sample import TruthSample,Group
from array import array
from ErrorFloat import ErrorFloat
from hist import IntegralError

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.LoadMacro("../../../draw/AtlasStyle.C")
ROOT.gROOT.LoadMacro("../../../draw/AtlasLabels.C")
ROOT.SetAtlasStyle()
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetFrameBorderSize(2)
ROOT.gStyle.SetMarkerSize(0.4)
ROOT.gStyle.SetOptFit(1111)
ROOT.gStyle.SetStatTextColor(1)
ROOT.gStyle.SetStatBorderSize(0)
ROOT.gStyle.SetStatH(0.12)
ROOT.gStyle.SetStatX(.99); ROOT.gStyle.SetStatY(.99)
ROOT.gStyle.SetHatchesLineWidth(3)

hf = ROOT.TFile('../output/PythiaReco.root')

c = ROOT.TCanvas()
c.UseCurrentStyle()
l = ROOT.TLegend(0.48,0.72,0.92,0.92)
l.SetTextSize(0.035)
l.SetBorderSize(0)
l.SetFillColor(0)
#l.SetNColumns(1)


hand = "LH"
cname = 'RecoPlot_%s_overlay' % hand

allp = hf.Get('PythReco_upsilon_all_0_%s' % hand)
p1n = hf.Get('PythReco_upsilon_all_4_%s' % hand)
p0n = hf.Get('PythReco_upsilon_all_3_%s' % hand)

nallp = allp.Integral()

pXn = allp.Clone('pXn')
pXn.Add(p1n,-1)
pXn.Add(p0n,-1)

for h in [allp,p1n,p0n,pXn]:
	print h.Integral()
	h.Rebin(2)
	h.Scale( 1./ nallp)
	h.SetLineWidth(2)
	if hand == 'LH': 
		h.SetLineColor(2)
		h.SetMarkerColor(2)
		h.SetFillColor(2)
	elif hand == 'RH': 
		h.SetLineColor(4)
		h.SetMarkerColor(4)
		h.SetFillColor(4)
	print h.Integral()



allp.SetXTitle("#Upsilon")
#allp.SetFillStyle(3004)
allp.SetFillColor(0)
p1n.SetFillStyle(3004)
p0n.SetFillStyle(3001)
#pXn.SetFillStyle(3005)
pXn.SetLineStyle(3)
pXn.SetFillColor(0)

c.Update()

allp.Draw()
for h in [p0n,p1n,pXn]:#,allp]:
	h.Draw("SAME,HIST")

l.AddEntry(allp, '%s Truth-matched, numTrack=1' % hand, "P")
l.AddEntry(p0n,"%s #tau#rightarrow#pi#nu" % hand )
l.AddEntry(p1n,"%s #tau#rightarrow#pi#pi^{0}#nu" % hand)
l.AddEntry(pXn,"%s other decay modes" % hand)
l.Draw()


c.SaveAs("../plots/%s.eps" % cname)

