import sys,os
import ROOT
import streams
from Drawer import *
import hist
import run
import Groups
from ErrorFloat import ErrorFloat
from printing import row
fcol = 30; col = 20
avg_p = -0.1414

# Atlas Style
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.LoadMacro("../AtlasStyle.C")
ROOT.gROOT.LoadMacro("../AtlasLabels.C")
ROOT.SetAtlasStyle()
ROOT.gStyle.SetFrameBorderSize(1)

stream = 'el'


def compare_signals( entry, hand = '', norm = True):
    rf = ROOT.TFile('../output/%s/%s_regions.root' % (stream + '_Zpt_NoWeight',stream))
    wf = ROOT.TFile('../output/%s/%s_regions.root' % (stream+'_Zpt',stream))
    # Canvas
    cname = '%s_%s' % (stream,entry.name)
    if hand: cname += '_%s' % hand
    if norm: cname += '_%s' % norm
    c = ROOT.TCanvas(cname,cname,1200,1000)
    # Split for ratio plot
    c.Divide(1,2)
    c.cd(1).SetPad(0.0,0.2,1.0,1.0)
    c.cd(2).SetPad(0.0,0.0,1.0,0.2)
    c.cd(2).SetTopMargin(0)
    c.cd(1)
    leg = ROOT.TLegend( 0.65,0.75,0.90,0.92)
    #leg = ROOT.TLegend( 0.22,0.75,0.42,0.85)
    leg.SetTextSize(0.034)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    if not hand:
        hname = 'h_%s_Ztautau_SR_OS' % entry.name
    elif hand:
        hname = 'h_%s_%s_Ztautau_SR_OS' % (entry.name, hand)

    h1 = rf.Get(hname).Clone('h1')
    h2 = wf.Get(hname).Clone('h2')
    print h2.Integral(0,90)
    print h2.Integral(90,h2.GetNbinsX()+1)

    for h in [h1,h2]:
        h.Rebin(entry.rebin)
        if norm: h.Scale(1./h.Integral())
        
    maxi = max(h1.GetMaximum(), h2.GetMaximum())

    h2.SetLineColor(2)
    h2.SetMarkerColor(2)

    if not hand: hh = 'Inclusive'
    else: 
        hh = hand
    leg.AddEntry(h1,'Unweighted %s' % hh)
    leg.AddEntry(h2,'ZpT reweighted')

    h1.SetXTitle(entry.xtitle)
    h1.SetMaximum(1.4*maxi)



    if norm: h1.SetYTitle('Normalized')
    h1.Draw()
    h2.Draw("SAME")
    leg.Draw()

    tm = ROOT.TPaveText(0.17,0.85,0.21,0.90,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    if stream == 'mu':
        tm.AddText('#mu#tau_{had}')
    elif stream == 'el':
        tm.AddText('e#tau_{had}')
    tm.Draw()

    c.cd(2)
    ROOT.gPad.SetTopMargin(0)

    ratio = h2.Clone('ratio')
    ratio.Divide(h1)
    ratio.GetXaxis().SetLabelSize(0.14)
    #ratio.GetYaxis().SetLabelSize(0.16)
    ratio.GetYaxis().SetTitleSize(0.16)
    ratio.SetMarkerSize(0.2)
    ratio.SetYTitle('Weighted/Unweighted')
    err = h1.Clone('err')
    err.Divide(err)
    err.SetFillStyle(3354)
    err.SetLineColor(1)
    err.SetFillColor(12)
    err.SetMarkerStyle(1)

    xmin = ratio.GetXaxis().GetXmin()
    xmax = ratio.GetXaxis().GetXmax()
    cl = ROOT.TLine(xmin,1,xmax,1)

    rmax = 1.5; rmin = 0.5
    if norm: rmax = 1.1; rmin = 0.9

    ratio.SetMaximum(rmax)
    ratio.SetMinimum(rmin)

    ratio.Draw()
    cl.Draw("SAME")
    err.Draw("SAME,E2")



    c.SaveAs('plots/ZptWeight/%s/%s.gif' % (stream,cname))


def compare_scales( ):
    rf = ROOT.TFile('../output/%s/%s_regions.root' % (stream+'_Zpt_NoWeight',stream))
    wf = ROOT.TFile('../output/%s/%s_regions.root' % (stream+'_Zpt',stream))


    tname = 'h_upsilon_Ztautau_SR_OS'
    lname = 'h_upsilon_LH_Ztautau_SR_OS'
    rname = 'h_upsilon_RH_Ztautau_SR_OS'

    for name in [tname,lname,rname]:
        t1 = rf.Get(name).Clone('t1')
        t2 = wf.Get(name).Clone('t2')

        n1 = ErrorFloat(*hist.IntegralError(t1))
        n2 = ErrorFloat(*hist.IntegralError(t2))
        print n1
        print n2
        print n2/n1



compare_scales()
sys.exit()

#for entry in [truZpt,tau_et, tau_leadTrkPt, visMass, tau_eta, lep_pt, lep_eta]:
for entry in [upsilon]:
    #for hand in ['']:#,'LH','RH']:
    for hand in ['','LH','RH']:
        compare_signals(entry, hand = hand, norm = False)
        compare_signals(entry, hand = hand, norm = True)


