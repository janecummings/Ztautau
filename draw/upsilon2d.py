import ROOT
import hist
import run
import os
import streams
from ErrorFloat import ErrorFloat
from Drawer import *
from printing import row
import math 
from array import array

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.LoadMacro("AtlasStyle.C")
ROOT.gROOT.LoadMacro("AtlasLabels.C")
ROOT.SetAtlasStyle()
ROOT.gStyle.SetFrameBorderSize(1)
ROOT.gStyle.SetPalette(55)

pythia = 'TestXX'
modes = {1:'All SP',3:'#tau#rightarrow#pi#nu',4:'#tau#rightarrow#pi#pi^{0}#nu'}


def draw_preselection(mode,hand):
    hf = ROOT.TFile('~/analysis/truth/cxx/thesis/Pythia%s.root' % pythia)
    cname = 'preselect_%s_%s' % (mode,hand)
    c = ROOT.TCanvas(cname,cname,1500,1000)
    c.UseCurrentStyle()
    h = hf.Get('upsilon_upsilon20_%s_%s' % (mode,hand))
    h.SetXTitle('Truth #Upsilon')
    h.SetYTitle('Reconstructed #Upsilon')
    h.GetYaxis().SetTitleOffset(1.1)
    #h.SetContour(40)
    #h.Scale(1/h.Integral())
    
    h.Draw('COLZ')
    tm = ROOT.TPaveText(0.16,0.75,0.36,0.90,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetShadowColor(0)
    tm.SetTextSize(0.045)
    tm.SetTextAlign(10)
    tm.AddText(hand)
    tm.AddText(modes[mode])
    tm.AddText('Correlation: %.3f' % h.GetCorrelationFactor())
    tm.Draw()
    c.SetRightMargin(0.13)
    c.Print('../plots/upsilon2d/%s.eps'%cname)


def draw_ldPt(mu, var = tau_leadTrkPt_20, split = False):


    hf = ROOT.TFile('../output/prod4/leadpt/%s_regions.root' % mu )
    cname = 'upsilon_%s_%s' % (mu,var.name)
    c = ROOT.TCanvas(cname,cname,1500,1000)
    if split:
        c = ROOT.TCanvas(cname,cname,1500,1200)
        c.Divide(1,2)
        c.cd(1).SetPad(0.0,0.9,1.0,1.0)
        c.cd(2).SetPad(0.0,0.0,1.0,0.9)
    c.UseCurrentStyle()
    h = hf.Get('h_upsilon_%s_Ztautau_SR_OS' % var.name)
    h.SetXTitle('#Upsilon')
    h.SetYTitle(var.xtitle)
    h.GetYaxis().SetTitleOffset(1.1)
    #h.SetContour(40)
    #h.Scale(1/h.Integral())
    
    h.Draw('COLZ')
    if split:
        c.cd(1)
        tm = ROOT.TPaveText(0.0,0.0,0.36,0.6,"NDC")
        tm.SetTextSize(0.44)
    else:
        tm = ROOT.TPaveText(0.0,0.80,0.36,0.90,"NDC")
        tm.SetTextSize(0.044)
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetShadowColor(0)
    tm.SetTextAlign(10)
    if mu == 'mu':
        tm.AddText('#mu#tau_{had}')
    elif mu == 'el':
        tm.AddText('e#tau_{had}')
    tm.AddText('Correlation: %.2g' % h.GetCorrelationFactor())
    tm.Draw()

    c.SetRightMargin(0.13)
    c.Print('../plots/upsilon2d/%s.eps'%cname)


def draw_select(channel,hand):
    hf = ROOT.TFile('~/panalysis/output/%s_upsilon2d/%s_regions.root' % (channel,channel))
    cname = 'select_%s_%s' % (channel,hand)
    c = ROOT.TCanvas(cname,cname,1500,1000)
    c.UseCurrentStyle()

    h = hf.Get('h_upsi_upsilon_%s_Ztautau_SR_OS' % (hand))
    h.SetXTitle('Truth #Upsilon')
    h.SetYTitle('Reconstructed #Upsilon')
    h.GetYaxis().SetTitleOffset(1.1)
    h.RebinX(2)
    h.RebinY(2)
    h.SetContour(61)
    #h.Scale(1/h.Integral())
    print h.GetCorrelationFactor()
    h.Draw('COLZ')
    tm = ROOT.TPaveText(0.16,0.78,0.34,0.90,"NDC")
    tm.SetShadowColor(0)
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.045)
    tm.SetTextAlign(10)
    #tm.AddText(hand)
    if 'mu' in channel: ch = '#mu'
    elif 'el' in channel: ch = 'e'
    tm.AddText('%s %s#tau_{had} SR' % (hand,ch))
    tm.AddText('Correlation: %.3f' % h.GetCorrelationFactor())
    tm.Draw()
    c.SetRightMargin(0.13)
    c.Update()
    c.Print('../plots/upsilon2d/%s.eps'%cname)


def draw_upsilon(mu):
    rf = ROOT.TFile('~/panalysis/output/%s_thesis/fit/%s_regions.root' % (mu.name,mu.name))   
    entry = upsilon_20

    cname = 'upsilon_%s' % (mu.name)
    c = ROOT.TCanvas(cname,cname,1500,1000)
    c.UseCurrentStyle()
    l = ROOT.TLegend(.7,.7,.9,.9)
    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    total = rf.Get('h_upsilon_Ztautau_SR_OS')
    lh = rf.Get('h_upsilon_LH_Ztautau_SR_OS')
    rh = rf.Get('h_upsilon_RH_Ztautau_SR_OS')

    chi2 = ROOT.Double(0)
    ndf = ROOT.Long(0)
    igoot= ROOT.Long(0)
    print lh.Chi2TestX(rh,chi2,ndf,igoot,"WW")
    print chi2
    print ndf
    print igoot
    #print rh.Chi2Test(lh,"WW")

    lh.SetLineColor(mu.LH.fill)
    lh.SetMarkerColor(mu.LH.fill)
    rh.SetLineColor(mu.RH.fill)
    rh.SetMarkerColor(mu.RH.fill)

    l.AddEntry(total, 'Z#rightarrow#tau#tau')
    l.AddEntry(lh, 'Left-Handed')
    l.AddEntry(rh,'Right-Handed')

    maxi = max(total.GetMaximum(),lh.GetMaximum(),rh.GetMaximum())
    total.SetMaximum(1.3*maxi)
    total.Draw("PE")
    lh.Draw("SAME,PE")
    rh.Draw("SAME,PE")

    l.Draw()



    tm = ROOT.TPaveText(0.19,0.78,0.22,0.90,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.05)
    tm.SetTextAlign(10)
    if mu.name == 'mu':
        tm.AddText('#mu#tau_{had}')
    elif mu.name == 'el':
        tm.AddText('e#tau_{had}')
    tm.AddText('#chi^{2}/n.d.f. = %i/%i' % (chi2,ndf))
    tm.Draw()

    c.Print('../plots/note/%s.eps' % cname)


def draw2d():
    # draw_preselection(1,'LH')
    # draw_preselection(1,'RH')
    # draw_preselection(4,'LH')
    # draw_preselection(4,'RH')
    #############
    draw_select('mu','LH')
    draw_select('mu','RH')
    draw_select('el','LH')
    draw_select('el','RH')

##########################################################################################################
draw2d()
#draw_ldPt('el', met )
#draw_ldPt('el', met, split = True )
#draw_upsilon(streams.mu)
#draw_upsilon(streams.el)