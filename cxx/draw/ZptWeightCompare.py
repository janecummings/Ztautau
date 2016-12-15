import ROOT
import sys
from Sample import TruthSample,Group
from array import array
from ErrorFloat import ErrorFloat
from hist import IntegralError
import math

draw_dir = '/afs/cern.ch/user/c/cummings/analysis/draw'

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.LoadMacro("%s/AtlasStyle.C" % draw_dir)
ROOT.gROOT.LoadMacro("%s/AtlasLabels.C" % draw_dir)
ROOT.SetAtlasStyle()
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetFrameBorderSize(2)
ROOT.gStyle.SetMarkerSize(0.4)
ROOT.gStyle.SetOptFit(1111)
ROOT.gStyle.SetStatTextColor(1)
ROOT.gStyle.SetStatBorderSize(0)
ROOT.gStyle.SetStatH(0.12)
ROOT.gStyle.SetStatX(.99); ROOT.gStyle.SetStatY(.99)

plotdir = 'plots'
TrueNp0 = TruthSample(name='Np0', ls = 'mcevt*Np0', xsec = 718.87*1.18*1.0)#, ls = 'Np0')
TrueNp1 = TruthSample(name='Np1', ls = 'mcevt*Np1', xsec = 175.75*1.18*1.0)
TrueNp2 = TruthSample(name='Np2', ls = 'mcevt*Np2', xsec = 58.856*1.18*1.0)
TrueNp3 = TruthSample(name='Np3', ls = 'mcevt*Np3', xsec = 15.667*1.18*1.0)
TrueNp4 = TruthSample(name='Np4', ls = 'mcevt*Np4', xsec = 4.0121*1.18*1.0)
TrueNp5 = TruthSample(name='Np5', ls = 'mcevt*Np5', xsec = 1.256*1.18*1.0)

group = Group('TrueZtautau', [TrueNp0, TrueNp1, TrueNp2, TrueNp3, TrueNp4, TrueNp5])

XSEC = 0.
for sample in group:
    XSEC+=sample.xsec
rf = ROOT.TFile('../output/ZptReweight2016.root')



def compare_Inclusive( entry, xtitle, mode, hand = "", flavor = "all"):
    cname = 'ZptComparisonA'
    pname = 'Gen_Zpt_%s_%s' % (entry,hand)
    if mode not in ["0"]: cname+="_%s" % mode
    if hand: cname += "_%s" % hand

    c = ROOT.TCanvas()

    if entry in ['x1','x2']:
        c.SetLogx()
    # Split for ratio plot
    c.Divide(1,2)
    c.cd(1).SetPad(0.0,0.2,1.0,1.0)
    c.cd(2).SetPad(0.0,0.0,1.0,0.2)
    c.cd(1)
    if entry in ['Zmass','scale']:
        ROOT.gPad.SetLogy()

    l = ROOT.TLegend(0.62,0.77,0.92,0.92)
    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)

    axis = None
    if entry in ['Zmass','scale']:
        axis = [60,65,70]; axis.extend( [72 + 2*i for i in xrange( (80-70)/2 -1 )])
        axis.extend( [80 + i for i in xrange(19)] )
        axis.extend(  [100 + 2*i for i in xrange (4) ])
        axis.extend([110,115,120,125])

    h1 = None
    h2 = None
    for i,sample in enumerate(group):
        if hand:
            hia = rf.Get("%s_%s_%s_%s_weighted" % (sample.name, entry,  mode, hand)).Clone("hb")
            hib = rf.Get("%s_%s_%s_%s" % (sample.name, entry,  mode, hand)).Clone("ha")
        else:
            #print "%s_%s_%s_%s" % (sample.name, entry,  mode, "LH")
            hia = rf.Get("%s_%s_%s_%s_weighted" % (sample.name, entry,  mode, "LH")).Clone("hia")
            hiRa =rf.Get("%s_%s_%s_%s_weighted" % (sample.name, entry,  mode, "RH")).Clone("hiRa")
            hia.Add(hiRa)
            hib = rf.Get("%s_%s_%s_%s"  % (sample.name, entry,  mode, "LH")).Clone("hib")
            hiRb =rf.Get("%s_%s_%s_%s" % (sample.name, entry,  mode, "RH")).Clone("hiRb")
            hib.Add(hiRb)
        na = hia.GetEntries(); nb = hib.GetEntries()
        if entry in ['Zpt']:
            hia.GetXaxis().SetRangeUser(0,90)
            hib.GetXaxis().SetRangeUser(0,90)
            na = hia.Integral(0,90)
            nb = hib.Integral(0,90)
        elif entry in ['Zmass']:
            hia.GetXaxis().SetRangeUser(60,120)
            hib.GetXaxis().SetRangeUser(60,120)
            na = hia.Integral(60,120)
            nb = hib.Integral(60,120)
        hia.Scale( sample.xsec/(XSEC*na) )#
        hib.Scale( sample.xsec/(XSEC*nb) )#

        if i == 0: 
            h2 = hia
            h1 = hib
        else: 
            h2.Add( hia )
            h1.Add( hib )
    h2.SetLineColor(46)
    h2.SetMarkerColor(46)
    h1.SetXTitle( xtitle )

    maxi = max(h1.GetMaximum(),h2.GetMaximum())
    draw = ''

    h1.SetMaximum( 1.4 * maxi )
    h1.Draw(draw)
    h2.Draw("SAME")
    l.AddEntry(h1,"Alpgen Z#rightarrow#tau#tau")
    l.AddEntry(h2, 'Z p_{T} Reweighted')



    tm = ROOT.TPaveText(0.18,0.84,0.28,0.9,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    tmtext = ''
    if hand == 'LH': tmtext = 'Left-Handed'
    elif hand == 'RH': tmtext = 'Right-Handed'
    tm.AddText( tmtext )
    tm.Draw() 

    l.Draw()

    c.cd(2)
    
    ROOT.gPad.SetTopMargin(0)
    ratio = h2.Clone("ratio")
    ratio.SetLabelSize(.09,"y")
    ratio.Divide(h1)
    ratio.SetMinimum(0.8)

    ratio.Draw()
    xmin = ratio.GetXaxis().GetXmin()
    xmax = ratio.GetXaxis().GetXmax()
    cl = ROOT.TLine(xmin,1,xmax,1)
    cl.SetLineStyle(3)
    cl.Draw("SAME")
    c.SaveAs('%s/%s.eps' % (plotdir,pname))


def compare_polarization( entry, xtitle, mode,):
    cname = 'ZptPolarization'
    if mode not in ["0"]: cname+="_%s" % mode


    c = ROOT.TCanvas()

    if entry in ['x1','x2']:
        c.SetLogx()
    # Split for ratio plot
    c.Divide(1,2)
    c.cd(1).SetPad(0.0,0.2,1.0,1.0)
    c.cd(2).SetPad(0.0,0.0,1.0,0.2)
    c.cd(1)


    l = ROOT.TLegend(0.62,0.77,0.92,0.92)
    if entry == 'Zpt':
        l = ROOT.TLegend(0.52,0.77,0.67,0.92)

    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)

    axis = None
    if entry in ['Zmass','scale']:
        axis = [60,65,70]; axis.extend( [72 + 2*i for i in xrange( (80-70)/2 -1 )])
        axis.extend( [80 + i for i in xrange(19)] )
        axis.extend(  [100 + 2*i for i in xrange (4) ])
        axis.extend([110,115,120,125])

    h1 = None
    hL1 = None
    hR1 = None
    h2 = None
    hL2 = None
    hR2 = None
    for i,sample in enumerate(group):
        hia = rf.Get("%s_%s_%s_%s_weighted" % (sample.name, entry,  mode, "LH")).Clone("hia")
        hiLa = rf.Get("%s_%s_%s_%s_weighted" % (sample.name, entry,  mode, "LH")).Clone("hiLa")
        hiRa =rf.Get("%s_%s_%s_%s_weighted" % (sample.name, entry,  mode, "RH")).Clone("hiRa")

        hia.Add(hiRa)

        hib = rf.Get("%s_%s_%s_%s"  % (sample.name, entry,  mode, "LH")).Clone("hib")
        hiLb =rf.Get("%s_%s_%s_%s" % (sample.name, entry,  mode, "LH")).Clone("hiLb")
        hiRb =rf.Get("%s_%s_%s_%s" % (sample.name, entry,  mode, "RH")).Clone("hiRb")
        hib.Add(hiRb)

        na = hia.Integral(); nb = hib.Integral()
        if entry in ['Zpt']:
            for h in [hia,hiLa,hiRa,hib,hiLb,hiRb]:
                h.Rebin(2)
                h.GetXaxis().SetRangeUser(0,90)
            na = hia.Integral(0,90)
            nb = hib.Integral(0,90)
        elif entry in ['Zmass']:
            for h in [hia,hiLa,hiRa,hib,hiLb,hiRb]:
                #h.Rebin(2)
                h.GetXaxis().SetRangeUser(60,120)
            na = hia.Integral(60,120)
            nb = hib.Integral(60,120)

        for h in [hia,hiLa,hiRa]:
            #h.Rebin(2)
            h.Scale(sample.xsec / (XSEC*na))
        for h in [hib,hiLb,hiRb]:
            #h.Rebin(2)
            h.Scale(sample.xsec / (XSEC*nb))
        if i == 0:
            h1 = hia; h2 = hib
            hL1 = hiLa; hL2 = hiLb
            hR1 = hiRa; hR2 = hiRb
        else:
            h1.Add(hia); h2.Add(hib)
            hL1.Add(hiLa); hL2.Add(hiLb)
            hR1.Add(hiRa); hR2.Add(hiRb)

    if entry == 'Zmass':
        blo = 60; bhi = 120
        nL1 = ErrorFloat(*IntegralError(hL1, of=False,binlo=blo,binhi=bhi))
        nL2 = ErrorFloat(*IntegralError(hL2, of=False,binlo=blo,binhi=bhi))
        nR1 = ErrorFloat(*IntegralError(hR1, of=False,binlo=blo,binhi=bhi))
        nR2 = ErrorFloat(*IntegralError(hR2, of=False,binlo=blo,binhi=bhi))

        pol1 = (nR1-nL1)/(nR1+nL1)
        pol2 = (nR2-nL2)/(nR2+nL2)
        print '%s -> %s' % (pol1,pol2)

    p1 = hR1.Clone('p1')
    p1.Add(hL1,-1)
    p1.Divide(h1)

    p2 = hR2.Clone('p2')
    p2.Add(hL2,-1)
    p2.Divide(h2)
    p2.SetLineColor(46)
    p2.SetMarkerColor(46)
    p1.SetXTitle( xtitle )

    maxi = max(p1.GetMaximum(),p2.GetMaximum())

    draw = ''

    p1.SetYTitle("Polarization")
    p1.Draw(draw)
    p2.Draw("SAME")
    l.AddEntry(p1,"Inclusive Alpgen")
    l.AddEntry(p2,"ZpT Reweighted")

    tm = ROOT.TPaveText(0.18,0.8,0.25,0.85,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    tm.AddText( '' )
    #tm.Draw()
    l.Draw()

    c.cd(2)
    
    ROOT.gPad.SetTopMargin(0)
    diff = p2.Clone("diff")
    diff.SetLabelSize(.1,"y")
    diff.Add(p1,-1)
    diff.Draw()
    xmin = 0
    xmax = 90
    cl = ROOT.TLine(xmin,0,xmax,0)
    cl.SetLineStyle(3)
    cl.Draw("SAME")

    c.SaveAs( '%s/%s.eps' % (plotdir,cname))

######################################################################################################################
entries = {"Zpt": "p_{T}^{Z} [GeV]", "Zmass": "m_{Z} [GeV]",'taupt':"p_{T}^{#tau} [GeV]",
            'vistaupt':"p_{T,vis}^{#tau} [GeV]", "taueta":"#eta^{#tau}", "vistaueta":"#eta_{vis}^{#tau}",
            'ptlep':"p_{T}^{l} [GeV]",'upsilon':"#Upsilon"}
for entry in entries:
    mode = 0
    if entry == 'upsilon':
        #continue
        mode=4
    compare_Inclusive( entry, entries[entry], mode, hand = 'RH' )
    compare_Inclusive( entry, entries[entry], mode, hand = 'LH' )

    

