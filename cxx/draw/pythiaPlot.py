import sys
import ROOT
from Sample import TruthSample,Group
from array import array
from ErrorFloat import ErrorFloat
from hist import IntegralError

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.LoadMacro("../../draw/AtlasStyle.C")
ROOT.gROOT.LoadMacro("../../draw/AtlasLabels.C")
ROOT.SetAtlasStyle()
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetFrameBorderSize(2)
ROOT.gStyle.SetMarkerSize(0.4)
ROOT.gStyle.SetOptFit(1111)
ROOT.gStyle.SetStatTextColor(1)
ROOT.gStyle.SetStatBorderSize(0)
ROOT.gStyle.SetStatH(0.12)
ROOT.gStyle.SetStatX(.99); ROOT.gStyle.SetStatY(.99)

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

plotdir = 'plots/TruthStudy/Pythia'

entries = ["Zpt","Zmass","TScostheta","cosa","xfrac","upsilon","taupt","vistaupt","taueta","vistaueta","CosDirect"]#,"x1","x2","scale"]
xtitles = ["Z p_{T} [GeV]", "m_{Z} [GeV]", "TauSpinner cos#theta*", "Tau RF cos#theta", "x = p_{T}^{vis}/p_{T}^{#tau}",
            "#Upsilon","#tau p_{T}","#tau p_{T}^{vis}",'#tau #eta','#tau_{vis} #eta',"Direct cos#theta*"]#, "x1", "x2", "PDF scale"]

plots = dict(zip(entries,xtitles))

hf = ROOT.TFile('PythiaTruth.root')
hA = ROOT.TFile('TruthStudy.root')
flavor = 'all'
def fitXfrac( entry, xtitle, mode, hand = ""):
    cname = "Pythia_%s" % entry
    if mode not in ["0"]: cname+= "_%s" % mode
    if hand: cname += "_%s" % hand
    c = ROOT.TCanvas()

    l = ROOT.TLegend(0.77,0.77,0.92,0.92)
    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)

    h = hf.Get("%s_%s_%s_%s_%s" % ('Pythia8', entry, flavor,mode, "LH")).Clone()
    nl = ErrorFloat(*IntegralError(h,of=False))
    hR =hf.Get("%s_%s_%s_%s_%s" % ('Pythia8', entry, flavor, mode, "RH")).Clone()
    nr = ErrorFloat(*IntegralError(hR,of=False))
    h.Add(hR)

    h.Rebin(5)
    h.GetXaxis().SetRange(1,19)
    h.Scale( 1/ h.Integral())
    h.GetXaxis().SetRange( 3,19 )
    ptau = (nr - nl )/ (nr + nl)
    print ptau
    h.Fit('pol1')
    fit = h.GetFunction('pol1')
    p0 = fit.GetParameter(0)
    p1 = fit.GetParameter(1)
    p0e = fit.GetParError(0)
    p1e = fit.GetParError(1)

    p0 = ErrorFloat(  p0,  p0e  )
    p1 = ErrorFloat( p1, p1e )
    two = ErrorFloat( 2, 0)
    print ( p1 ) / ( p1 + two * p0)


    h.SetXTitle(xtitle)
    h.Draw()
    l.AddEntry(h,'Pythia8 Z#tau#tau')
    l.Draw()
    c.SaveAs('%s/%s.pdf' % (plotdir,cname))

def drawKinematics( entry, xtitle, mode, hand = "" ):

    cname = "Pythia_%s" % entry 
    if mode not in ["0"]: cname+= "_%s" % mode
    if hand: cname += "_%s" % hand
    c = ROOT.TCanvas()

    l = ROOT.TLegend(0.77,0.77,0.92,0.92)
    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)

    h = hf.Get("%s_%s_%s_%s_%s" % ('Pythia8', entry, flavor,mode, "LH")).Clone()
    hL = hf.Get("%s_%s_%s_%s_%s" % ('Pythia8', entry, flavor,mode, "LH")).Clone()
    nl = ErrorFloat(*IntegralError(h,of=False))
    hR =hf.Get("%s_%s_%s_%s_%s" % ('Pythia8', entry, flavor, mode, "RH")).Clone()
    nr = ErrorFloat(*IntegralError(hR,of=False))
    h.Add(hR)

    h.Rebin(2)
    hL.Rebin(2);hR.Rebin(2)
    h.Scale( 1/ h.Integral())
    hL.Scale( 1/ hL.Integral() )
    hR.Scale( 1/ hR.Integral() )
    hL.SetMarkerColor( 2 ); hL.SetLineColor(2)
    hR.SetMarkerColor( 4 ); hR.SetLineColor(4)

    maxi = max( h.GetMaximum(), hL.GetMaximum(), hR.GetMaximum())
    h.SetMaximum(1.4*maxi)
    h.SetXTitle(xtitle)
    h.Draw()
    l.AddEntry(h,'Pythia8 Z#tau#tau')
    if 'L' in hand:
        hL.Draw("SAME")
        l.AddEntry( hL, "Left")
    if 'R' in hand:
        hR.Draw("SAME")
        l.AddEntry( hR, "Right")
    l.Draw()
    c.SaveAs('%s/%s.gif+1' % (plotdir,cname))

def drawPolarization( entry, xtitle, mode ):
    cname = "Pythia_Pol" 
    if mode not in ["0"]: cname+= "_%s" % mode

    c = ROOT.TCanvas()

    l = ROOT.TLegend(0.77,0.77,0.92,0.92)
    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)

    h = hf.Get("%s_%s_%s_%s_%s" % ('Pythia8', entry, flavor,mode, "RH")).Clone()
    hL = hf.Get("%s_%s_%s_%s_%s" % ('Pythia8', entry, flavor,mode, "LH")).Clone()
    nl = ErrorFloat(*IntegralError(h,of=False))
    hR =hf.Get("%s_%s_%s_%s_%s" % ('Pythia8', entry, flavor, mode, "RH")).Clone()
    nr = ErrorFloat(*IntegralError(hR,of=False))


    axis = None
    if entry in ['Zmass','scale']:
        axis = [60,65,70]; axis.extend( [72 + 2*i for i in xrange( (80-70)/2 -1 )])
        axis.extend( [80 + i for i in xrange(19)] )
        axis.extend(  [100 + 2*i for i in xrange (4) ])
        axis.extend([110,115,120,125])
        axis.extend([130+i*20 for i in xrange(3)])
    if axis:
        h = h.Rebin(len(axis) -1 , "rebin_h", array('d',axis))
        hL = hL.Rebin(len(axis) -1 , "rebin_hL", array('d',axis))
        hR = hR.Rebin(len(axis) -1 , "rebin_hR", array('d',axis))
    else:
        h.Rebin(2);hL.Rebin(2);hR.Rebin(2)

    hR.Add(hL)
    h.Add( hL, -1)
    h.Divide( hR )
    h.SetMaximum(1.4*h.GetMaximum())
    h.SetXTitle(xtitle)
    h.Draw()
    l.AddEntry(h,'Pythia8 Z#tau#tau')
    l.Draw()
    c.SaveAs('%s/%s.gif+1' % (plotdir,cname))

def comparePolarization( entry, xtitle, mode, diff = True):
    cname = "Pythia_Alp_Pol" 
    if mode not in ["0"]: cname+= "_%s" % mode
    if diff: cname+='_diff'
    c = ROOT.TCanvas()

    l = ROOT.TLegend(0.77,0.77,0.92,0.92)
    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)

    h = hf.Get("%s_%s_%s_%s_%s" % ('Pythia8', entry, flavor,mode, "RH")).Clone()
    hL = hf.Get("%s_%s_%s_%s_%s" % ('Pythia8', entry, flavor,mode, "LH")).Clone()
    nl = ErrorFloat(*IntegralError(h,of=False))
    hR =hf.Get("%s_%s_%s_%s_%s" % ('Pythia8', entry, flavor, mode, "RH")).Clone()
    nr = ErrorFloat(*IntegralError(hR,of=False))


    axis = None
    if entry in ['Zmass','scale']:
        axis = [60,65,70]; axis.extend( [72 + 2*i for i in xrange( (80-70)/2 -1 )])
        axis.extend( [80 + i for i in xrange(19)] )
        axis.extend(  [100 + 2*i for i in xrange (4) ])
        axis.extend([110,115,120,125])
        axis.extend([130+i*20 for i in xrange(3)])
    if axis:
        h = h.Rebin(len(axis) -1 , "rebin_h", array('d',axis))
        hL = hL.Rebin(len(axis) -1 , "rebin_hL", array('d',axis))
        hR = hR.Rebin(len(axis) -1 , "rebin_hR", array('d',axis))
    else:
        h.Rebin(2);hL.Rebin(2);hR.Rebin(2)

    hR.Add(hL)
    h.Add( hL, -1)
    h.Divide( hR )
    h.SetMaximum(1.4*h.GetMaximum())
    h.SetXTitle(xtitle)


    for i,sample in enumerate(group):
        hLA = hA.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "LH")).Clone()
        hRA = hA.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "RH")).Clone()

        if axis:
            hLA = hLA.Rebin(len(axis) -1 , "rebin_hLA", array('d',axis))
            hRA = hRA.Rebin(len(axis) -1 , "rebin_hRA", array('d',axis))
        else:
            hLA.Rebin( 2 )
            hRA.Rebin( 2 )
        nT = hLA.Integral() + hRA.Integral()
        hLA.Scale( sample.xsec / (nT * XSEC) )
        hRA.Scale( sample.xsec / (nT * XSEC)  )
        if i==0: 
            hLeft = hLA.Clone("hLeft")
            hRight = hRA.Clone("hRight")
        else:
            hLeft.Add(hLA); hRight.Add(hRA)

    hT = hRight.Clone("hTotal")
    hT.Add( hLeft )
    hRight.Add( hLeft, -1)
    hRight.Divide(hT)

    h.SetLineColor( 48 ); h.SetMarkerColor( 48)
    hRight.SetLineColor( 38 ); hRight.SetMarkerColor( 38 )

    if diff:
        h.SetYTitle("Pythia P_{#tau} - Alpgen P_{#tau}")
        h.Add( hRight, -1)
        h.Draw()
    else:
        h.SetYTitle("Polarization")

        l.AddEntry(h, 'Pythia8')
        l.AddEntry(hRight,'Alpgen')
        h.Draw()
        hRight.Draw("SAME")
        l.Draw()
    c.SaveAs('%s/%s.gif+1' % (plotdir,cname))

###################################################################################################################
for entry,xtitle in zip(entries,xtitles):
    mode = '0'
    #hand = 'LR'
    #if not entry in ['CosDirect']: continue
    #drawKinematics( entry, xtitle, mode, hand=hand)
    #drawPolarization( entry, xtitle, mode)
    comparePolarization( entry, xtitle, mode)