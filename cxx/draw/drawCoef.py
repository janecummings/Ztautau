import ROOT
import sys,os
from Sample import TruthSample,Group
from array import array
import math
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

AlpgenOn = True
PythiaOn = False

TrueNp0 = TruthSample(name='Np0', ls = 'mcEvt*Np0', xsec = 718.87*1.18*1.0)#, ls = 'Np0')
TrueNp1 = TruthSample(name='Np1', ls = 'mcEvt*Np1', xsec = 175.75*1.18*1.0)
TrueNp2 = TruthSample(name='Np2', ls = 'mcEvt*Np2', xsec = 58.856*1.18*1.0)
TrueNp3 = TruthSample(name='Np3', ls = 'mcEvt*Np3', xsec = 15.667*1.18*1.0)
TrueNp4 = TruthSample(name='Np4', ls = 'mcEvt*Np4', xsec = 4.0121*1.18*1.0)
TrueNp5 = TruthSample(name='Np5', ls = 'mcEvt*Np5', xsec = 1.256*1.18*1.0)

group1 = Group('Alpgen', [TrueNp0, TrueNp1, TrueNp2, TrueNp3, TrueNp4, TrueNp5])
Pythia = TruthSample(name = 'Pythia')
group2 = Group('Pythia',[Pythia])
group3 = Group('AlpgenZpTReweighted', [TrueNp0, TrueNp1, TrueNp2, TrueNp3, TrueNp4, TrueNp5])
group4 = Group('AlpgenZmassReweighted', [TrueNp0, TrueNp1, TrueNp2, TrueNp3, TrueNp4, TrueNp5])

pf = ROOT.TFile('CSPythia.root')#,'UPDATE')
af = ROOT.TFile('CSAlpgen.root')#,'UPDATE')
afrw = ROOT.TFile('CSAlpgen_Reweight_Fix3.root')#,'UPDATE')
if (AlpgenOn ):
    group = group1
    hf = af
    XSEC = 0.
    for sample in group:
        XSEC+=sample.xsec
if (PythiaOn):
    group = group2
    hf = pf

xtitles = {'Zmass':'m_{Z} [GeV]', 'ZpT':'p_{T}^{Z} [GeV]', 'taupt': 'p_{T}^{#tau} [GeV]',
            'leppt':'p_{T}^{lep}', 'CosTheta':'cos#theta', 'ZpT_AZ':'p_{T}^{Z} [GeV]','CosTheta_AZ':'cos#theta' }


def get_coef( entry , coef, sample = False):
    c = ROOT.TCanvas()
    l = ROOT.TLegend(0.77,0.77,0.9,0.92)
    l.SetTextSize(0.034)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    hists = []
    maxi = 0.05
    mini = 0
    #flavor = 'all'
    #flavor = 'A0'
    if 'ZpT' in entry:
        axis = [0]
        axis.extend([10*i for i in xrange(1,6)])
        axis.extend([60+20*i for i in xrange(3)])
        axis.extend([120,200])
    elif 'Zmass' in entry:
        axis = [60,65,70]; axis.extend( [72 + 2*i for i in xrange( (80-70)/2 -1 )])
        axis.extend( [80 + i for i in xrange(19)] )
        axis.extend(  [100 + 2*i for i in xrange (4) ])
        axis.extend([110,115,120,125])
        axis.extend([130+i*20 for i in xrange(3)])

    hf.cd()
    w = None
    wA = None
    for i,sample in enumerate(group):
        hA = hf.Get("%s_%s_%s" % (sample.name, "A", entry)).Clone("hA")
        h = hf.Get("%s_%s_%s" % (sample.name, coef, entry)).Clone("h")
        hA = hA.Rebin(len(axis)-1, "rebin_hA", array('d',axis))
        h = h.Rebin(len(axis)-1, "rebin_h", array('d',axis))
        NORM = hA.Integral(0,hA.GetNbinsX()+1)
        #NORM = 1.
        for hh in [hA,h]:
            hh.SetBinContent( hh.GetNbinsX(), hh.GetBinContent(hh.GetNbinsX()) +hh.GetBinContent(hh.GetNbinsX()+1))
            hh.SetBinError( hh.GetNbinsX(), math.sqrt( pow(hh.GetBinError(hh.GetNbinsX()),2) + pow(hh.GetBinError(hh.GetNbinsX()+1),2)))
            hh.SetBinContent( hh.GetNbinsX()+1,0)
            hh.SetBinError(hh.GetNbinsX()+1,0)
            hh.Scale(sample.xsec/(XSEC*NORM))
        if not w: 
            w = h
            wA = hA
        else:
            w.Add(h)
            wA.Add(hA)

    if 'ZpT' in entry:
        w.SetXTitle("Z p_{T} [GeV]")
    elif 'Zmass' in entry:
        w.SetXTitle("m_{Z} [GeV]")
    if not entry == 'CosTheta':
        w.Divide( wA )
        w.Scale(4)
    for x in xrange( h.GetNbinsX()+1):
        print w.GetBinLowEdge(x), w.GetBinContent(x), w.GetBinError(x)
    print w.GetBinContent(1)

    w.SetName('%s_%s' % (sample.name,coef))
    l.AddEntry(w,group.name)
    hists.append(w)
    maxi = max(maxi,w.GetMaximum())
    mini = min(mini,w.GetMinimum())

    draw = ''

    w.SetYTitle(coef)

    hists.reverse()
    for i,h in enumerate(hists):
        h.SetMinimum(1.7*mini)
        h.SetMaximum( 1.7 * maxi )
        h.Draw(draw)
        if i==0: draw+='SAME'
    tm = ROOT.TPaveText(0.22,0.8,0.48,0.85,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    if 'AZ' in entry:
        tt = '80 GeV < m_{Z} < 100 GeV'
    else: tt = ''
    tm.AddText( tt )
    tm.Draw()
    l.Draw()
    c.SaveAs( 'plots/%s_%sCoef_%s_fix.pdf' % (group.name, coef,entry))

    wf = ROOT.TFile('PythiaWeights.root', "UPDATE")
    wf.cd()
    h.SetName('%s_%s_%s_fix' % (group.name, coef, entry))
    h.Write()
    wf.Write()
    wf.Close()


def compare_coef( entry, coef):
    pname = 'Coef'
    c = ROOT.TCanvas()
    l = ROOT.TLegend(0.77,0.77,0.9,0.92)
    #l.SetNColumns(2)
    l.SetTextSize(0.034)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    hists = []
    maxi = 0.05
    mini = 0
    if 'ZpT' in entry:
        axis = [0]
        axis.extend([10*i for i in xrange(1,6)])
        axis.extend([60+20*i for i in xrange(3)])
        axis.extend([120,200])
    elif 'Zmass' in entry:
        axis = [60,65,70]; axis.extend( [72 + 2*i for i in xrange( (80-70)/2 -1 )])
        axis.extend( [80 + i for i in xrange(19)] )
        axis.extend(  [100 + 2*i for i in xrange (4) ])
        axis.extend([110,115,120,125])
        axis.extend([130+i*20 for i in xrange(3)])

    af = afrw
    af.cd()
    w = None
    wA = None
    for i,sample in enumerate(group1):# + group1[2:7]):
        wt = ''
        if af == afrw: wt = '_rwZpt'
        hA = af.Get("%s_%s_%s%s" % (sample.name, "A", entry,wt)).Clone("hA")
        h = af.Get("%s_%s_%s%s" % (sample.name, coef, entry,wt)).Clone("h")
        print hA.Integral()
        hA = hA.Rebin(len(axis)-1, "rebin_hA", array('d',axis))
        h = h.Rebin(len(axis)-1, "rebin_h", array('d',axis))
        NORM = hA.Integral(0,hA.GetNbinsX()+1)
        for hh in [hA,h]:
            hh.SetBinContent( hh.GetNbinsX(), hh.GetBinContent(hh.GetNbinsX()) +hh.GetBinContent(hh.GetNbinsX()+1))
            hh.SetBinError( hh.GetNbinsX(), math.sqrt( pow(hh.GetBinError(hh.GetNbinsX()),2) + pow(hh.GetBinError(hh.GetNbinsX()+1),2)))
            hh.SetBinContent( hh.GetNbinsX()+1,0)
            hh.SetBinError(hh.GetNbinsX()+1,0)
            hh.Scale(sample.xsec/(XSEC*NORM))
        if not w: 
            w = h
            wA = hA
        else:
            w.Add(h)
            wA.Add(hA)
    w.Divide(wA)
    w.Scale(4)

    pf.cd()
    for i,sample in enumerate(group2):
        hA = pf.Get("%s_%s_%s" % (sample.name, "A", entry)).Clone("hAPyth")
        h = pf.Get("%s_%s_%s" % (sample.name, coef, entry)).Clone("hPyth")
        hA = hA.Rebin(len(axis)-1, "rebin_hA_pyth", array('d',axis))
        h = h.Rebin(len(axis)-1, "rebin_h_pyth", array('d',axis))
        for hh in [hA,h]:
            hh.SetBinContent( hh.GetNbinsX(), hh.GetBinContent(hh.GetNbinsX()) +hh.GetBinContent(hh.GetNbinsX()+1))
            hh.SetBinError( hh.GetNbinsX(), math.sqrt( pow(hh.GetBinError(hh.GetNbinsX()),2) + pow(hh.GetBinError(hh.GetNbinsX()+1),2)))
            hh.SetBinContent( hh.GetNbinsX()+1,0)
            hh.SetBinError(hh.GetNbinsX()+1,0)
            #hh.Scale(sample.xsec/XSEC)
    h.Divide(hA)
    h.Scale(4)
    if 'ZpT' in entry:
        w.SetXTitle("Z p_{T} [GeV]")
    elif 'Zmass' in entry:
        w.SetXTitle("m_{Z} [GeV]")

    h.SetLineColor(4)
    h.SetMarkerColor(4)
    l.AddEntry( w, 'Alpgen' )
    l.AddEntry( h, 'Pythia')
    #hists.append(w)
    maxi = max(maxi,w.GetMaximum(),h.GetMaximum())
    mini = min(mini,w.GetMinimum(),h.GetMinimum())
    print mini

    draw = ''

    w.SetYTitle(coef)
    hists = [w,h]
    for i,h in enumerate(hists):
        h.SetMinimum(1.7*mini)
        h.SetMaximum( 1.7 * maxi )
        h.Draw(draw)
        if i==0: draw+='SAME'
    tm = ROOT.TPaveText(0.22,0.8,0.48,0.85,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    if 'AZ' in entry:
        tt = '80 GeV < m_{Z} < 100 GeV'
    else: tt = ''
    tm.AddText( tt )
    tm.Draw()
    l.Draw()
    c.SaveAs( 'plots/Coef/%s_%sCoef_%s.gif' % ('Pyth_AlpRW', coef,entry))
    

def draw_kin( entry, groups , hand = 'all', rebin = 0, pol = False, mode = None):
    c = ROOT.TCanvas()
    axis = None
    if 'Zmass' == entry: 
        if not pol:
            c.SetLogy()
        else:
            axis = [60,65,70]; axis.extend( [72 + 2*i for i in xrange( (80-70)/2 -1 )])
            axis.extend( [80 + i for i in xrange(19)] )
            axis.extend(  [100 + 2*i for i in xrange (4) ])
            axis.extend([110,115,120,125])
            axis.extend([130+i*20 for i in xrange(3)])
    if 'ZpT' in entry and pol:
        axis = [0]
        axis.extend([10*i for i in xrange(1,6)])
        axis.extend([60+20*i for i in xrange(3)])
        axis.extend([120,200])

    l = ROOT.TLegend(0.77,0.77,0.9,0.92)
    l.SetTextSize(0.034)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    hists = []
    maxi = 0.0
    mini = 0
    
    for j,group in enumerate(groups):
        print group.name
        if group == group1: hf = af
        if group == group2: hf = pf
        if group == group3: hf = afrw
        w = None
        wL = None
        wR = None
        hL = None
        if pol: hand = 'LR'
        for i,sample in enumerate(group):
            wt = ''
            md = ''
            if mode: md = '_%s' % mode
            if 'AlpgenZmassReweighted' == group.name:
                wt = '_rwZmass'
            elif 'AlpgenZpTReweighted' == group.name:
                wt = '_rwZpt'
            if not hand == 'LR':
                h = hf.Get("%s_%s%s_%s%s" % (sample.name, entry, md, hand,wt)).Clone("h%s" % sample.name)
            else:
                print "%s_%s%s_%s%s" % (sample.name, entry, md,  "LH",wt)
                hL = hf.Get("%s_%s%s_%s%s" % (sample.name, entry, md,  "LH",wt)).Clone("hL%s" % sample.name)
                hR = hf.Get("%s_%s%s_%s%s" % (sample.name, entry, md, "RH",wt)).Clone("hR%s" % sample.name)
                h = hf.Get("%s_%s%s_%s%s" % (sample.name, entry, md, "all",wt)).Clone("hA%s" % sample.name)
            if axis:
                h = h.Rebin(len(axis)-1, "rebin_h%s" % sample.name, array('d',axis))
                if hL:
                    hL = hL.Rebin(len(axis)-1, "rebin_hL" , array('d',axis))                     
                    hR = hR.Rebin(len(axis)-1, "rebin_hR" , array('d',axis))                     
            elif rebin:
                h.Rebin(rebin)
                if hL:
                    for hh in [hL,hR]:
                        hh.Rebin(rebin)

            norm = h.Integral()
            h.Scale(1./norm)
            if 'Alpgen' in group.name:
                h.Scale(sample.xsec/XSEC)
            if hL:
                for hh in [hL,hR]:
                    hh.Scale(1./norm)
                    if 'Alpgen' in group.name:
                        hh.Scale(sample.xsec/(XSEC))
            if not w: 
                w = h
                if hL and not wL:
                    wL = hL
                    wR = hR 
            else:
                w.Add(h)
                if wL:
                    wL.Add(hL)
                    wR.Add(hR)

        w.SetXTitle(xtitles[entry])
        if not pol:
            if group.name == 'AlpgenZpTReweighted':
                l.AddEntry(w,'Alpgen Z#tau#tau (p_{T}^{Z} RW)')
            elif group.name == 'AlpgenZmassReweighted':
                l.AddEntry(w,'Alpgen Z#tau#tau (m_{Z} RW)')
            else:
                l.AddEntry(w,group.name)
            if j > 0: 
                w.SetLineColor(9)
                w.SetMarkerColor(9)
            hists.append(w)
            if wL:
                wL.SetLineColor(2)
                wL.SetMarkerColor(2)
                wR.SetLineColor(4)
                wR.SetMarkerColor(4)
                hists.extend([wL,wR])
                l.AddEntry(wL, 'LH')
                l.AddEntry(wR,'RH')
        else:

            if False:
                for xlo,xhi in [ (89,92), (80,100), (70,110),(60,120),(60,180)]:
                    blo = wR.FindBin(xlo)
                    bhi = wR.FindBin(xhi)
                    nR = ErrorFloat(*IntegralError( wR, binlo = blo,binhi = bhi))
                    nL = ErrorFloat(*IntegralError( wL, binlo = blo,binhi = bhi))
                    print '[%.0f, %.0f] %s' %(wR.GetBinLowEdge(blo), wR.GetBinLowEdge(bhi), (nR-nL)/(nR+nL))
                    


            wp = wR.Clone("wp")
            wp.Add(wL,-1)
            wp.Divide(wL+wR)
            wp.SetXTitle(xtitles[entry])
            wp.SetYTitle('Polarization')
            if j>0:
                wp.SetLineColor(9)
                wp.SetMarkerColor(9)
            l.AddEntry(wp,group.name)
            print wp.GetNbinsX()

            hists.append(wp)

    maxi = 0; mini = 0
    for h in hists:
        maxi = max(maxi,h.GetMaximum())
        mini = min(mini,h.GetMinimum())

    draw = ''
    for i,h in enumerate(hists):
        #if maxi == 0: maxi = 0.1
        print h.Integral()
        if not entry == 'Zmass':
            h.SetMinimum(1.4*mini)
            h.SetMaximum( 1.4 * maxi )
        #h.SetMinimum(-0.1)
        #h.SetMaximum(0.1)
        h.Draw(draw)
        if i==0: draw+='SAME'
    tm = ROOT.TPaveText(0.22,0.8,0.48,0.85,"NDC")
    if pol: tm = ROOT.TPaveText(0.22,0.2,0.48,0.35,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    if 'AZ' in entry:
        tt = '80 GeV < m_{Z} < 100 GeV'
    else: tt = ''
    tm.AddText( tt )
    if tt:
        tm.Draw()
    l.Draw()
    if not os.path.exists('plots/%s' % entry): os.mkdir('plots/%s' % entry)
    gname = groups[0].name
    if len(groups)>1: gname += '_%s' % groups[1].name
    gname += '_%s' % entry
    if pol: gname += '_Pol'
    else: gname += '_%s' % hand
    c.SaveAs( 'plots/%s/%s.gif' % (entry,gname))


