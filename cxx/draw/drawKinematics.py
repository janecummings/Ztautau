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

#plotdir = 'plots/TruthStudy/CSReweight'
plotdir = '../plots/'

TrueNp0 = TruthSample(name='Np0', ls = 'mcevt*Np0', xsec = 718.87*1.18*1.0)#, ls = 'Np0')
TrueNp1 = TruthSample(name='Np1', ls = 'mcevt*Np1', xsec = 175.75*1.18*1.0)
TrueNp2 = TruthSample(name='Np2', ls = 'mcevt*Np2', xsec = 58.856*1.18*1.0)
TrueNp3 = TruthSample(name='Np3', ls = 'mcevt*Np3', xsec = 15.667*1.18*1.0)
TrueNp4 = TruthSample(name='Np4', ls = 'mcevt*Np4', xsec = 4.0121*1.18*1.0)
TrueNp5 = TruthSample(name='Np5', ls = 'mcevt*Np5', xsec = 1.256*1.18*1.0)

group = Group('TrueZtautau', [TrueNp0, TrueNp1, TrueNp2, TrueNp3, TrueNp4, TrueNp5])

hf = ROOT.TFile('../output/TruthStudy.root')
#hf = ROOT.TFile('CS_Reweight.root')
flavors = ["up_qq","down_qq","none_qq","up_qg","down_qg","none_qg","up_gg","down_gg","none_gg","all"]
entries = ["Zpt","Zmass","TScostheta","cosa","xfrac","upsilon","taupt","vistaupt","taueta","vistaueta","CosDirect","x1","x2","scale"]
xtitles = ["Z p_{T} [GeV]", "m_{Z} [GeV]", "TauSpinner cos#theta*", "Tau RF cos#theta", "x = p_{T}^{vis}/p_{T}^{#tau}",
            "#Upsilon","#tau p_{T}","#tau p_{T}^{vis}",'#tau #eta','#tau_{vis} #eta',"Direct cos#theta*", "x1", "x2", "PDF scale"]

xtit = dict(zip(entries,xtitles))

modes = {"0":"","1":"1p","3":"#tau#rightarrow#pi#nu","4":"#tau#rightarrow#pi#pi^{0}#nu"}

XSEC = 0.
for sample in group:
    XSEC+=sample.xsec


## Compare Kinematics in Np Samples
def compare_Np( entry, xtitle, mode, hand = "", flavor = "all"):
    cname = "Np_%s" % entry
    if mode not in ["0"]: cname+="_%s" % mode
    if hand: cname += "_%s" % hand
    if not flavor == 'all': cname += '_%s' % flavor 
    c = ROOT.TCanvas()
    if entry in ['Zmass']:
        c.SetLogy()
    if entry in ['x1','x2']:
        c.SetLogx()

    l = ROOT.TLegend(0.77,0.77,0.92,0.92)
    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    l.SetNColumns(2)
    hists = []
    maxi = 0
    pol = ErrorFloat(0,0)
    for i,sample in enumerate(group):
        if hand:
            h = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, hand)).Clone()
        else:
            h = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "LH")).Clone()
            nl = ErrorFloat(*IntegralError(h,of=True))
            hR =hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "RH")).Clone()
            nr = ErrorFloat(*IntegralError(hR,of=True))
            h.Add(hR)

        h.SetLineColor(i+1)
        h.SetMarkerColor(i+1)
        h.SetXTitle(xtitle)
        ptau = (nr-nl)/(nr+nl)
        print ptau
        pol += ErrorFloat(sample.xsec/XSEC,0)*(nr-nl)/(nl+nr)

        h.Rebin(5)
        if not h.Integral(): continue
        if i == 0: 
            hT = h.Clone( "hT")
            hT.Scale( sample.xsec / hT.Integral() )
        else: hT.Add( h, sample.xsec / h.Integral() )
        if entry in ['scale']: h.GetXaxis().SetRangeUser(0,250)

        h.GetXaxis().SetRange(1,19)

        h.Scale( 1 / h.Integral() )
        l.AddEntry(h, sample.name)
        hists.append(h)

        maxi = max(maxi,h.GetMaximum())
        print '**** %s ****' % sample.name
        #h.GetXaxis().SetRangeUser(3,4)
        h.GetXaxis().SetRange( 3,19 )
        # for ib in xrange(h.GetNbinsX()):
        #     print ib,h.GetBinLowEdge(ib)
        # sys.exit()
  
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

    print pol
    return 0

    print "TauSpinner Polarization: %s" % pol

    draw = ''
    hists.reverse()
    for i,h in enumerate(hists):
        h.SetMaximum( 1.4 * maxi )
        h.SetMinimum(0.)
        h.Draw(draw)
        if i==0: draw+='SAME'

    tm = ROOT.TPaveText(0.18,0.8,0.25,0.85,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    tm.AddText( '' )
    tm.Draw() 
    l.Draw()

    #c.SaveAs( '%s/%s.gif+1' % (plotdir,cname))

def compare_Project( entry1, entry2, xtitle, mode, hand = ""):

    entry = entry1
    cname = "Np_%s_%s" % (entry1,entry2)
    if mode not in ["0"]: cname+="_%s" % mode
    if hand: cname += "_%s" % hand
    if not flavor == 'all': cname += '_%s' % flavor 
    c = ROOT.TCanvas()
    if entry in ['Zmass']:
        c.SetLogy()
    if entry in ['x1','x2']:
        c.SetLogx()

    l = ROOT.TLegend(0.77,0.77,0.92,0.92)
    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    l.SetNColumns(2)
    hists = []
    maxi = 0
    pol = ErrorFloat(0,0)

    lo = 60; hi = 120
    cname+="_%s" % lo

    results = []
    ts = []


    for i,sample in enumerate(group):
        if hand:
            h = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry1, entry2, mode, hand)).Clone()
        else:
            h = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry1, entry2, mode, "LH")).Clone()
            h = h.ProjectionX( "p%iL" % i, lo, hi )

            nl = ErrorFloat(*IntegralError(h,of=False))
            hR =hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry1, entry2, mode, "RH")).Clone()
            hR = hR.ProjectionX( "p%iR" % i, lo, hi )
            nr = ErrorFloat(*IntegralError(hR,of=False))
            h.Add(hR)

        h.SetLineColor(i+1)
        h.SetMarkerColor(i+1)
        h.SetXTitle(xtitle)
        ptau = (nr-nl)/(nr+nl)
        print ptau
        ts.append( ( ptau.val, ptau.err))
        pol += ErrorFloat(sample.xsec/XSEC,0)*(nr-nl)/(nl+nr)

        h.Rebin(5)
        if not h.Integral(): continue
        if i == 0: 
            hT = h.Clone( "hT")
            hT.Scale( sample.xsec / hT.Integral() )
        else: hT.Add( h, sample.xsec / h.Integral() )
        if entry in ['scale']: h.GetXaxis().SetRangeUser(0,250)

        h.GetXaxis().SetRange(1,19)

        h.Scale( 1 / h.Integral() )
        l.AddEntry(h, sample.name)
        hists.append(h)

        maxi = max(maxi,h.GetMaximum())
        #print '**** %s ****' % sample.name
        #h.GetXaxis().SetRangeUser(3,4)
        h.GetXaxis().SetRange( 3,19 )
        # for ib in xrange(h.GetNbinsX()):
        #     print ib,h.GetBinLowEdge(ib)
        # sys.exit()
  
        # h.Fit('pol1')
        # fit = h.GetFunction('pol1')
        # p0 = fit.GetParameter(0)
        # p1 = fit.GetParameter(1)
        # p0e = fit.GetParError(0)
        # p1e = fit.GetParError(1)

        # p0 = ErrorFloat(  p0,  p0e  )
        # p1 = ErrorFloat( p1, p1e )
        # two = ErrorFloat( 2, 0)
        # r = ( p1 ) / ( p1 + two * p0)
        ##results.append(( p1 ) / ( p1 + two * p0))
        #results.append( (r.val, r.err))
    print pol
    #return 0

    print "TauSpinner Polarization: %s" % pol

    draw = ''
    hists.reverse()
    for i,h in enumerate(hists):
        h.SetMaximum( 1.4 * maxi )
        h.SetMinimum(0.)
        h.Draw(draw)
        if i==0: draw+='SAME'

    tm = ROOT.TPaveText(0.18,0.8,0.25,0.85,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    tm.AddText( '' )
    tm.Draw() 
    l.Draw()


   # print ts
   # print results


    # for r in results:
    #     print r

    hT.GetXaxis().SetRange(1,19)
    hT.Scale( 1 / hT.Integral() )
    hT.GetXaxis().SetRange( 3,19 )
    hT.Fit('pol1')
    fit = hT.GetFunction('pol1')
    p0 = fit.GetParameter(0)
    p1 = fit.GetParameter(1)
    p0e = fit.GetParError(0)
    p1e = fit.GetParError(1)

    p0 = ErrorFloat(  p0,  p0e  )
    p1 = ErrorFloat( p1, p1e )
    two = ErrorFloat( 2, 0)
    print ( p1 ) / ( p1 + two * p0)

    #c.SaveAs( '%s/%s.gif+1' % (plotdir,cname))

def compare_Inclusive( entry, xtitle, mode, hand = "", flavor = "all"):
    cname = "Inclusive_%s" % flavor
    if mode not in ["0"]: cname+="_%s" % mode
    if hand: cname += "_%s" % hand

    c = ROOT.TCanvas()
    if entry in ['Zmass','scale']:
        c.SetLogy()
    if entry in ['x1','x2']:
        c.SetLogx()

    l = ROOT.TLegend(0.77,0.82,0.92,0.92)
    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)


    if hand:
        h0 = hf.Get("%s_%s_%s_%s_%s" % ("Np0", entry, flavor, mode, hand)).Clone()
    else:
        h0 = hf.Get("%s_%s_%s_%s_%s" % ("Np0", entry, flavor, mode, "LH")).Clone()
        h0R =hf.Get("%s_%s_%s_%s_%s" % ("Np0", entry, flavor, mode, "RH")).Clone()
        h0.Add(h0R)
    h0.Scale( 1./h0.Integral())

    h = None
    for i,sample in enumerate(group):
        if hand:
            hi = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, hand)).Clone()
        else:
            hi = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "LH")).Clone()
            hiR =hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "RH")).Clone()
            hi.Add(hiR)
        hi.Scale( sample.xsec *  hi.Integral( 0,hi.GetNbinsX()+1  )   )
        if not h: 
            h = hi
        else: h.Add( hi )
        if not h.Integral(): continue
        if entry in ['scale']: h.GetXaxis().SetRangeUser(0,250)
    h.Scale( 1 / h.Integral() )
    h.SetLineColor(46)
    h.SetMarkerColor(46)
    h.SetXTitle( xtitle )

    maxi = max(h0.GetMaximum(),h.GetMaximum())

    draw = ''
    h.SetMaximum( 1.4 * maxi )
    h.Draw(draw)
    h0.Draw("SAME")
    l.AddEntry(h0,"Np0")
    l.AddEntry(h,"Inclusive")

    tm = ROOT.TPaveText(0.18,0.8,0.25,0.85,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    tm.AddText( '' )
    #tm.Draw() 
    l.Draw()

    c.SaveAs( '%s/%s.gif+1' % (plotdir,cname))

def compare_polarization( entry, xtitle, mode, flavor = 'all'):
    cname = "Np_Pol_%s_diff" % flavor
    c = ROOT.TCanvas()
    if entry in ['x1','x2',]:
        c.SetLogx()

    l = ROOT.TLegend(0.77,0.77,0.92,0.92)
    if entry in ['cosa','xfrac','upsilon']:
        l = ROOT.TLegend( 0.42, 0.77, 0.57, 0.92 )
    elif entry in ['Zpt', 'x1','x2']:
        l = ROOT.TLegend( 0.77, 0.23, 0.92, 0.38 )
    elif entry in ['taupt']:
        l = ROOT.TLegend( 0.23, 0.23, 0.38, 0.38 )
    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    l.SetNColumns(2)
    hists = []
    maxi = 0
    mini = 0
    axis = None
    if entry in ['Zpt']:
        axis = [ 10*i for i in xrange(20)]
    if entry in ['Zmass','scale']:
        axis = [60,65,70]; axis.extend( [72 + 2*i for i in xrange( (80-70)/2 -1 )])
        axis.extend( [80 + i for i in xrange(19)] )
        axis.extend(  [100 + 2*i for i in xrange (4) ])
        axis.extend([110,115,120,125])
        axis.extend([130+i*20 for i in xrange(3)])


    for i,sample in enumerate(group):
            hL = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "LH")).Clone()
            h = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "RH")).Clone()
            hT = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "RH")).Clone()
            hT.Add(hL)
            h.Add( hL, -1 )
            #h.Rebin(2)
            #hT.Rebin(2)
            if axis:
                h = h.Rebin(len(axis) -1 , "rebin_h", array('d',axis))
                hT = hT.Rebin(len(axis) -1 , "rebin_hT", array('d',axis))
            else:
                try:
                    h.Rebin( 6 )
                    hT.Rebin( 6 )
                except:
                    h.Rebin( 5)
                    hT.Rebin( 5 )

            h.Divide(hT)

            if i==0: h0 = h.Clone("h0")

            h.Add(h0,-1)
            h.SetLineColor(i+1)
            h.SetMarkerColor(i+1)
            h.SetXTitle(xtitle)
            h.SetYTitle("P_{#tau} - Np0 P_{#tau}")
            #if not h.Integral(): continue
            if entry in ['scale']: h.GetXaxis().SetRangeUser(0,250)
            if entry in ['xfrac']: 
                print 'here'
                h.GetXaxis().SetRangeUser(0,1.0)
            if entry in ['x1','x2']: 
                h.GetXaxis().SetRangeUser(0.,0.5)
                #h.SetMinimum(-0.3)
            #h.Scale( 1 / h.Integral() )
            l.AddEntry(h, sample.name)
            hists.append(h)
            maxi = max(maxi,h.GetMaximum())
            mini = min(mini,h.GetMinimum())
    if maxi == 0: maxi = 0.1
    draw = ''
    hists.reverse()
    for i,h in enumerate(hists):
        h.SetMaximum( 1.4 * maxi )
        h.SetMinimum( 1.4 * mini )
        #h.SetMinimum(0.)
        h.Draw(draw)
        if i==0: draw+='SAME'
    tm = ROOT.TPaveText(0.18,0.8,0.25,0.85,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    tm.AddText( '' )
    #tm.Draw() 
    l.Draw()

    c.SaveAs( '%s/%s.gif+1' % (plotdir,cname))    

def draw_Inclusive( entry, xtitle, mode, hand = "", flavor = "all"):
    cname = "Incl_%s" % flavor
    if mode not in ["0"]: cname+="_%s" % mode
    if hand: cname += "_%s" % hand

    c = ROOT.TCanvas()
    if entry in ['Zmass','scale']:
        c.SetLogy()
    if entry in ['x1','x2']:
        c.SetLogx()

    l = ROOT.TLegend(0.77,0.77,0.92,0.92)
    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)

    h = None;
    for i,sample in enumerate(group):

        hi = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "LH")).Clone()
        hiL = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "LH")).Clone()
        hiR =hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "RH")).Clone()
        hi.Add(hiR)
        nT = hi.Integral()
        hi.Scale( sample.xsec / (XSEC * hi.Integral( 0,hi.GetNbinsX()+1)  ))
        hiL.Scale( sample.xsec / (XSEC * nT))
        hiR.Scale( sample.xsec / (XSEC * nT))

        if entry in ['Zmass']:
            for hst in [hi,hiL,hiR]:
                hst.GetXaxis().SetRange( 60, 200)


        if not h: 
            h = hi; hL = hiL; hR = hiR
        else: 
            h.Add( hi ); hL.Add(hiL); hR.Add(hiR)
        if not h.Integral(): continue
        if entry in ['scale']: h.GetXaxis().SetRangeUser(0,250)
    #h.Scale( 1 / h.Integral() )
    h.SetLineColor(1)
    h.SetMarkerColor(1)
    hL.SetMarkerColor(2); hL.SetLineColor(2)
    hR.SetMarkerColor(4); hR.SetLineColor(4)
    h.SetXTitle( xtitle )


    draw = ''
    maxi = max( h.GetMaximum(), hL.GetMaximum(), hR.GetMaximum())
    h.SetMaximum( 1.4 * maxi)
    h.Draw(draw)
    hR.Draw("SAME")
    hL.Draw("SAME")

    l.AddEntry(h,"Z#rightarrow#tau#tau")
    l.AddEntry(hL, 'LH')
    l.AddEntry(hR, 'RH' )
    tm = ROOT.TPaveText(0.18,0.8,0.25,0.85,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    tm.AddText( '' )
    #tm.Draw() 
    l.Draw()

    c.SaveAs( '%s/%s.gif+1' % (plotdir,cname))

def inclusive_polarization( entry, xtitle, mode, flavor = 'all', low = None, high = None):
    cname = "Incl_Pol_%s" % flavor
    c = ROOT.TCanvas()
    if entry in ['x1','x2',]:
        c.SetLogx()

    l = ROOT.TLegend(0.77,0.77,0.92,0.92)
    if entry in ['cosa','xfrac','upsilon']:
        l = ROOT.TLegend( 0.42, 0.77, 0.57, 0.92 )
    elif entry in ['Zpt', 'x1','x2']:
        l = ROOT.TLegend( 0.77, 0.23, 0.92, 0.38 )
    elif entry in ['taupt']:
        l = ROOT.TLegend( 0.23, 0.23, 0.38, 0.38 )
    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    l.SetNColumns(2)
    axis = None
    if entry in ['Zpt']:
        axis = [ 10*i for i in xrange(20)]
    if entry in ['Zmass','scale']:
        axis = [60,65,70]; axis.extend( [72 + 2*i for i in xrange( (80-70)/2 -1 )])
        axis.extend( [80 + i for i in xrange(19)] )
        axis.extend(  [100 + 2*i for i in xrange (4) ])
        axis.extend([110,115,120,125])
        axis.extend([130+i*20 for i in xrange(3)])
        axis.append(200)


    for i,sample in enumerate(group):
            hL = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "LH")).Clone()
            hR = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "RH")).Clone()


            if axis:
                hL = hL.Rebin(len(axis) -1 , "rebin_hL", array('d',axis))
                hR = hR.Rebin(len(axis) -1 , "rebin_hR", array('d',axis))
            else:
                hL.Rebin( 4 )
                hR.Rebin( 4 )
            nT = hL.Integral() + hR.Integral()
            hL.Scale( sample.xsec / (nT * XSEC) )
            hR.Scale( sample.xsec / (nT * XSEC)  )

            if i==0: 
                hLeft = hL.Clone("hLeft")
                hRight = hR.Clone("hRight")
            else:
                hLeft.Add(hL); hRight.Add(hR)
    draw = ''
    blo = 1; bhi = hLeft.GetNbinsX()
    if low:
        blo = hLeft.FindBin( low );
    if high:
        bhi = hLeft.FindBin( high );
    nL = ErrorFloat(*IntegralError(hLeft,of=False, binlo = blo, binhi = bhi))
    nR = ErrorFloat(*IntegralError(hRight,of=False, binlo = blo, binhi = bhi))
    print '[%s,%s] = %s' % ( hLeft.GetBinLowEdge(blo), hRight.GetBinLowEdge( bhi + 1),( nR-nL )/ (nL+nR))


    hT = hRight.Clone("hTotal")
    hT.Add( hLeft )
    hRight.Add( hLeft, -1)
    hRight.Divide(hT)

    hRight.SetYTitle("Polarization")
    hRight.SetXTitle(xtitle)

    hRight.Draw()

    tm = ROOT.TPaveText(0.18,0.8,0.25,0.85,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    tm.AddText( '' )
    #tm.Draw() 
    #l.Draw()

    c.SaveAs( '%s/%s.gif+1' % (plotdir,cname))    

def compare_Q2_polarization( entry, xtitle, mode, lo,hi,flavor = 'all'):
    cname = "Np_PolQ2_%s" % flavor
    c = ROOT.TCanvas()
    if entry in ['x1','x2',]:
        c.SetLogx()

    l = ROOT.TLegend(0.77,0.77,0.92,0.92)
    if entry in ['cosa','xfrac','upsilon']:
        l = ROOT.TLegend( 0.42, 0.77, 0.57, 0.92 )
    elif entry in ['Zpt', 'x1','x2']:
        l = ROOT.TLegend( 0.77, 0.23, 0.92, 0.38 )
    elif entry in ['taupt']:
        l = ROOT.TLegend( 0.23, 0.23, 0.38, 0.38 )
    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    l.SetNColumns(2)
    hists = []
    maxi = 0
    mini = 0
    axis = None
    if entry in ['Zpt']:
        axis = [ 10*i for i in xrange(20)]
    if entry in ['Zmass','scale']:
        axis = [60,65,70]; axis.extend( [72 + 2*i for i in xrange( (80-70)/2 -1 )])
        axis.extend( [80 + i for i in xrange(19)] )
        axis.extend(  [100 + 2*i for i in xrange (4) ])
        axis.extend([110,115,120,125])
        axis.extend([130+i*20 for i in xrange(3)])


    for i,sample in enumerate(group):
            if lo < 100 and i > 3: continue
            entryname = 'scale_%s' % entry
            hL = hf.Get("%s_%s_%s_%s" % (sample.name, entryname, mode, "LH")).Clone()
            h = hf.Get("%s_%s_%s_%s" % (sample.name, entryname, mode, "RH")).Clone()
            hT = hf.Get("%s_%s_%s_%s" % (sample.name, entryname, mode, "RH")).Clone()
            hL = hL.ProjectionY( "pL%i" % i, lo, hi  )
            h = h.ProjectionY( "p%i" %i, lo, hi)
            hT = hT.ProjectionY( "pT%i" %i, lo, hi)

            hT.Add(hL)
            h.Add( hL, -1 )
            #h.Rebin(2)
            #hT.Rebin(2)
            if axis:
                h = h.Rebin(len(axis) -1 , "rebin_h", array('d',axis))
                hT = hT.Rebin(len(axis) -1 , "rebin_hT", array('d',axis))
            else:
                h.Rebin( 4 )
                hT.Rebin( 4 )

            h.Divide(hT)
            h.SetLineColor(i+1)
            h.SetMarkerColor(i+1)
            h.SetXTitle(xtitle)
            h.SetYTitle("Polarization")
            if not h.Integral(): continue
            if entry in ['scale']: h.GetXaxis().SetRangeUser(0,250)
            if entry in ['x1','x2']: 
                h.GetXaxis().SetRangeUser(0.,0.5)
                h.SetMinimum(-0.3)
            #h.Scale( 1 / h.Integral() )
            l.AddEntry(h, sample.name)
            hists.append(h)
            maxi = max(maxi,h.GetMaximum())
            mini = min(mini,h.GetMinimum())

    draw = ''
    hists.reverse()
    for i,h in enumerate(hists):
        h.SetMaximum( 1.4 * maxi )
        h.SetMinimum( 1.4 * mini )
        #h.SetMinimum(0.)
        h.Draw(draw)
        if i==0: draw+='SAME'
    tm = ROOT.TPaveText(0.28,0.8,0.45,0.85,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    tm.AddText( 'Q^2 = [%s,%s]'%(lo,hi) )
    tm.Draw() 
    l.Draw()

    c.SaveAs( '%s/%s.gif+1' % (plotdir,cname))    

def project_Q2( entry, xtitle, mode, lo, hi, flavor = 'all', hand = '' ):
    cname = "Np_all_Q2" #% flavor
    if mode not in ["0"]: cname+="_%s" % mode
    if hand: cname += "_%s" % hand

    c = ROOT.TCanvas()
    if entry in ['Zmass']:
        c.SetLogy()
    if entry in ['x1','x2']:
        c.SetLogx()

    l = ROOT.TLegend(0.77,0.77,0.92,0.92)
    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    l.SetNColumns(2)
    hists = []
    maxi = 0
    for i,sample in enumerate(group):
        #if i>3: continue
        if lo < 100 and i > 3: continue
        entryname = 'scale_%s' % entry
        if hand:
            h = hf.Get("%s_%s_%s_%s" % (sample.name, entryname, mode, hand)).Clone()
        else:
            h = hf.Get("%s_%s_%s_%s" % (sample.name, entryname,  mode, "LH")).Clone()
            hR =hf.Get("%s_%s_%s_%s" % (sample.name, entryname, mode, "RH")).Clone()
            h.Add(hR)
        h = h.ProjectionY( "p%i" % i, lo, hi )
        h.SetLineColor(i+1)
        h.SetMarkerColor(i+1)
        h.SetXTitle(xtitle)
        if not h.Integral(): continue
        if entry in ['scale']: h.GetXaxis().SetRangeUser(0,250)
        #h.Rebin(4)
        h.Scale( 1 / h.Integral() )
        l.AddEntry(h, sample.name)
        hists.append(h)
        maxi = max(maxi,h.GetMaximum())

    draw = ''
    hists.reverse()
    for i,h in enumerate(hists):
        h.SetMaximum( 1.4 * maxi )
        #h.SetMinimum(0.)
        h.Draw(draw)
        if i==0: draw+='SAME'
    tm = ROOT.TPaveText(0.28,0.8,0.35,0.85,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    tm.AddText( 'Q^2 = [%s,%s]' % (lo,hi) )
    tm.Draw() 
    l.Draw()

    c.SaveAs( '%s/%s.gif+1' % (plotdir,cname))

def plot_TH2( entry1, entry2, xtitle, ytitle, mode, hand = ''):

    for i,sample in enumerate(group):
        cname = "%s_%s" % (entry1,entry2)

        if mode not in ["0"]: cname+="_%s" % mode
        if hand: cname += "_%s" % hand

        c = ROOT.TCanvas()
        c.SetMargin(.2,.2,.2,.1)
        #if entry1 in ['Zmass']:
        #    c.SetLogy()
        #if entry1 in ['x1','x2']:
        #    c.SetLogx()

        l = ROOT.TLegend(0.85,0.9,0.92,0.96)
        l.SetTextSize(0.035)
        l.SetBorderSize(1)
        l.SetFillColor(0)

        if hand:
            h = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry1, entry2, mode, hand)).Clone()
        else:
            h = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry1, entry2, mode, "LH")).Clone()
            hR =hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry1, entry2, mode, "RH")).Clone()
            h.Add(hR)
        h.SetLineColor(i+1)
        h.SetMarkerColor(i+1)
        h.SetXTitle(xtitle)
        h.SetYTitle(ytitle)
        if not h.Integral(): continue
        if entry1 in ['scale']: h.GetXaxis().SetRangeUser(0,250)
        h.Scale( 1 / h.Integral() )
        l.AddEntry(h, sample.name)
        


        h.Draw("COLZ")
        tm = ROOT.TPaveText(0.18,0.8,0.25,0.85,"NDC")
        tm.SetBorderSize(0)
        tm.SetFillColor(0)

        tm.AddText( '' )
        #tm.Draw() 
        l.Draw()

        c.SaveAs( '%s/TH2/%s.gif+1' % (plotdir,cname))

def ud_ratio( entry, xtitle, mode , hand=""):
    cname = "Np_qqud" 
    if hand: cname += "_%s" % hand
    c = ROOT.TCanvas()
    if entry in ['Zmass','scale']:
        c.SetLogy()
    if entry in ['x1','x2']:
        c.SetLogx()
        c.SetLogy()

    l = ROOT.TLegend(0.77,0.77,0.92,0.92)
    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    l.SetNColumns(2)
    hists = []
    maxi = 0
    axis = None
    if entry in ['Zmass','scale']:
        axis = [60,65,70]; axis.extend( [72 + 2*i for i in xrange( (80-70)/2 -1 )])
        axis.extend( [80 + i for i in xrange(19)] )
        axis.extend(  [100 + 2*i for i in xrange (4) ])
        axis.extend([110,115,120,125])
        axis.extend([130+i*20 for i in xrange(3)])
    for i,sample in enumerate(group):

        if hand:
            h = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, 'up_qq', "0", hand)).Clone()
            h_down = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, 'down_qq', "0", hand)).Clone()
        else:

            h = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, 'up_qq', "0", "LH")).Clone()
            hR_up =hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, 'up_qq', "0", "RH")).Clone()
            h.Add(hR_up)
            h_down = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, 'down_qq', "0", "LH")).Clone()
            hR_down =hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, 'down_qq', "0", "RH")).Clone()
            h_down.Add(hR_down)
        if axis:
            h_down = h_down.Rebin( len(axis)-1, "h_down_rb", array("d",axis))
            h = h.Rebin( len(axis)-1, "h_rb", array("d",axis))
        else:
            h_down.Rebin(4)
            h.Rebin(4)
        h.Divide(h_down)
        h.SetLineColor(i+1)
        h.SetMarkerColor(i+1)
        h.SetXTitle(xtitle)
        if not h.Integral(): continue
        if entry in ['scale']: h.GetXaxis().SetRangeUser(60,150)
        #h.Scale( 1 / h.Integral() )

        l.AddEntry(h, sample.name)
        hists.append(h)
        maxi = max(maxi,h.GetMaximum())

    draw = ''
    hists.reverse()
    for i,h in enumerate(hists):
        h.SetMaximum( 1.4 * maxi )
        #h.SetMinimum(0.)
        h.Draw(draw)
        if i==0: draw+='SAME'
    tm = ROOT.TPaveText(0.2,0.82,0.25,0.86,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.045)
    tm.AddText( hand )
    tm.Draw() 
    l.Draw()

    c.SaveAs( '%s/%s.gif+1' % (plotdir,cname))

### Compare flavor distributions for fixed Np sample
def compare_flavors( entry, xtitle = '', modes = modes, hand = "" ):
    cname = '%s_flavors' % entry
    if hand: cname += '_%s' % hand
    for sample in group:
        print "### Sample %s" % sample.name
        c = ROOT.TCanvas()
        #if 'mass' not in entry: continue
        #c.SetLogy()
        l = ROOT.TLegend(0.8,0.7,0.9,0.9)
        l.SetTextSize(0.035)
        l.SetBorderSize(0)
        l.SetFillColor(0)
        hists = []
        maxi = 0
        for i,flavor in enumerate(flavors):
            if hand:
                h = hf.Get("%s_%s_%s_%s" % (sample.name, entry, flavor, "0", hand)).Clone()
            else:
                #print "%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", "LH")
                h = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", "LH")).Clone()
                hR =hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", "RH")).Clone()
                h.Add(hR)
            hT = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, "all", "0", "LH")).Clone()
            hRT =hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, "all", "0", "RH")).Clone()
            hT.Add(hRT)
            nT = hT.Integral(0,hT.GetNbinsX()+1)
            nf = h.Integral(0,h.GetNbinsX()+1)
            print "%s: %s, %.1f" % (flavor, h.Integral(0,h.GetNbinsX()+1), 100.*nf/nT)
            continue
            if flavor == 'all': i = 5
            h.SetLineColor(i+1)
            h.SetMarkerColor(i+1)
            h.SetXTitle(xtitle)
            if not h.Integral(): continue
            h.Scale( 1 / h.Integral() )
            l.AddEntry(h, flavor )
            hists.append(h)
            maxi = max(maxi,h.GetMaximum())
        c.SetLogy()
        draw = ''
        hists.reverse()
        for i,h in enumerate(hists):
            #h.SetMaximum( 1.4 * maxi )
            #h.SetMinimum(0.)
            h.Draw(draw)
            if i==0: draw+='SAME'
        tm = ROOT.TPaveText(0.18,0.8,0.25,0.85,"NDC")
        tm.SetBorderSize(0)
        tm.SetFillColor(0)
        tm.SetTextSize(0.04)
        tm.AddText( sample.name )
        tm.Draw() 
        l.Draw()

        c.SaveAs( 'plots/%s.gif+1' % cname)

def get_flavor_groups ( sample, entry, hand = '', what = 'up', mode = "0"):
    hist = None
    for i,flavor in enumerate(flavors):
        print what, flavor
        if what not in flavor: continue 
        if hand:
            h = hf.Get("%s_%s_%s_%s" % (sample.name, entry, flavor, mode, hand)).Clone()
        else:
            print "%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "LH")
            h = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "LH")).Clone()
            hR =hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "RH")).Clone()
            h.Add(hR)
        if not hist: hist = h
        else: hist.Add(h)
    return hist

def get_flavor_pola ( sample, entry, hand = '', what = 'up', mode = "0"):
    hLeft = None; hRight = None
    for i,flavor in enumerate(flavors):
        if what not in flavor: continue 
        hL = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "LH")).Clone()
        hR =hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, mode, "RH")).Clone()
        rb = 2
        hL.Rebin(rb)
        hR.Rebin(rb)
        if not hLeft:
            hLeft = hL
            hRight = hR
        else:
            hLeft.Add(hL)
            hRight.Add(hR)
    hSum = hRight.Clone("sum")
    hDiff = hRight.Clone("diff")
    hDiff.Add( hLeft, -1 )
    hSum.Add(hLeft)
    hDiff.Divide(hSum)
    return hDiff

def compare_flavor_groups( entry, xtitle = '', mode = "0", hand = ""):
    cname = '%s_flavor_pola' % entry
    if hand: cname += '_%s' % hand
    for sample in group:
        c = ROOT.TCanvas()
        #if 'mass' in entry: continue
        #c.SetLogy()
        l = ROOT.TLegend(0.8,0.72,0.9,0.92)
        l.SetTextSize(0.035)
        l.SetBorderSize(0)
        l.SetFillColor(0)
        hists = []
        maxi = 0; mini = 0;
        for i,flavor in enumerate(['up','down','qq','qg','gg']):
        #for i,flavor in enumerate(['all']):
            #h = get_flavor_groups( sample, entry, hand = hand, what = flavor , mode = mode)
            h = get_flavor_pola( sample, entry, hand = hand, what = flavor , mode = mode)
            #if flavor == 'all': i = 5
            if not h.Integral(): continue
            h.SetLineColor(i+1)
            h.SetMarkerColor(i+1)
            h.SetXTitle(xtitle)
            h.SetYTitle("Polarization")

            #h.Scale( 1 / h.Integral() )
            l.AddEntry(h, flavor )
            

            h.GetXaxis().SetRangeUser(0,100)
            hists.append(h)
            maxi = max(maxi,h.GetMaximum())
            mini = min(mini, h.GetMinimum())

        #c.SetLogy()
        draw = ''
        hists.reverse()
        for i,h in enumerate(hists):
            h.SetMaximum( 1.5 * maxi )
            h.SetMinimum( 1.2 * mini )
            h.Draw(draw)
            if i==0: draw+='SAME'
        tm = ROOT.TPaveText(0.18,0.8,0.25,0.92,"NDC")
        tm.SetBorderSize(0)
        tm.SetFillColor(0)
        tm.SetTextSize(0.04)
        tm.AddText( sample.name )
        if not mode == '0':
            tm.AddText( modes[mode] )
        tm.Draw() 
        l.Draw()
        c.SaveAs( 'plots/%s.gif+1' % cname)

### Compare Np samples for fixed flavor 
def compare_NpFlavors( entry, xtitle = '', modes = modes, hand = "" ):
    cname = 'NpFlavors'
    if hand: cname += '_%s' % hand
    for flavor in flavors:
        c = ROOT.TCanvas()
        l = ROOT.TLegend(0.8,0.7,0.9,0.92)
        l.SetTextSize(0.034)
        l.SetBorderSize(0)
        l.SetFillColor(0)
        hists = []
        maxi = 0

        for i,sample in enumerate(group):
            if hand:
                h = hf.Get("%s_%s_%s_%s" % (sample.name, entry, flavor, "0", hand)).Clone()
            else:
                print "%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", "LH")
                h = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", "LH")).Clone()
                hR =hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", "RH")).Clone()
                h.Add(hR)
            h.SetLineColor(i+1)
            h.SetMarkerColor(i+1)
            h.SetXTitle(xtitle)
            if not h.Integral(): continue
            h.Scale( 1 / h.Integral() )
            l.AddEntry(h, sample.name )
            hists.append(h)
            maxi = max(maxi,h.GetMaximum())

        draw = ''
        hists.reverse()
        for i,h in enumerate(hists):
            h.SetMaximum( 1.4 * maxi )
            h.Draw(draw)
            if i==0: draw+='SAME'
        tm = ROOT.TPaveText(0.18,0.8,0.25,0.85,"NDC")
        tm.SetBorderSize(0)
        tm.SetFillColor(0)
        tm.SetTextSize(0.04)
        tm.AddText( flavor )
        if not flavor == 'all': tm.Draw() 
        l.Draw()

        c.SaveAs( 'plots/%s.gif+1' % cname)
 
def compare_Np_UpDown( entry, xtitle = '', modes = modes, hand = "" ):
    flavor = 'u/d'
    cname = 'Np_%s' % 'ud'
    if hand: cname += '_%s' % hand
    c = ROOT.TCanvas()
    #c.SetLogx()

    l = ROOT.TLegend(0.8,0.7,0.9,0.92)
    l.SetTextSize(0.034)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    hists = []
    maxi = 0
    mini = 0
    axis = [60+4*i for i in xrange(15)]
    axis.extend([120 + 20*i for i in xrange(5)] )
    #print axis   
    for i,sample in enumerate(group):
        h_up = get_flavor_groups( sample, entry, hand = hand, what = 'up' , mode = '0')
        h_down = get_flavor_groups( sample, entry, hand = hand, what = 'down' , mode = '0')
        h_up = h_up.Rebin( len(axis) - 1, "rebin_Hup_%s" % sample.name, array('d',axis)  )
        h_down = h_down.Rebin( len(axis) - 1, "rebin_Hdown_%s" % sample.name, array('d',axis)  )
        #if not h.Integral(): continue
        #total = h_up.Integral() + h_down.Integral()
        h_up.Divide(h_down)
        h = h_up
        h.SetLineColor(i+1)
        h.SetMarkerColor(i+1)
        h.SetXTitle(xtitle)

        #h = h.Rebin( len(axis) - 1, "rebin_%s" % sample.name, array('d',axis)  )
        #h.GetXaxis().SetRangeUser(60,200)
        #h.Scale( 1 / total )
        l.AddEntry(h, sample.name )
        hists.append(h)
        maxi = max(maxi,h.GetMaximum())

    draw = ''
    #hists.reverse()


    for i,h in enumerate(hists):
        h.SetMaximum( 1.4 * maxi )
        h.SetMinimum( mini )
        #h.GetXaxis().SetRangeUser(60,200)
        
        
        h.Draw(draw)
        if i==0: draw+='SAME'
    tm = ROOT.TPaveText(0.18,0.8,0.25,0.85,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    tm.AddText( flavor )
    if not flavor == 'all': tm.Draw() 
    l.Draw()

    c.SaveAs( 'plots/%s.gif+1' % cname)

def flavor_ratios( entry, xtitle = '', hand = '' ):
    cname = 'FlavorRatio'
    if hand: cname += '_%s' % hand
    c = ROOT.TCanvas()
    l = ROOT.TLegend(0.8,0.7,0.9,0.92)
    l.SetTextSize(0.034)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    hists = []
    maxi = 0
    h_up = None; h_down = None
    for i,sample in enumerate(group):
        for flavor in flavors:
            if hand:
                if 'up' in flavor:
                    if not h_up : h_up = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", hand)).Clone()
                    else: 
                        h_up.Add(hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", hand)).Clone())
                if 'down' in flavor:
                    if not h_down : h_down = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", hand)).Clone()
                    else: 
                        h_down.Add(hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", hand)).Clone())
            else:   
                print flavor
                if 'up' in flavor:
                    if not h_up:
                        #print "%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", "LH")
                        h_up = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", "LH")).Clone()
                        h_up.Add(hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", "RH")).Clone())
                        #print h_up
                    else:
                        h_up.Add(hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", "RH")).Clone())
                        h_up.Add(hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", "LH")).Clone())
                if 'down' in flavor:
                    if not h_down:
                        h_down = hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", "LH")).Clone()
                        h_down.Add(hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", "RH")).Clone())
                    else:
                        h_down.Add(hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", "RH")).Clone())
                        h_down.Add(hf.Get("%s_%s_%s_%s_%s" % (sample.name, entry, flavor, "0", "LH")).Clone())
        h = h_up.Clone(sample.name); h.Divide(h_down)
        #print i,sample.name
        h.SetLineColor(i+1)
        h.SetMarkerColor(i+1)
        h.SetXTitle(xtitle)
            
        #h.Scale( 1 / h.Integral() )
        l.AddEntry(h, sample.name )
        hists.append(h)
        maxi = max(maxi,h.GetMaximum())

    draw = ''
    hists.reverse()
    for i,h in enumerate(hists):
        h.SetMaximum( 1.1 * maxi )
        h.Draw(draw)
        if i==0: draw+='SAME'
    tm = ROOT.TPaveText(0.18,0.8,0.25,0.85,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    tm.AddText( flavor )
    if not flavor == 'all': tm.Draw() 
    l.Draw()

    c.SaveAs( 'plots/%s.gif+1' % cname)