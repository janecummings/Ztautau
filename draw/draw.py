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

#### Atlas Style ####
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.LoadMacro("AtlasStyle.C")
ROOT.gROOT.LoadMacro("AtlasLabels.C")
ROOT.SetAtlasStyle()
ROOT.gStyle.SetFrameBorderSize(1)

## mode fill style
modefill = {'1p0n':3003, '1p1n':3012,'1pXn':3001,'Xp':1}

def plot_stat(spec,entry,region,test='',ratio=True,k_entry = upsilon, verbose = False, SF = 'WQ',noHands=True, syst = None, crspec = None):
    ## Spec
    stream = spec.stream
    PLOTDIR = spec.plotdir 
    if not crspec: crspec = spec
    ## ROOT Files
    rf = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,stream.name))
    kCR = ROOT.TFile('%s/%s_regions.root' % (crspec.OUTPUTDIR,stream.name))
    if syst:
        sysf = ROOT.TFile('%s/%s_%s.root' % (spec.OUTPUTDIR,stream.name, syst))

    # for scaling left and right handed distributions
    kZL = ErrorFloat(1.,0.); kZR = ErrorFloat(1.,0)
    up = 'UP'; down = 'DOWN'
    #up = up.lower(); down = down.lower()

    # scale factors
    kW_OS = ErrorFloat(1.,0.); kW_OS_UP = ErrorFloat(1.,0.); kW_OS_DOWN = ErrorFloat(1.,0.) 
    kW_SS = ErrorFloat(1.,0.); kW_SS_UP = ErrorFloat(1.,0.); kW_SS_DOWN = ErrorFloat(1.,0.) 
    rQCD = ErrorFloat(1.,0.); rQCD_UP = ErrorFloat(1.,0.); rQCD_DOWN = ErrorFloat(1.,0.) 
    if 'W' in SF:
        kW_OS = hist.calc_kX(kCR,stream,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, kZL = kZL, kZR = kZR )
        kW_SS = hist.calc_kX(kCR,stream,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose)
        if syst:
            kW_OS_UP = hist.calc_kX(rf,stream,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, kZL = kZL, kZR = kZR, sys = syst, sf = sysf, ud = up )
            kW_SS_UP = hist.calc_kX(rf,stream,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose, sys = syst, sf =sysf, ud = up)
            kW_OS_DOWN = hist.calc_kX(rf,stream,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, kZL = kZL, kZR = kZR , sys = syst, sf =sysf, ud = down)
            kW_SS_DOWN = hist.calc_kX(rf,stream,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose, sys = syst, sf = sysf, ud = down)
    if 'Q' in SF:
        rQCD = hist.calc_rQCD(kCR,stream,k_entry,'QCD',kW_OS = kW_OS, kW_SS = kW_SS, verbose = verbose)
        if syst:
            rQCD_UP = hist.calc_rQCD(kCR,stream,k_entry,'QCD',kW_OS = kW_OS_UP, kW_SS = kW_SS_UP, verbose = verbose, sys = syst, sf = sysf, ud = up)
            rQCD_DOWN = hist.calc_rQCD(kCR,stream,k_entry,'QCD',kW_OS = kW_OS_DOWN, kW_SS = kW_SS_DOWN, verbose = verbose,  sys = syst, sf = sysf, ud = down)

    if syst:
        print 'kW_OS: %s' % kW_OS
        print 'kW_OS_UP: %s' % kW_OS_UP
        print 'kW_OS_DOWN: %s' % kW_OS_DOWN
        print 'kW_SS: %s' % kW_SS
        print 'kW_SS_UP: %s' % kW_SS_UP
        print 'kW_SS_DOWN: %s' % kW_SS_DOWN
        print 'rQCD: %s' % rQCD
        print 'rQCD_UP: %s' % rQCD_UP
        print 'rQCD_DOWN: %s' % rQCD_DOWN


    # rebin histos
    rb = entry.rebin
    
    # Canvas
    cname = '%s_%s_%s' % (stream.name,entry.name,region)
    pname = '%s_%s' % (stream.name,region.rstrip('_OS'))
    if syst:
        cname += '_%s' % syst
        pname += '_%s' % syst
    c = ROOT.TCanvas(cname,cname,1200,1000)
    if entry.blind: ratio = False

    #Split for ratio plot
    c.Divide(1,2)
    c.cd(1).SetPad(0.0,0.2,1.0,1.0)
    c.cd(2).SetPad(0.0,0.0,1.0,0.2)
    c.cd(1)

    # data	
    #print 'h_%s_%s_%s' % (entry.name,stream.data.name,region)
    data = rf.Get('h_%s_%s_%s' % (entry.name,stream.data.name,region)).Clone()
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.2)
    data.SetLabelSize(0.04,"y")
    data.Rebin(rb)
    if entry.xmax:
        data.SetAxisRange(0,entry.xmax)

    # stat. unc. on stack
    stacksum = data.Clone("tmp")
    stacksum.Reset()
    if syst:
        stacksum_UP = data.Clone("tmp")
        stacksum_UP.Reset()
        stacksum_DOWN = data.Clone("tmp")
        stacksum_DOWN.Reset()
    
    # Legend
    x0 = 0.57; y0 = 0.56; step = 0.045
    if entry.name == 'upsilon': 
        x0 = .62; y0 = 0.65
    if region == 'WCR_OS' and entry.name == 'met':
        x0 = .2; y0 = 0.6
    nl = 3 + len(stream.MC) - 1*entry.blind + noHands*1 
    l = ROOT.TLegend(x0,y0,x0+.25,y0+nl*step)

    if entry.leg:
        leg = entry.leg
        l = ROOT.TLegend(leg[0],leg[1],leg[2],leg[3])
    l.SetTextSize(0.032)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    if not entry.blind:
        l.AddEntry(data,'Data 2012 (%.0f)' % data.Integral(0,data.GetNbinsX()+1))

    # Stack
    st = ROOT.THStack(cname,"%s; %s; " % (cname,entry.xtitle))

    # Same Sign
    SS = hist.get_SS(rf,entry,stream.data,region,kSS = rQCD)
    SS.SetLineWidth(0)
    SS.SetLineColor(3)
    SS.SetFillColor(3)
    SS.SetMarkerStyle(0)
    SS.Rebin(rb)

    if syst:
        SS_UP = hist.get_SS(rf,entry,stream.data,region,kSS = rQCD_UP)
        SS_DOWN = hist.get_SS(rf,entry,stream.data,region,kSS = rQCD_DOWN)
        SS_UP.Rebin(rb); SS_DOWN.Rebin(rb)



    if entry.xmax:
        SS.SetAxisRange(0,entry.xmax)
    
    st.Add(SS)
    stacksum.Add(SS)
    if syst:
        stacksum_DOWN.Add(SS_DOWN)
        stacksum_UP.Add(SS_UP)

    if verbose:
        print 'Data: %.0f +/- %.0f' % (hist.IntegralError(data)[0],hist.IntegralError(data)[1])
        print 'NSS: %.0f +/- %.0f' % (hist.IntegralError(SS)[0],hist.IntegralError(SS)[1])
    if noHands:
        groups = stream.MC + [stream.signal]
    else:
        groups = stream.MC + [stream.LH,stream.RH]
        cname += '_LRstack'
        pname += '_LRstack'
    mc_leg = []
    for group in groups:
        kOS = ErrorFloat(1.0,0.0); kSS = ErrorFloat(1.0,0.0)
        kOS_UP = ErrorFloat(1.,0.); kOS_DOWN = ErrorFloat(1.,0.)
        kSS_UP = ErrorFloat(1.,0.); kSS_DOWN = ErrorFloat(1.,0.) 
        if group.name == 'Wlnu': 
            kOS = kW_OS; kSS = kW_SS * rQCD
            if syst:
                kOS_UP = kW_OS_UP; kOS_DOWN = kW_OS_DOWN
                kSS_UP = kW_SS_UP; kSS_DOWN = kW_SS_DOWN
        elif group.name == 'Ztautau_LH': kOS = kZL; kSS=kOS * rQCD
        elif group.name == 'Ztautau_RH': kOS = kZR; kSS=kOS * rQCD
        else: 
            if syst: 
                kSS_UP = kSS_UP * rQCD_UP
                kSS_DOWN = kSS_DOWN * rQCD_DOWN
            else:
                kSS = kSS * rQCD 
        h = hist.get_OS_SS(rf,entry,group,region,kOS=kOS,kSS=kSS)
        if syst:
            h_up = hist.get_OS_SS(sysf,entry,group, region+'_%s_UP' %syst, kOS = kOS_UP, kSS = kSS_UP)
            h_down = hist.get_OS_SS(sysf,entry,group, region+'_%s_DOWN' %syst, kOS = kOS_DOWN, kSS = kSS_DOWN)
            h_up.Rebin(rb)
            h_down.Rebin(rb)

        h.Rebin(rb)
        if entry.xmax:
            h.SetAxisRange(0,entry.xmax)
        h.SetLineColor(entry[group.name]['fill'])
        h.SetLineWidth(0)
        h.SetMarkerStyle(0)
        h.SetFillColor(entry[group.name]['fill'])
        if 'real' in group.name: h.SetFillStyle(3001)
        if 'fake' in group.name: h.SetFillStyle(3004)
        h.SetLabelSize(0.04,"y")
    
        st.Add(h)
        stacksum.Add(h)
        if syst:
            stacksum_UP.Add(h_up)
            stacksum_DOWN.Add(h_down)
        mc_leg.append( (h,group.legend,hist.IntegralError(h)[0],hist.IntegralError(h)[1]))
        #l.AddEntry(h,'%s (OS-SS): %.0f +/- %.0f' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1]))
        if verbose:
            print '%s: %.0f +/- %.0f' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1])

    for li in reversed(mc_leg):
        l.AddEntry(li[0], '%s (OS-SS): %.0f +/- %.0f' % (li[1],li[2],li[3]))

    l.AddEntry(SS,'SS: %.0f +/- %.0f' % (hist.IntegralError(SS)[0],hist.IntegralError(SS)[1]))

    data.SetXTitle(entry.xtitle)
    data.GetXaxis().SetTitleOffset(1.25)
    data.SetTitle('')

    stacksum.SetFillStyle(3354)
    stacksum.SetLineColor(1)
    stacksum.SetFillColor(12)
    stacksum.SetMarkerStyle(1)

    if syst:
        sys_error = hist.get_sys_error(stacksum, stacksum_DOWN, stacksum_UP)
        sys_error.SetFillColor(12)
        sys_error.SetLineWidth(2)
        sys_error.SetFillStyle(3019)


    maxi = max(data.GetMaximum(),st.GetMaximum())
    mini = 0.
    data.SetMinimum(1.2*mini)
    data.SetMaximum(1.5*maxi)

    if syst:
        if not entry.blind:
            data.Draw()
            st.Draw("HIST,SAME")
            sys_error.Draw("E2,SAME")
            data.Draw("SAME, PE")
        else:
            data.Reset()
            data.Draw()
            st.SetTitle('')
            
            st.SetMinimum(1.2*st.GetMinimum())
            st.SetMaximum(1.5*st.GetMaximum())
            st.Draw("HIST,SAME")
            sys_error.Draw("E2,SAME")

    else:
        if not entry.blind:
            data.Draw("PE")
            st.Draw("HIST,SAME")
            #st.Draw("HIST")
            stacksum.Draw("E2,SAME")
            data.Draw("SAME, PE")
        else:
            data.Reset()
            data.Draw()
            st.SetTitle('')
            
            st.SetMinimum(1.2*st.GetMinimum())
            st.SetMaximum(1.5*st.GetMaximum())
            st.Draw("HIST,SAME")
            stacksum.Draw("E2,SAME")


    c.RedrawAxis()
    l.Draw()

    if ratio:
        c.cd(2)
        ROOT.gPad.SetTopMargin(0)
        dataratio = data.Clone('ratio')
        if entry.xmax:
            dataratio.SetAxisRange(0,entry.xmax)
        dataratio.GetYaxis().SetTitle('Data / MC')
        dataratio.GetYaxis().SetTitleOffset(0.25)
        dataratio.GetYaxis().SetTitleSize(0.16)
        dataratio.SetLabelSize(.1,"y")
        dataratio.Divide(stacksum)
        dataratio.SetMaximum(1.5)
        dataratio.SetMinimum(0.5)
        
        err = stacksum.Clone('err')
        err.Divide(err)
        err.SetFillStyle(3354)
        err.SetLineColor(1)
        err.SetFillColor(12)
        err.SetMarkerStyle(1)
        xmin = dataratio.GetXaxis().GetXmin()
        xmax = entry.xmax or dataratio.GetXaxis().GetXmax()
        cl = ROOT.TLine(xmin,1,xmax,1)
        cl.SetLineStyle(3)
        dataratio.GetXaxis().SetLabelSize(0.04)
        dataratio.Draw('PE')
        cl.Draw("SAME")
        err.Draw("SAME,E2")
        c.Update()

    c.Update()
    c.Print('%s/%s.gif' % (PLOTDIR,cname))
    c.Print('%s/%s.gif+1' % (PLOTDIR,pname))

def plot_stat_qcd(spec,entry,region,test='',ratio=True,k_entry = upsilon_CR, verbose = False, SF = 'WQ',noHands=True, syst = None):
    stream = spec.stream
    PLOTDIR = spec.plotdir 
    if syst: PLOTDIR = PLOTDIR +'/syst'
    kspec = spec
    rf = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,spec.stream.name))
    kCR = ROOT.TFile('%s/%s_regions.root' % (kspec.OUTPUTDIR,kspec.stream.name))
    sysf = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,spec.stream.name))
    up = 'UP'
    down = 'DOWN'
    #up = up.lower(); down = down.lower()

    # for scaling left and right handed distributions
    kZL = ErrorFloat(1.,0.); kZR = ErrorFloat(1.,0)
    sys_var = False
    if not syst and region not in ['SR_OS','SR_SS']:
        sys_var = True

    # scale factors
    kW_OS = ErrorFloat(1.,0.); kW_OS_UP = ErrorFloat(1.,0.); kW_OS_DOWN = ErrorFloat(1.,0.) 
    kW_SS = ErrorFloat(1.,0.); kW_SS_UP = ErrorFloat(1.,0.); kW_SS_DOWN = ErrorFloat(1.,0.) 
    rQCD = ErrorFloat(1.,0.); rQCD_UP = ErrorFloat(1.,0.); rQCD_DOWN = ErrorFloat(1.,0.) 

    if 'W' in SF:
        kW_OS = hist.calc_kX(kCR,stream,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, kZL = kZL, kZR = kZR )
        kW_SS = hist.calc_kX(kCR,stream,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose)
        if syst or sys_var:
            kW_OS_UP = hist.calc_kX(rf,stream,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, kZL = kZL, kZR = kZR, sys = syst, sf = sysf, ud = up )
            kW_SS_UP = hist.calc_kX(rf,stream,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose, sys = syst, sf =sysf, ud = up)
            kW_OS_DOWN = hist.calc_kX(rf,stream,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, kZL = kZL, kZR = kZR , sys = syst, sf =sysf, ud = down)
            kW_SS_DOWN = hist.calc_kX(rf,stream,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose, sys = syst, sf = sysf, ud = down)
    if 'Q' in SF:
        rQCD = hist.calc_rQCD(kCR,stream,k_entry,'QCD',kW_OS = kW_OS, kW_SS = kW_SS, verbose = verbose)
        if syst or sys_var:
            rQCD_UP = hist.calc_rQCD(rf,stream,k_entry,'QCD',kW_OS = kW_OS_UP, kW_SS = kW_SS_UP, verbose = verbose, sys = syst, sf = sysf, ud = up)
            rQCD_DOWN = hist.calc_rQCD(rf,stream,k_entry,'QCD',kW_OS = kW_OS_DOWN, kW_SS = kW_SS_DOWN, verbose = verbose,  sys = syst, sf = sysf, ud = down)

    rb = entry.rebin

    # Canvas
    cname = '%s_%s_%s_qcd' % (stream.name,entry.name,region)
    pname = '%s_%s_qcd' % (stream.name,region.rstrip('_OS'))
    if syst:
        cname += '_%s' % syst
        pname += '_%s' % syst
    c = ROOT.TCanvas(cname,cname,1200,1000)
    if entry.blind: ratio = False

    #Split for ratio plot
    c.Divide(1,2)
    c.cd(1).SetPad(0.0,0.2,1.0,1.0)
    c.cd(2).SetPad(0.0,0.0,1.0,0.2)
    c.cd(1)

    # data  
    #print 'h_%s_%s_%s' % (entry.name,stream.data.name,region)
    if 'SR_OS' in region: dregion = 'SR_OS'
    if 'SR_SS' in region: dregion = 'SR_SS'
    data = rf.Get('h_%s_%s_%s' % (entry.name,stream.data.name,dregion)).Clone()
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.2)
    data.SetLabelSize(0.04,"y")
    data.Rebin(rb)
    if entry.xmax:
        data.SetAxisRange(0,entry.xmax)

    # stat. unc. on stack
    stacksum = data.Clone("tmp")
    stacksum.Reset()
    if syst:
        stacksum_UP = data.Clone("tmp")
        stacksum_UP.Reset()
        stacksum_DOWN = data.Clone("tmp")
        stacksum_DOWN.Reset()
    
    # Legend
    x0 = 0.72; y0 = 0.55; step = 0.045
    #if entry.name == 'upsilon': 
    #    x0 = .62; y0 = 0.5 + 0.1*entry.blind
    if region == 'WCR_OS' and entry.name == 'met':
        x0 = .2; y0 = 0.6

    nl = 3 + len(stream.MC) - 1*entry.blind + noHands*1 
    l = ROOT.TLegend(x0,y0,x0+.15,y0+nl*step)
    if entry.leg:
        leg = entry.leg
        l = ROOT.TLegend(leg[0],leg[1],leg[2],leg[3])

    l.SetTextSize(0.032)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    if not entry.blind:
        l.AddEntry(data,'Data 2012')# (%.0f)' % data.Integral(0,data.GetNbinsX()+1))

    tm = ROOT.TPaveText(0.18,0.85,0.22,0.90,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.04)
    if stream.name == 'mu':
        tm.AddText('#mu#tau_{had}')
    elif stream.name == 'el':
        tm.AddText('e#tau_{had}')
    #tm.Draw()

    # Stack
    st = ROOT.THStack(cname,"%s; %s; " % (cname,entry.xtitle))

    # QCD ESTIMATE
    #if region == 'SR_SS': rQCD = ErrorFloat(1.,0.)
    qcd = hist.get_QCD(rf,stream,entry,'SR_SS', kW = kW_SS, rqcd = rQCD)
    qcd.SetLineWidth(0)
    qcd.SetLineColor(3)
    qcd.SetFillColor(3)
    qcd.SetMarkerStyle(0)
    qcd.SetMarkerColor(3)
    qcd.Rebin(rb)

    if syst:
        #print rQCD_UP, rQCD_DOWN, kW_SS_UP, kW_SS_DOWN
        qcd_UP = hist.get_QCD(rf,stream,entry,'SR_SS_%s_%s' % (syst,up), kW = kW_SS_UP, rqcd = rQCD_UP, sysf = sysf)
        qcd_DOWN = hist.get_QCD(rf,stream,entry,'SR_SS_%s_%s' %(syst,down) ,kW = kW_SS_DOWN, rqcd = rQCD_DOWN, sysf = sysf)
        qcd_UP.Rebin(rb); qcd_DOWN.Rebin(rb)

    if entry.xmax:
        qcd.SetAxisRange(0,entry.xmax)
    
    st.Add(qcd)
    stacksum.Add(qcd)
    if syst:
        stacksum_DOWN.Add(qcd_DOWN)
        stacksum_UP.Add(qcd_UP)

    if verbose:
        total_est = ErrorFloat(*hist.IntegralError(qcd))
        print 'Data: %.0f +/- %.0f' % (hist.IntegralError(data)[0],hist.IntegralError(data)[1])
        print 'NQCD: %.0f +/- %.0f' % (hist.IntegralError(qcd)[0],hist.IntegralError(qcd)[1])
    if noHands:
        groups = stream.MC + [stream.signal]
    else:
        groups = stream.MC + [stream.LH,stream.RH]
        cname += '_LRstack'
        pname += '_LRstack'
    mc_leg = []

    n_total = ErrorFloat(*hist.IntegralError(hist.get_OS(rf,entry,stream.signal,region)))
    if syst:
        n_total_up = ErrorFloat(*hist.IntegralError(hist.get_OS(sysf,entry,stream.signal,region+'_%s_%s' %(syst,up))))  
        n_total_down = ErrorFloat(*hist.IntegralError(hist.get_OS(sysf,entry,stream.signal,region+'_%s_%s' %(syst,down))))  

    for group in groups:
        kOS = ErrorFloat(1.0,0.0);kSS = ErrorFloat(1.0,0.0);
        kOS_UP = ErrorFloat(1.,0.); kOS_DOWN = ErrorFloat(1.,0.)
        kSS_UP = ErrorFloat(1.,0.); kSS_DOWN = ErrorFloat(1.,0.)
        if group.name == 'Wlnu': 
            kOS = kW_OS;
            kSS = kW_SS
            if syst:
                kOS_UP = kW_OS_UP; kOS_DOWN = kW_OS_DOWN
                kSS_UP = kW_SS_UP; kSS_DOWN = kW_SS_DOWN
        #elif 'LH' in group.name:
        #    kOS = ErrorFloat(  0.5 *( 1- avg_p)); kOS_UP = kOS; kOS_DOWN = kOS
        #elif 'RH' in group.name:
        #    kOS = ErrorFloat(  0.5 *( 1+ avg_p)); kOS_UP = kOS; kOS_DOWN = kOS
        if 'SS' not in region:
            h = hist.get_OS(rf,entry,group,region,kOS=kOS)
            if syst:
                h_up = hist.get_OS(sysf,entry,group, region+'_%s_%s' %(syst,up), kOS = kOS_UP)
                h_down = hist.get_OS(sysf,entry,group, region+'_%s_%s' %(syst,down), kOS = kOS_DOWN)
                h_up.Rebin(rb)
                h_down.Rebin(rb)
        else:
            h = hist.get_SS(rf,entry,group,region,kSS=kSS)
            if syst:
                h_up = hist.get_SS(sysf,entry,group, region+'_%s_%s' %(syst,up), kSS = kSS_UP)
                h_down = hist.get_SS(sysf,entry,group, region+'_%s_%s' %(syst,down), kSS = kSS_DOWN)
                h_up.Rebin(rb)
                h_down.Rebin(rb)

        if 'L' in group.name:
            nleft = ErrorFloat(*hist.IntegralError(h))
            #eff_left = nleft / stream.Nleft
        elif 'R' in group.name:
            nright = ErrorFloat(*hist.IntegralError(h))
            #eff_right = nright / stream.Nright

        h.Rebin(rb)
        if entry.xmax:
            h.SetAxisRange(0,entry.xmax)
        h.SetLineColor(entry[group.name]['fill'])
        h.SetLineWidth(0)
        h.SetMarkerStyle(0)
        h.SetMarkerColor(entry[group.name]['fill'])
        h.SetFillColor(entry[group.name]['fill'])
        if 'real' in group.name: h.SetFillStyle(3001)
        if 'fake' in group.name: h.SetFillStyle(3004)
        h.SetLabelSize(0.04,"y")
        if 'L' in group.name: hl = h
        elif 'R' in group.name: hr = h
        else:   
            st.Add(h)
            stacksum.Add(h)
        if syst:
            stacksum_UP.Add(h_up)
            stacksum_DOWN.Add(h_down)
        mc_leg.append( (h,group.legend,hist.IntegralError(h)[0],hist.IntegralError(h)[1]))
        #l.AddEntry(h,'%s (OS-SS): %.0f +/- %.0f' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1]))
        if verbose:
            print '%s: %.0f +/- %.0f' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1])
            total_est += ErrorFloat(*hist.IntegralError(h))
            if syst:
                print '%s: %.0f +/- %.0f' % ('up',hist.IntegralError(h_up)[0],hist.IntegralError(h_up)[1])
                print '%s: %.0f +/- %.0f' % ('down',hist.IntegralError(h_down)[0],hist.IntegralError(h_down)[1])
        # if 'LH' in group.name and syst:
        #     left_up = ErrorFloat(*hist.IntegralError( h_up ))
        #     left_down = ErrorFloat(*hist.IntegralError( h_down ))
        #     eff_left_up = ErrorFloat(*hist.IntegralError( h_up )) / stream.Nleft
        #     eff_left_down = ErrorFloat(*hist.IntegralError( h_down )) / stream.Nleft
        #     #print row(group.name, ['%.0f +/- %.0f' % (hist.IntegralError(h_up)[0],hist.IntegralError(h_up)[1]),'%.0f +/- %.0f' % (hist.IntegralError(h_down)[0],hist.IntegralError(h_down)[1])], fcol=fcol,col=col)       
        #     #print row('' , [eff_left_up,eff_left_down], fcol =fcol, col=col )
        # if 'RH' in group.name and syst:
        #     right_up = ErrorFloat(*hist.IntegralError( h_up )) 
        #     right_down = ErrorFloat(*hist.IntegralError( h_down )) 
        #     eff_right_up = ErrorFloat(*hist.IntegralError( h_up )) / stream.Nright
        #     eff_right_down = ErrorFloat(*hist.IntegralError( h_down )) / stream.Nright

            #print row('' , [eff_right_up,eff_right_down], fcol =fcol, col=col )
    if verbose:
        print 'Total Estimated: %.0f +/- %.0f' % (total_est.val,total_est.err)

    #l.AddEntry(np0, 'Np0 Template'  )
    for li in reversed(mc_leg):
        #l.AddEntry(li[0], '%s (OS): %.0f +/- %.0f' % (li[1],li[2],li[3]))
        l.AddEntry(li[0], '%s'% (li[1]))

    #l.AddEntry(qcd,'QCD: %.0f +/- %.0f' % (hist.IntegralError(qcd)[0],hist.IntegralError(qcd)[1]))
    l.AddEntry(qcd,'Multijet')#: %.0f +/- %.0f' % (hist.IntegralError(qcd)[0],hist.IntegralError(qcd)[1]))


    #print row( syst, [ kW_OS_UP, kW_OS_DOWN, kW_SS_UP, kW_SS_DOWN, rQCD_UP, rQCD_DOWN, eff_left_up, eff_left_down, eff_right_up, eff_right_down, p_up, p_down], fmt = 'csv')
    if syst:
        eff_l_nom = ErrorFloat(0.1079,0.001389)
        eff_r_nom = ErrorFloat( 0.1129,0.001881)

        nl_up = ErrorFloat(.5725,0) * eff_l_nom / eff_left_up
        nr_up  = ErrorFloat(.4275,0) * eff_r_nom / eff_right_up
        p_up = (nr_up - nl_up)/(nr_up + nl_up)
        nl_down = ErrorFloat(.5725,0) * eff_l_nom / eff_left_down
        nr_down  = ErrorFloat(.4275,0) * eff_r_nom / eff_right_down
        p_down = (nr_down - nl_down)/(nr_down + nl_down)
        #p_down = (ErrorFloat(2,0)*n_total - left_down - right_down) / (right_down - left_down)


        #print row( syst, [ kW_OS_UP, kW_SS_UP, rQCD_UP, eff_left_up, eff_right_up, p_up], fmt = 'csv')
        #print row( syst, [ kW_OS_DOWN, kW_SS_DOWN, rQCD_DOWN, eff_left_down,eff_right_down,p_down], fmt = 'csv')
        #print row( syst, [p_up, p_down], fmt = 'twiki')
        if round(eff_l_nom.val,4) == round(eff_left_up.val,4) and round(eff_r_nom.val,4) == round(eff_right_up.val,4):
            pcup = 0
        else: 
            pcup = ( -0.145 - p_up.val ) / 0.145
        if round(eff_l_nom.val,4) == round(eff_left_down.val,4) and round(eff_r_nom.val,4) == round(eff_right_down.val,4):
            pcdown = 0
        else: 
            pcdown = (-0.145 - p_down.val) / 0.145
        print row( syst,[eff_left_up, eff_right_up, eff_left_down, eff_right_down,pcup,pcdown], fmt = 'twiki')

    elif not noHands:

        ptau = (ErrorFloat(2,0)*n_total- nleft - nright) / (nright - nleft)
        
        hl.Scale( 0.5 * (1-ptau.val))
        hr.Scale(0.5*(1+ptau.val))
        st.Add(hl); st.Add(hr)
        stacksum.Add(hl); stacksum.Add(hr)

        #print row( 'nominal', [kW_OS, kW_SS, rQCD, eff_left, eff_right, ptau], fcol = fcol, col=col   )

    data.SetXTitle(entry.xtitle)
    data.GetXaxis().SetTitleOffset(1.25)
    #data.SetTitle('')

    stacksum.SetFillStyle(3354)
    stacksum.SetLineColor(1)
    stacksum.SetFillColor(12)
    stacksum.SetMarkerStyle(1)

    if syst:
        sys_error = hist.get_sys_error(stacksum, stacksum_DOWN, stacksum_UP)
        sys_error.SetFillColor(1)
        sys_error.SetLineWidth(2)
        sys_error.SetFillStyle(0)
        #sys_error.SetFillStyle(3019)
        l2 = ROOT.TLegend( 0.25,0.75,0.4,0.85)
        l2.SetTextSize(0.032)
        l2.SetBorderSize(0)
        l2.SetFillColor(0)
        l2.AddEntry(stacksum_UP, '%s UP' % syst)
        l2.AddEntry(stacksum_DOWN, '%s DOWN' % syst)


    maxi = max(data.GetMaximum(),st.GetMaximum())
    mini = 0.
    data.SetMinimum(1.2*mini)
    data.SetMaximum(1.4*maxi)

    bw = data.GetXaxis().GetBinWidth(1)
    print bw
    data.GetYaxis().SetTitle('Events / %s %s' % (bw, entry.units))
    if syst:
        if not entry.blind:
            data.Draw()
            st.Draw("HIST,SAME")
            sys_error.Draw("E2,SAME")
            data.Draw("SAME, PE")
        else:
            data.Reset()
            data.Draw()
            #st.SetTitle('')
            
            st.SetMinimum(1.2*st.GetMinimum())
            st.SetMaximum(1.8*st.GetMaximum())
            st.Draw("HIST,SAME")
            stacksum_UP.Draw("HIST,SAME")
            stacksum_DOWN.SetLineStyle(2)
            stacksum_DOWN.Draw("HIST,SAME")
            #sys_error.Draw("E2,SAME")
            l2.Draw()
    else:
        if not entry.blind:
            data.Draw("PE")
            st.Draw("HIST,SAME")
            stacksum.Draw("E2,SAME")
            data.Draw("SAME, PE")
        else:
            data.Reset()
            data.Draw()
            st.SetTitle('')
            
            st.SetMinimum(1.2*st.GetMinimum())
            st.SetMaximum(1.5*st.GetMaximum())
            st.Draw("HIST,SAME")
            #np0.SetLineStyle(2)
            #np0.SetMarkerSize(0)
            #np0.Draw("HIST,SAME")
            #np0.Draw("P2,SAME")
            stacksum.Draw("E2,SAME")


    c.RedrawAxis()
    l.Draw()
    tm.Draw()

    if ratio:
        c.cd(2)
        ROOT.gPad.SetTopMargin(0)
        dataratio = data.Clone('ratio')
        if entry.xmax:
            dataratio.SetAxisRange(0,entry.xmax)
        dataratio.GetYaxis().SetTitle('Data / MC')
        dataratio.GetYaxis().SetTitleOffset(0.25)
        dataratio.GetYaxis().SetTitleSize(0.16)
        dataratio.SetLabelSize(.1,"y")
        dataratio.Divide(stacksum)
        dataratio.SetMaximum(1.5)
        dataratio.SetMinimum(0.5)
        
        err = stacksum.Clone('err')
        err.Divide(err)
        err.SetFillStyle(3354)
        err.SetLineColor(1)
        err.SetFillColor(12)
        err.SetMarkerStyle(1)
        xmin = dataratio.GetXaxis().GetXmin()
        xmax = entry.xmax or dataratio.GetXaxis().GetXmax()
        cl = ROOT.TLine(xmin,1,xmax,1)
        cl.SetLineStyle(3)
        dataratio.GetXaxis().SetLabelSize(0.04)
        dataratio.Draw('PE')
        cl.Draw("SAME")
        err.Draw("SAME,E2")
        c.Update()

    c.Update()
    #c.Print('%s/%s.gif' % (PLOTDIR,cname))
    c.Print('%s/%s.eps' % (PLOTDIR,cname))

def plot_QCD_sys(spec,entry,test='', k_entry = upsilon_CR, syst = None, SF = 'WQ', verbose=False):
    stream = spec.stream
    PLOTDIR = spec.plotdir 
    if syst: PLOTDIR = PLOTDIR + '/syst'
    region = 'SR_SS'
    rf = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,spec.stream.name))
    kCR = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,spec.stream.name))
    sysf = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,spec.stream.name))
    #sysf = ROOT.TFile('%s/%s_syst.root' % (spec.OUTPUTDIR,spec.stream.name))
    down = 'DOWN'
    up = 'UP'
    #down = down.lower(); up = up.lower()

    # for scaling left and right handed distributions
    kZL = ErrorFloat(1.,0.); kZR = ErrorFloat(1.,0)

    # scale factors
    kW_OS = ErrorFloat(1.,0.); kW_OS_UP = ErrorFloat(1.,0.); kW_OS_DOWN = ErrorFloat(1.,0.) 
    kW_SS = ErrorFloat(1.,0.); kW_SS_UP = ErrorFloat(1.,0.); kW_SS_DOWN = ErrorFloat(1.,0.) 
    rQCD = ErrorFloat(1.,0.); rQCD_UP = ErrorFloat(1.,0.); rQCD_DOWN = ErrorFloat(1.,0.) 
    if 'W' in SF:
        kW_OS = hist.calc_kX(kCR,stream,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, kZL = kZL, kZR = kZR )
        kW_SS = hist.calc_kX(kCR,stream,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose)
        if syst:
            kW_OS_UP = hist.calc_kX(rf,stream,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, kZL = kZL, kZR = kZR, sys = syst, sf = sysf, ud = up )
            kW_SS_UP = hist.calc_kX(rf,stream,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose, sys = syst, sf =sysf, ud = up)
            kW_OS_DOWN = hist.calc_kX(rf,stream,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, kZL = kZL, kZR = kZR , sys = syst, sf =sysf, ud = down)
            kW_SS_DOWN = hist.calc_kX(rf,stream,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose, sys = syst, sf = sysf, ud = down)
    if 'Q' in SF:
        rQCD = hist.calc_rQCD(kCR,stream,k_entry,'QCD',kW_OS = kW_OS, kW_SS = kW_SS, verbose = verbose)
        if syst:
            rQCD_UP = hist.calc_rQCD(rf,stream,k_entry,'QCD',kW_OS = kW_OS_UP, kW_SS = kW_SS_UP, verbose = verbose, sys = syst, sf = sysf, ud = up)
            rQCD_DOWN = hist.calc_rQCD(rf,stream,k_entry,'QCD',kW_OS = kW_OS_DOWN, kW_SS = kW_SS_DOWN, verbose = verbose,  sys = syst, sf = sysf, ud = down)

    if syst:
        # print 'kW_OS: %s' % kW_OS
        # print 'kW_OS_UP: %s' % kW_OS_UP
        # print 'kW_OS_DOWN: %s' % kW_OS_DOWN
        # print 'kW_SS: %s' % kW_SS
        # print 'kW_SS_UP: %s' % kW_SS_UP
        # print 'kW_SS_DOWN: %s' % kW_SS_DOWN
        # print 'rQCD: %s' % rQCD
        # print 'rQCD_UP: %s' % rQCD_UP
        # print 'rQCD_DOWN: %s' % rQCD_DOWN
        print row('',[kW_OS_UP,kW_OS_DOWN,kW_SS_UP,kW_SS_DOWN] )

    # rebin histos 
    rb = entry.rebin
    
    # Canvas
    cname = '%s_%s_%s_qcd' % (stream.name,entry.name,region)
    pname = '%s_%s_QCD' % (stream.name,region.rstrip('_OS'))
    if syst:
        cname += '_%s' % syst
        pname += '_%s' % syst
    c = ROOT.TCanvas(cname,cname,1200,1000)
    if entry.blind: ratio = False

    # data  
    #print 'h_%s_%s_%s' % (entry.name,stream.data.name,region)
    data = rf.Get('h_%s_%s_%s' % (entry.name,stream.data.name,region)).Clone()
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.2)
    data.SetLabelSize(0.04,"y")
    data.Rebin(rb)
    if entry.xmax:
        data.SetAxisRange(0,entry.xmax)
    
    # Legend
    x0 = 0.57; y0 = 0.6; step = 0.045
    if entry.name == 'upsilon': 
        x0 = .64; y0 = 0.55 + 0.1*entry.blind
    if region == 'WCR_OS' and entry.name == 'met':
        x0 = .2; y0 = 0.6
    nl = 3 + len(stream.MC) - 1*entry.blind #+ noHands*1 
    l = ROOT.TLegend(x0,y0,x0+.25,y0+nl*step)

    if entry.leg:
        leg = entry.leg
        l = ROOT.TLegend(leg[0],leg[1],leg[2],leg[3])
    l.SetTextSize(0.032)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    if not entry.blind:
        l.AddEntry(data,'Data 2012 (%.0f)' % data.Integral(0,data.GetNbinsX()+1))


    # QCD ESTIMATE
    rQCD = ErrorFloat(1.,0.)
    qcd = hist.get_QCD(rf,stream,entry,'SR_SS', kW = kW_SS, rqcd = rQCD)
    qcd.SetLineWidth(0)
    qcd.SetLineColor(3)
    qcd.SetFillColor(3)
    qcd.SetMarkerStyle(0)
    qcd.Rebin(rb)

    if syst:
        qcd_UP = hist.get_QCD(rf,stream,entry,'SR_SS_%s_%s' % (syst,up), kW = kW_SS_UP, rqcd = rQCD, sysf = sysf)
        qcd_DOWN = hist.get_QCD(rf,stream,entry,'SR_SS_%s_%s' %(syst,down) ,kW = kW_SS_DOWN, rqcd = rQCD, sysf = sysf)
        qcd_UP.Rebin(rb); qcd_DOWN.Rebin(rb)
    qcd.SetMaximum(1.5 * qcd.GetMaximum())

    #diff = hist.get_sys_error(qcd, qcd_DOWN, qcd_UP)
    #diff.SetFillColor(1)
    #diff.SetLineColor(1)
    #diff.SetLineWidth(2)
    #diff.SetFillStyle(0)

    qcd.SetXTitle(entry.xtitle)

    qcd.Draw("HIST")
    qcd_UP.Draw("HIST,SAME")
    qcd_DOWN.SetLineStyle(2)
    qcd_DOWN.Draw("HIST,SAME")

    l2 = ROOT.TLegend( 0.62,0.75,0.82,0.85)
    l2.SetTextSize(0.032)
    l2.SetBorderSize(0)
    l2.SetFillColor(0)
    l2.AddEntry(qcd_UP, '%s UP' % syst)
    l2.AddEntry(qcd_DOWN, '%s DOWN' % syst)
    l2.Draw()

    c.SaveAs('%s/%s.eps' % (PLOTDIR,pname))

def plot_CR(spec,entry,region,test='',ratio=True,k_entry = tau_et, verbose = False, SF='W'):

    stream = spec.stream
    PLOTDIR = spec.plotdir + '/control'
    rf = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,spec.stream.name))
    kCR = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,spec.stream.name))

    kW_OS = ErrorFloat(1.,0.) 
    kW_SS = ErrorFloat(1.,0.) 
    rQCD = ErrorFloat(1.,0.) 
    if 'W' in SF:
        kW_OS = hist.calc_kX(kCR,stream,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose)
        kW_SS = hist.calc_kX(kCR,stream,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose)
    if 'Q' in SF:
        rQCD = hist.calc_rQCD(kCR,stream,k_entry,'QCD',kW_OS = kW_OS, kW_SS = kW_SS, verbose = verbose)
    # rebin histos 
    rb = entry.rebin
    
    # Canvas
    cname = '%s_%s_%s' % (stream.name,entry.name,region)
    pname = '%s_%s' % (stream.name,region)
    c = ROOT.TCanvas(cname,cname,1200,1000)
    if entry.blind: ratio = False

    # Split for ratio plot
    c.Divide(1,2)
    c.cd(1).SetPad(0.0,0.2,1.0,1.0)
    c.cd(2).SetPad(0.0,0.0,1.0,0.2)
    c.cd(2).SetTopMargin(0)
    c.cd(1)

    # data  
    data = rf.Get('h_%s_%s_%s' % (entry.name,stream.data.name,region)).Clone()
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.2)
    data.SetLabelSize(0.04,"y")
    data.Rebin(rb)
    if entry.xmax:
        data.SetAxisRange(0,entry.xmax)

    # stat. unc. on stack
    stacksum = data.Clone("tmp")
    stacksum.Reset()
    
    # Legend
    x0 = 0.6; y0 = 0.6; step = 0.05
    if entry.name in ['met','upsilon']: 
        x0 = .2; y0 = 0.6
    nl = 3 + len(stream.MC) 
    l = ROOT.TLegend(x0,y0,x0+.25,y0+nl*step)


    if entry.leg:
        leg = entry.leg
        l = ROOT.TLegend(leg[0],leg[1],leg[2],leg[3])
    l.SetTextSize(0.032)
    l.SetBorderSize(0)
    l.SetFillColor(0)

    # Stack
    st = ROOT.THStack(cname,"%s; %s; " % (cname,entry.xtitle))
    # Same Sign

    if verbose:
        print 'Data: %.0f +/- %.0f' % (hist.IntegralError(data)[0],hist.IntegralError(data)[1])

    l.AddEntry(data,'Data: %.0f +/- %.0f' % (hist.IntegralError(data)[0],hist.IntegralError(data)[1]))

    groups = stream.MC + [stream.signal]
    mc_leg = []
    for group in groups:
        kOS = ErrorFloat(1.0,0.0); kSS = ErrorFloat(1.0,0.0) 
        if group.name == 'Wlnu': kOS = kW_OS; kSS = kW_SS 
        if 'OS' in region:
            h = hist.get_OS(rf,entry,group,region,kOS=kOS)
        elif 'SS' in region:
            h = hist.get_SS(rf,entry,group,region,kSS=kSS)
        h.Rebin(rb)
        if entry.xmax:
            h.SetAxisRange(0,entry.xmax)
        h.SetLineColor(entry[group.name]['fill'])
        h.SetLineWidth(0)
        h.SetMarkerStyle(0)
        h.SetFillColor(entry[group.name]['fill'])
        h.SetLabelSize(0.04,"y")
    
        st.Add(h)
        stacksum.Add(h)
        mc_leg.append( (h,group.legend,hist.IntegralError(h)[0],hist.IntegralError(h)[1]))
        #l.AddEntry(h,'%s: %.0f +/- %.0f' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1]))
        if verbose:
            print '%s: %.0f +/- %.0f' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1])

    for li in reversed(mc_leg):
        l.AddEntry(li[0], '%s: %.0f +/- %.0f' % (li[1],li[2],li[3]))

    data.SetXTitle(entry.xtitle)
    data.SetTitle('')
    data.GetXaxis().SetTitleOffset(1.25)

    stacksum.SetFillStyle(3354)
    stacksum.SetLineColor(1)
    stacksum.SetFillColor(12)
    stacksum.SetMarkerStyle(1)

    maxi = max(data.GetMaximum(),st.GetMaximum())
    mini = 0.
    data.SetMinimum(1.2*mini)
    data.SetMaximum(1.5*maxi)

    data.Draw()
    st.Draw("HIST,SAME")
    stacksum.Draw("E2,SAME")
    data.Draw("SAME, PE")

    c.RedrawAxis()
    l.Draw()

    if ratio:
        c.cd(2)

        dataratio = data.Clone('ratio')
        if entry.xmax:
            dataratio.SetAxisRange(0,entry.xmax)
        dataratio.GetYaxis().SetTitle('Data / MC')
        dataratio.GetYaxis().SetTitleOffset(0.25)
        dataratio.GetYaxis().SetTitleSize(0.16)
        dataratio.Divide(stacksum)
        dataratio.SetMaximum(1.5)
        dataratio.SetMinimum(0.5)
        
        err = stacksum.Clone('err')
        err.Divide(err)
        err.SetFillStyle(3354)
        err.SetLineColor(1)
        err.SetFillColor(12)
        err.SetMarkerStyle(1)
        xmin = dataratio.GetXaxis().GetXmin()
        xmax = entry.xmax or dataratio.GetXaxis().GetXmax()
        cl = ROOT.TLine(xmin,1,xmax,1)
        cl.SetLineStyle(3)
        dataratio.GetXaxis().SetLabelSize(0.04)
        dataratio.Draw('PE')
        cl.Draw("SAME")
        err.Draw("SAME,E2")
        c.Update()


    c.Update()
    c.Print('%s/%s.gif' % (PLOTDIR,cname))
    c.Print('%s/%s.gif+1' % (PLOTDIR,pname))

def plot_modes(spec,entry,region,test='',ratio=True,k_entry = tau_et, verbose = False, SF = 'WQ', hand = False):
    stream = spec.stream
    PLOTDIR = spec.plotdir + '/modes'
    rf = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,spec.stream.name))
    kCR = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,spec.stream.name))
    kZL = ErrorFloat(1.,0.); kZR = ErrorFloat(1.,0)

    # scale factors
    kW_OS = ErrorFloat(1.,0.) 
    kW_SS = ErrorFloat(1.,0.) 
    rQCD = ErrorFloat(1.,0.) 
    if 'W' in SF:
        kW_OS = hist.calc_kX(kCR,stream,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, kZL = kZL, kZR = kZR )
        kW_SS = hist.calc_kX(kCR,stream,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose)
    if 'Q' in SF:
        rQCD = hist.calc_rQCD(kCR,stream,k_entry,'QCD',kW_OS = kW_OS, kW_SS = kW_SS, verbose = verbose)

    # rebin histos 
    rb = entry.rebin
    
    # Canvas
    cname = '%s_%s_%s' % (stream.name,entry.name,region)
    pname = '%s_%s_modes' % (stream.name,region.rstrip('_OS'))
    if hand: pname += '_%s' % hand
    c = ROOT.TCanvas(cname,cname,1200,1000)
    if entry.blind: ratio = False

    # Split for ratio plot
    c.Divide(1,2)
    c.cd(1).SetPad(0.0,0.2,1.0,1.0)
    c.cd(2).SetPad(0.0,0.0,1.0,0.2)
    c.cd(1)

    # data  
    #print 'h_%s_%s_%s' % (entry.name,stream.data.name,region)
    data = rf.Get('h_%s_%s_%s' % (entry.name,stream.data.name,region)).Clone()
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.2)
    data.SetLabelSize(0.04,"y")
    data.Rebin(rb)
    if entry.xmax:
        data.SetAxisRange(0,entry.xmax)

    # stat. unc. on stack
    stacksum = data.Clone("tmp")
    stacksum.Reset()
    
    # Legend
    x0 = 0.6; y0 = 0.5; step = 0.04
    if entry.name == 'upsilon': 
        x0 = .60; y0 = 0.60
    nl = 6 + len(stream.MC) - 1*entry.blind
    if hand: 
        nl = nl - len(stream.MC) - 1
        step = 0.05
        y0 += .1
    l = ROOT.TLegend(x0,y0,x0+.25,y0+nl*step)

    if entry.leg:
        leg = entry.leg
        l = ROOT.TLegend(leg[0],leg[1],leg[2],leg[3])
    l.SetTextSize(0.032)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    if not entry.blind:
        l.AddEntry(data,'Data 2012 (%.0f)' % data.Integral(0,data.GetNbinsX()+1))

    # Stack
    st = ROOT.THStack(cname,"%s; %s; " % (cname,entry.xtitle))
    # Same Sign
    SS = hist.get_SS(rf,entry,stream.data,region,kSS = rQCD)
    SS.SetLineWidth(0)
    SS.SetLineColor(3)
    SS.SetFillColor(3)
    SS.SetMarkerStyle(0)
    SS.Rebin(rb)
    if entry.xmax:
        SS.SetAxisRange(0,entry.xmax)
    if not hand:
        st.Add(SS)
        stacksum.Add(SS)

    if verbose:
        print 'Data: %.0f +/- %.0f' % (hist.IntegralError(data)[0],hist.IntegralError(data)[1])
        print 'NSS: %.0f +/- %.0f' % (hist.IntegralError(SS)[0],hist.IntegralError(SS)[1])


    if hand == 'LH': signal = stream.LH
    elif hand == 'RH': signal = stream.RH
    else: signal = stream.signal 
    groups = stream.MC + [signal]

    mc_leg = []
    for group in groups:
        kOS = ErrorFloat(1.0,0.0); kSS = ErrorFloat(1.0,0.0) 
        if group.name == 'Wlnu': kOS = kW_OS; kSS = kW_SS * rQCD
        elif group.name == 'Ztautau_LH': kOS = kZL; kSS=kOS * rQCD
        elif group.name == 'Ztautau_RH': kOS = kZR; kSS=kOS * rQCD
        else: kSS = kSS * rQCD 

        if group.name == signal.name:
            ## charged pions
            for mr in ['Xp','1pXn','1p0n','1p1n']:
                mregion = mr + '_' + region
                hm = hist.get_OS_SS(rf,entry,group,mregion,kOS=kOS,kSS=kSS).Clone()
                hm.SetFillColor(entry[group.name]['fill'])
                hm.SetFillStyle(modefill[mr])
                hm.SetLineColor(entry[group.name]['fill'])
                hm.SetLineWidth(0)
                hm.SetMarkerStyle(0)
                hm.Rebin(rb)
                st.Add(hm)
                stacksum.Add(hm)
                mc_leg.append( (hm,mr,hist.IntegralError(hm)[0],hist.IntegralError(hm)[1]))
        else:
            if hand: continue
            h = hist.get_OS_SS(rf,entry,group,region,kOS=kOS,kSS=kSS)

            h.Rebin(rb)
            if entry.xmax:
                h.SetAxisRange(0,entry.xmax)
            h.SetLineColor(entry[group.name]['fill'])
            h.SetLineWidth(0)
            h.SetMarkerStyle(0)
            h.SetFillColor(entry[group.name]['fill'])
            h.SetLabelSize(0.04,"y")
            
            st.Add(h)
            stacksum.Add(h)
            mc_leg.append( (h,group.legend,hist.IntegralError(h)[0],hist.IntegralError(h)[1]))
        #l.AddEntry(h,'%s (OS-SS): %.0f +/- %.0f' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1]))
        if verbose:
            print '%s: %.0f +/- %.0f' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1])

    for li in reversed(mc_leg):
        l.AddEntry(li[0], '%s (OS-SS): %.0f +/- %.0f' % (li[1],li[2],li[3]))
    if not hand:
        l.AddEntry(SS,'SS: %.0f +/- %.0f' % (hist.IntegralError(SS)[0],hist.IntegralError(SS)[1]))

    data.SetXTitle(entry.xtitle)
    data.SetTitle('')
    data.GetXaxis().SetTitleOffset(1.25)

    stacksum.SetFillStyle(3354)
    stacksum.SetLineColor(1)
    stacksum.SetFillColor(12)
    stacksum.SetMarkerStyle(1)


    maxi = max(data.GetMaximum(),st.GetMaximum())
    if hand: maxi = st.GetMaximum()
    mini = 0.
    data.SetMinimum(1.2*mini)
    data.SetMaximum(1.5*maxi)

    if not entry.blind and not hand:
        data.Draw()
        st.Draw("HIST,SAME")
        stacksum.Draw("E2,SAME")
        data.Draw("SAME, PE")
    else:
        data.Reset()
        data.Draw()
        st.SetTitle('')
        
        st.SetMinimum(1.2*st.GetMinimum())
        st.SetMaximum(1.5*st.GetMaximum())
        st.Draw("HIST,SAME")
        stacksum.Draw("E2,SAME")

    c.RedrawAxis()
    l.Draw()

    if ratio:
        c.cd(2)

        dataratio = data.Clone('ratio')
        if entry.xmax:
            dataratio.SetAxisRange(0,entry.xmax)
        dataratio.GetYaxis().SetTitle('Data / MC')
        dataratio.GetYaxis().SetTitleOffset(0.25)
        dataratio.GetYaxis().SetTitleSize(0.16)
        dataratio.Divide(stacksum)
        dataratio.SetMaximum(1.5)
        dataratio.SetMinimum(0.5)
        
        err = stacksum.Clone('err')
        err.Divide(err)
        err.SetFillStyle(3354)
        err.SetLineColor(1)
        err.SetFillColor(12)
        err.SetMarkerStyle(1)
        xmin = dataratio.GetXaxis().GetXmin()
        xmax = entry.xmax or dataratio.GetXaxis().GetXmax()
        cl = ROOT.TLine(xmin,1,xmax,1)
        cl.SetLineStyle(3)
        dataratio.GetXaxis().SetLabelSize(0.04)
        dataratio.Draw('PE')
        cl.Draw("SAME")
        err.Draw("SAME,E2")
        c.Update()


    #c.Print('%s/%s.gif' % (PLOTDIR,cname))
    c.Print('%s/%s.gif+1' % (PLOTDIR,pname))

def plot_LR(spec,entry,region,test='',verbose = False, hand = 'LR',  norm = False, ratio = False, comp = None):
    stream = spec.stream
    PLOTDIR = spec.plotdir 
    rf = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,spec.stream.name))
    # Canvas
    cname = '%s_%s_%s_%s' % (stream.name,entry.name,region,hand)
    pname = '%s_%s_%s' % (stream.name,region.rstrip('_OS'),hand)
    if norm: 
        cname+='_Norm'
        pname+='_Norm'

    if comp:
        cf = ROOT.TFile('%s/%s_regions.root' % (comp.OUTPUTDIR,comp.stream.name))
        cname += ''

    c = ROOT.TCanvas(cname,cname,1200,1000)
    #c.SetLogy()

    rb = entry.rebin
    # Split for ratio plot
    # c.Divide(1,2)
    # c.cd(1).SetPad(0.0,0.2,1.0,1.0)
    # c.cd(2).SetPad(0.0,0.0,1.0,0.2)
    # c.cd(2).SetTopMargin(0)
    # c.cd(1)
    # data  
    #rf.ls()
    #sys.exit()
    total = rf.Get('h_%s_%s_%s' % (entry.name, stream.signal.name, region))
    if 'L' in hand:
        left = rf.Get('h_%s_%s_%s' % (entry.name, stream.LH.name, region))
    if 'R' in hand:
        right = rf.Get('h_%s_%s_%s' % (entry.name, stream.RH.name, region))
    if 'O' in hand:
        #zoff = rf.Get('h_%s_%s_%s' % (entry.name, stream.Zoff.name, region))
        zlo = rf.Get('h_%s_%s_%s' % (entry.name, stream.Zlo.name, region))
        zhi = rf.Get('h_%s_%s_%s' % (entry.name, stream.Zhi.name, region))
    if comp:
        ctotal = cf.Get('h_%s_%s_%s' % (entry.name, comp.stream.signal.name, region))
        cleft = cf.Get('h_%s_%s_%s' % (entry.name, comp.stream.LH.name, region))
        cright = cf.Get('h_%s_%s_%s' % (entry.name, comp.stream.RH.name, region))

    #nleft = left.Integral(); nright = right.Integral(); ntotal = total.Integral()
    #pol = (2*ntotal - nright - nleft)/(nright-nleft)
    if comp:
        ncleft = cleft.Integral(); ncright = cright.Integral(); nctotal = ctotal.Integral()
        cpol = (2*nctotal - ncright - ncleft)/(ncright-ncleft)


    # if norm:
    #     left.Scale( 1/ left.Integral())
    #     right.Scale(1/right.Integral())
    #     if comp:
    #         cleft.Scale(1./cleft.Integral())
    #         cright.Scale(1./cright.Integral())
    # else:
    #     left.Scale( 0.5 * (1-pol))
    #     right.Scale( 0.5 * (1+pol))



    #for (h,group) in [(total, stream.signal),(zoff,stream.Zoff)]:

    #for (h,group) in [(total, stream.signal),(zhi,stream.Zhi),(zlo,stream.Zlo)]:
    for (h,group) in [ (left,stream.LH), (right,stream.RH), (total, stream.signal)]:
        entry.add_group(group)
        h.Rebin(rb)
        h.SetLineColor(entry[group.name]['fill'])
        h.SetLineWidth(3)
        h.SetMarkerSize(0)
        h.SetLineStyle(1)
        h.SetMarkerColor(entry[group.name]['fill'])
        h.SetFillColor(0)
        h.SetLabelSize(0.04,"y")
    if comp:
        for (h,group) in [ (cleft,comp.stream.LH), (cright,comp.stream.RH), (ctotal, comp.stream.signal)]:
            entry.add_group(group)
            h.SetLineColor(entry[group.name]['fill'])
            h.SetLineWidth(1)
            h.SetMarkerSize(0)
            h.SetLineStyle(2)
            h.SetMarkerColor(entry[group.name]['fill'])
            h.SetFillColor(0)
            h.SetLabelSize(0.04,"y")

    leg = ROOT.TLegend( 0.7,0.75,0.90,0.92)
    #leg = ROOT.TLegend( 0.22,0.75,0.42,0.85)
    leg.SetTextSize(0.034)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)

    draw = 'HIST'
    if 'L' in hand:
        left.Draw(draw)
        left.Draw("P2,SAME")
        left.SetMaximum(1.3*left.GetMaximum())
        leg.AddEntry(left, stream.LH.legend)
        draw+=',SAME'
        if comp: 
            cleft.Draw(draw)
            leg.AddEntry(cleft, 'Alpgen Np0')
            cleft.Draw("P2,SAME")
    if 'R' in hand: 
        right.SetMaximum(1.3*right.GetMaximum())
        right.Draw(draw)
        right.Draw("P2,SAME")
        leg.AddEntry(right, stream.RH.legend)
        if comp: 
            cright.Draw(draw+',SAME')
            cright.Draw("P2,SAME")
            leg.AddEntry(cright, 'Alpgen Np0')

    if 'T' in hand:
        total.SetMaximum(1.3*total.GetMaximum())
        total.SetXTitle( entry.xtitle   )
        total.Draw(draw)
        total.Draw("P2,SAME")
        print 'On peak: %f' % total.Integral()
        leg.AddEntry(total,stream.signal.legend)
        draw += ',SAME'

    if 'O' in hand:
        #zoff.Draw(draw)
        #zoff.Draw("P2,SAME")
        #print 'Off peak: %s' % zoff.Integral()
        #zhi.Draw(draw)
        zhi.Draw("HIST,SAME")
        zhi.Draw("P2,SAME")
        zlo.Draw("HIST,SAME")
        zlo.Draw("P2,SAME")
        #print 'Zlo: %f' % zlo.Integral()
        #print 'Zhi: %f' % zhi.Integral()
        #zlo.Draw("P2,SAME")
        leg.AddEntry(zlo,stream.Zlo.legend)
        leg.AddEntry(zhi,stream.Zhi.legend)
        #leg.AddEntry(zoff,stream.Zoff.legend)

    leg.Draw()



    c.Print('%s/%s.ps' % (PLOTDIR,cname))

def compare_QCD_shapes(spec,entry,k_entry = tau_et, verbose = False):
    # compare QCD shapes with EW subtracted
    stream = spec.stream
    PLOTDIR = spec.plotdir + '/control'
    rf = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,spec.stream.name))

    # rebin histos 
    rb = entry.rebin
    kZL = ErrorFloat(1.,0.); kZR = ErrorFloat(1.,0)
    kW_OS = hist.calc_kX(rf,stream,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, kZL = kZL, kZR = kZR )
    kW_SS = hist.calc_kX(rf,stream,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose)

    # Canvas
    cname = '%s_%s_%s' % (stream.name,entry.name,'QCD_shapes')
    pname = '%s_%s' % (stream.name,'QCD_shapes')
    c = ROOT.TCanvas(cname,cname,1200,1000)
    if entry.blind: ratio = False

    # Split for ratio plot
    c.Divide(1,2)
    c.cd(1).SetPad(0.0,0.2,1.0,1.0)
    c.cd(2).SetPad(0.0,0.0,1.0,0.2)
    c.cd(2).SetTopMargin(0)
    c.cd(1)

    # Legend
    l = ROOT.TLegend(0.7,0.7,0.90,0.90)
    l.SetTextSize(0.032)
    l.SetBorderSize(0)
    l.SetFillColor(0)

    qcd_os = hist.get_QCD(rf,stream,entry,'QCD_OS', kW = kW_OS)
    qcd_ss = hist.get_QCD(rf,stream,entry,'QCD_SS', kW = kW_SS)
    sr_ss = hist.get_QCD(rf,stream,entry,'SR_SS', kW = kW_SS)
    sr_ss.SetLineColor(3); sr_ss.SetMarkerColor(3)
    qcd_ss.SetLineColor(2); qcd_ss.SetMarkerColor(2)
    for h,name in [ (qcd_os, 'QCD OS'), (qcd_ss, 'QCD SS'), (sr_ss, 'SR_SS') ]:
        h.Rebin(rb)
        h.Scale(1./ h.Integral(0,h.GetNbinsX()+1))
        l.AddEntry(h,name)
    maxi = max(qcd_os.GetMaximum(),sr_ss.GetMaximum(),qcd_ss.GetMaximum())
    qcd_os.SetMaximum(maxi*1.2)
    qcd_os.SetXTitle(entry.xtitle)
    qcd_os.GetXaxis().SetTitleOffset(1.25)
    qcd_os.SetTitle('')

    qcd_os.Draw("C")
    qcd_ss.Draw("C,SAME")
    sr_ss.Draw("C,SAME")

    l.Draw()

    c.Update()
    c.Print('%s/%s.gif' % (PLOTDIR,cname))
    c.Print('%s/%s.gif+1' % (PLOTDIR,pname))

def plot_SS(spec,entry,region,k_entry = tau_et, verbose = False):
    stream = spec.stream
    PLOTDIR = spec.plotdir + '/control'
    rf = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,spec.stream.name))
    kCR = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,spec.stream.name))
       
    # for scaling left and right handed distributions
    kZL = ErrorFloat(1.,0.); kZR = ErrorFloat(1.,0)

    # scale factors
    kW_SS = hist.calc_kX(kCR,stream,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose)

    #  histos 
    rb = entry.rebin
    
    # Canvas
    cname = '%s_%s_%s' % (stream.name,entry.name,'SS_QCD')
    pname = '%s_%s' % (stream.name,region)
    c = ROOT.TCanvas(cname,cname,1200,1000)

    # Split for ratio plot
    c.Divide(1,2)
    c.cd(1).SetPad(0.0,0.2,1.0,1.0)
    c.cd(2).SetPad(0.0,0.0,1.0,0.2)
    c.cd(2).SetTopMargin(0)
    c.cd(1)

    # data  
    data = rf.Get('h_%s_%s_%s' % (entry.name,stream.data.name,region)).Clone()
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.2)
    data.SetLabelSize(0.04,"y")
    data.Rebin(rb)
    if entry.xmax:
        data.SetAxisRange(0,entry.xmax)

    # stat. unc. on stack
    stacksum = data.Clone("tmp")
    stacksum.Reset()
    
    # Legend
    x0 = 0.6; y0 = 0.6; step = 0.05
    if entry.name in ['upsilon']: 
        x0 = .65; y0 = 0.6
    nl = 3 + len(stream.MC) 
    l = ROOT.TLegend(x0,y0,x0+.25,y0+nl*step)
    if entry.leg:
        leg = entry.leg
        l = ROOT.TLegend(leg[0],leg[1],leg[2],leg[3])
    l.SetTextSize(0.032)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    if not entry.blind:
        l.AddEntry(data,'SS Data (%.0f)' % data.Integral(0,data.GetNbinsX()+1))

    # Stack
    st = ROOT.THStack(cname,"%s; %s; " % (cname,entry.xtitle))
    # qcd estimate
    qcd = hist.get_QCD(rf,stream,entry,'SR_SS', kW = kW_SS)
    qcd.SetLineWidth(0)
    qcd.SetLineColor(3)
    qcd.SetFillColor(3)
    qcd.SetMarkerStyle(0)
    qcd.Rebin(rb)
    if entry.xmax:
        qcd.SetAxisRange(0,entry.xmax)
    
    st.Add(qcd)
    stacksum.Add(qcd)

    if verbose:
        print 'Data: %.0f +/- %.0f' % (hist.IntegralError(data)[0],hist.IntegralError(data)[1])
        print 'NQCD: %.0f +/- %.0f' % (hist.IntegralError(qcd)[0],hist.IntegralError(qcd)[1])

    groups = stream.MC + [stream.signal]
    mc_leg = []
    for group in groups:
        kSS = ErrorFloat(1.0,0.0) 
        if group.name == 'Wlnu': kSS = kW_SS 
        h = hist.get_SS(rf,entry,group,region,kSS=kSS)

        h.Rebin(rb)
        if entry.xmax:
            h.SetAxisRange(0,entry.xmax)
        h.SetLineColor(entry[group.name]['fill'])
        h.SetLineWidth(0)
        h.SetMarkerStyle(0)
        h.SetFillColor(entry[group.name]['fill'])
        if 'real' in group.name: h.SetFillStyle(3001)
        if 'fake' in group.name: h.SetFillStyle(3004)
        h.SetLabelSize(0.04,"y")
    
        st.Add(h)
        stacksum.Add(h)
        mc_leg.append( (h,group.legend,hist.IntegralError(h)[0],hist.IntegralError(h)[1]))
        #l.AddEntry(h,'%s (OS-SS): %.0f +/- %.0f' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1]))
        if verbose:
            print '%s: %.0f +/- %.0f' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1])

    for li in reversed(mc_leg):
        l.AddEntry(li[0], '%s: %.0f +/- %.0f' % (li[1],li[2],li[3]))

    l.AddEntry(qcd,'QCD: %.0f +/- %.0f' % (hist.IntegralError(qcd)[0],hist.IntegralError(qcd)[1]))

    data.SetXTitle(entry.xtitle)
    data.SetTitle('')

    stacksum.SetFillStyle(3354)
    stacksum.SetLineColor(1)
    stacksum.SetFillColor(12)
    stacksum.SetMarkerStyle(1)


    maxi = max(data.GetMaximum(),st.GetMaximum())
    mini = 0.
    data.SetMinimum(1.2*mini)
    data.SetMaximum(1.5*maxi)

    data.Draw()
    st.Draw("HIST,SAME")
    stacksum.Draw("E2,SAME")
    data.Draw("SAME, PE")

    c.RedrawAxis()
    l.Draw()

    c.Update()
    c.Print('%s/%s.gif' % (PLOTDIR,cname))
    c.Print('%s/%s.gif+1' % (PLOTDIR,pname))

def plot_Np(spec,entry,region,hand='L'):

    stream = spec.stream
    PLOTDIR = spec.plotdir 
    rf = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,spec.stream.name))
    # rebin histos 
    rb = entry.rebin
    # Canvas
    cname = '%s_%s' % (stream.name,hand)
    c = ROOT.TCanvas(cname,cname,1500,1000)

    # Legend
    l = ROOT.TLegend(0.6,0.75,0.8,0.90)
    l.SetTextSize(0.028)
    l.SetBorderSize(0)
    l.SetFillColor(0)

    mc_leg = []

    groups = spec.groups
    draw = ''
    maxi = 0
    for group in groups:
        if 'Left' in group.name and 'L' not in hand: continue
        if 'Right' in group.name and 'R' not in hand: continue
        h = hist.get_OS(rf,entry,group,region)
        #Normalize
        h.Rebin(rb*2)
        h.Scale(1/h.Integral(0,h.GetNbinsX()+1))
        h.SetLineColor(entry[group.name]['fill'])
        if '0' in group.name or '5' in group.name:
            h.SetLineWidth(3)
        else: h.SetLineWidth(1)        
        h.SetMarkerStyle(0)
        h.SetFillColor(0)
        h.SetLabelSize(0.04,"y")
        h.SetXTitle(entry.xtitle)
        h.GetXaxis().SetTitleOffset(1.25)
        h.SetTitle('')
        maxi = max(maxi,h.GetMaximum())
        h.SetMaximum(1.4*maxi)
        #mc_leg.append(h,group.name)

        h.Draw(draw)
        l.AddEntry(h,group.name)
        draw += ',SAME'
    l.Draw()
    c.Print('%s/%s.gif+1' % (PLOTDIR,cname) )

def plot_toy(spec, SF='WQ'):
    stream = spec.stream
    PLOTDIR = spec.plotdir 
    rf = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,spec.stream.name))
    kCR = ROOT.TFile('%s/%s_regions.root' % (spec.OUTPUTDIR,spec.stream.name))
    # for scaling left and right handed distributions
    kZL = ErrorFloat(1.,0.); kZR = ErrorFloat(1.,0)
    # scale factors
    kW_OS = ErrorFloat(1.,0.); kW_OS_UP = ErrorFloat(1.,0.); kW_OS_DOWN = ErrorFloat(1.,0.) 
    kW_SS = ErrorFloat(1.,0.); kW_SS_UP = ErrorFloat(1.,0.); kW_SS_DOWN = ErrorFloat(1.,0.) 
    rQCD = ErrorFloat(1.,0.); rQCD_UP = ErrorFloat(1.,0.); rQCD_DOWN = ErrorFloat(1.,0.) 
    if 'W' in SF:
        kW_OS = hist.calc_kX(kCR,stream,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, kZL = kZL, kZR = kZR )
        kW_SS = hist.calc_kX(kCR,stream,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose)
    if 'Q' in SF:
        rQCD = hist.calc_rQCD(kCR,stream,k_entry,'QCD',kW_OS = kW_OS, kW_SS = kW_SS, verbose = verbose)
    # rebin histos 
    # rb = entry.rebin
    
    # Canvas
    cname = '%s_%s_%s' % (stream.name,entry.name,region)
    pname = '%s_%s' % (stream.name,region.rstrip('_OS'))
    if syst:
        cname += '_%s' % syst
        pname += '_%s' % syst
    c = ROOT.TCanvas(cname,cname,1200,1000)
    if entry.blind: ratio = False

    # Split for ratio plot
    c.Divide(1,2)
    c.cd(1).SetPad(0.0,0.2,1.0,1.0)
    c.cd(2).SetPad(0.0,0.0,1.0,0.2)
    c.cd(1)

    # toy data  
    wsf = ROOT.TFile("")
    data = rf.Get('h_%s_%s_%s' % (entry.name,stream.data.name,region)).Clone()
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.2)
    data.SetLabelSize(0.04,"y")
    data.Rebin(rb)
    if entry.xmax:
        data.SetAxisRange(0,entry.xmax)

    # stat. unc. on stack
    stacksum = data.Clone("tmp")
    stacksum.Reset()
    if syst:
        stacksum_UP = data.Clone("tmp")
        stacksum_UP.Reset()
        stacksum_DOWN = data.Clone("tmp")
        stacksum_DOWN.Reset()
    
    # Legend
    x0 = 0.57; y0 = 0.6; step = 0.045
    if entry.name == 'upsilon': 
        x0 = .62; y0 = 0.65
    if region == 'WCR_OS' and entry.name == 'met':
        x0 = .2; y0 = 0.6
    nl = 3 + len(stream.MC) - 1*entry.blind + noHands*1 
    l = ROOT.TLegend(x0,y0,x0+.25,y0+nl*step)

    if entry.leg:
        leg = entry.leg
        l = ROOT.TLegend(leg[0],leg[1],leg[2],leg[3])
    l.SetTextSize(0.032)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    if not entry.blind:
        l.AddEntry(data,'Data 2012 (%.0f)' % data.Integral(0,data.GetNbinsX()+1))

    # Stack
    st = ROOT.THStack(cname,"%s; %s; " % (cname,entry.xtitle))

    # Same Sign
    SS = hist.get_SS(rf,entry,stream.data,region,kSS = rQCD)
    SS.SetLineWidth(0)
    SS.SetLineColor(3)
    SS.SetFillColor(3)
    SS.SetMarkerStyle(0)
    SS.Rebin(rb)

    if syst:
        SS_UP = hist.get_SS(rf,entry,stream.data,region,kSS = rQCD_UP)
        SS_DOWN = hist.get_SS(rf,entry,stream.data,region,kSS = rQCD_DOWN)
        SS_UP.Rebin(rb); SS_DOWN.Rebin(rb)



    if entry.xmax:
        SS.SetAxisRange(0,entry.xmax)
    
    st.Add(SS)
    stacksum.Add(SS)
    if syst:
        stacksum_DOWN.Add(SS_DOWN)
        stacksum_UP.Add(SS_UP)

    if verbose:
        print 'Data: %.0f +/- %.0f' % (hist.IntegralError(data)[0],hist.IntegralError(data)[1])
        print 'NSS: %.0f +/- %.0f' % (hist.IntegralError(SS)[0],hist.IntegralError(SS)[1])
    if noHands:
        groups = stream.MC + [stream.signal]
    else:
        groups = stream.MC + [stream.LH,stream.RH]
        cname += '_LRstack'
        pname += '_LRstack'
    mc_leg = []
    for group in groups:
        kOS = ErrorFloat(1.0,0.0); kSS = ErrorFloat(1.0,0.0)
        kOS_UP = ErrorFloat(1.,0.); kOS_DOWN = ErrorFloat(1.,0.)
        kSS_UP = ErrorFloat(1.,0.); kSS_DOWN = ErrorFloat(1.,0.) 
        if group.name == 'Wlnu': 
            kOS = kW_OS; kSS = kW_SS * rQCD
            if syst:
                kOS_UP = kW_OS_UP; kOS_DOWN = kW_OS_DOWN
                kSS_UP = kW_SS_UP; kSS_DOWN = kW_SS_DOWN
        elif group.name == 'Ztautau_LH': kOS = kZL; kSS=kOS * rQCD
        elif group.name == 'Ztautau_RH': kOS = kZR; kSS=kOS * rQCD
        else: 
            if syst: 
                kSS_UP = kSS_UP * rQCD_UP
                kSS_DOWN = kSS_DOWN * rQCD_DOWN
            else:
                kSS = kSS * rQCD 
        h = hist.get_OS_SS(rf,entry,group,region,kOS=kOS,kSS=kSS)
        if syst:
            h_up = hist.get_OS_SS(sysf,entry,group, region+'_%s_UP' %syst, kOS = kOS_UP, kSS = kSS_UP)
            h_down = hist.get_OS_SS(sysf,entry,group, region+'_%s_DOWN' %syst, kOS = kOS_DOWN, kSS = kSS_DOWN)
            h_up.Rebin(rb)
            h_down.Rebin(rb)

        h.Rebin(rb)
        if entry.xmax:
            h.SetAxisRange(0,entry.xmax)
        h.SetLineColor(entry[group.name]['fill'])
        h.SetLineWidth(0)
        h.SetMarkerStyle(0)
        h.SetFillColor(entry[group.name]['fill'])
        if 'real' in group.name: h.SetFillStyle(3001)
        if 'fake' in group.name: h.SetFillStyle(3004)
        h.SetLabelSize(0.04,"y")
    
        st.Add(h)
        stacksum.Add(h)
        if syst:
            stacksum_UP.Add(h_up)
            stacksum_DOWN.Add(h_down)
        mc_leg.append( (h,group.legend,hist.IntegralError(h)[0],hist.IntegralError(h)[1]))
        #l.AddEntry(h,'%s (OS-SS): %.0f +/- %.0f' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1]))
        if verbose:
            print '%s: %.0f +/- %.0f' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1])

    for li in reversed(mc_leg):
        l.AddEntry(li[0], '%s (OS-SS): %.0f +/- %.0f' % (li[1],li[2],li[3]))

    l.AddEntry(SS,'SS: %.0f +/- %.0f' % (hist.IntegralError(SS)[0],hist.IntegralError(SS)[1]))

    data.SetXTitle(entry.xtitle)
    data.GetXaxis().SetTitleOffset(1.25)
    data.SetTitle('')

    stacksum.SetFillStyle(3354)
    stacksum.SetLineColor(1)
    stacksum.SetFillColor(12)
    stacksum.SetMarkerStyle(1)

    if syst:
        sys_error = hist.get_sys_error(stacksum, stacksum_DOWN, stacksum_UP)
        sys_error.SetFillColor(12)
        sys_error.SetLineWidth(2)
        sys_error.SetFillStyle(3019)


    maxi = max(data.GetMaximum(),st.GetMaximum())
    mini = 0.
    data.SetMinimum(1.2*mini)
    data.SetMaximum(1.5*maxi)

    if syst:
        if not entry.blind:
            data.Draw()
            st.Draw("HIST,SAME")
            sys_error.Draw("E2,SAME")
            data.Draw("SAME, PE")
        else:
            data.Reset()
            data.Draw()
            st.SetTitle('')
            
            st.SetMinimum(1.2*st.GetMinimum())
            st.SetMaximum(1.5*st.GetMaximum())
            st.Draw("HIST,SAME")
            sys_error.Draw("E2,SAME")

    else:
        if not entry.blind:
            data.Draw()
            st.Draw("HIST,SAME")
            stacksum.Draw("E2,SAME")
            data.Draw("SAME, PE")
        else:
            data.Reset()
            data.Draw()
            st.SetTitle('')
            
            st.SetMinimum(1.2*st.GetMinimum())
            st.SetMaximum(1.5*st.GetMaximum())
            st.Draw("HIST,SAME")
            stacksum.Draw("E2,SAME")


    c.RedrawAxis()
    l.Draw()

    if ratio:
        c.cd(2)
        ROOT.gPad.SetTopMargin(0)
        dataratio = data.Clone('ratio')
        if entry.xmax:
            dataratio.SetAxisRange(0,entry.xmax)
        dataratio.GetYaxis().SetTitle('Data / MC')
        dataratio.GetYaxis().SetTitleOffset(0.25)
        dataratio.GetYaxis().SetTitleSize(0.16)
        dataratio.SetLabelSize(.1,"y")
        dataratio.Divide(stacksum)
        dataratio.SetMaximum(1.5)
        dataratio.SetMinimum(0.5)
        
        err = stacksum.Clone('err')
        err.Divide(err)
        err.SetFillStyle(3354)
        err.SetLineColor(1)
        err.SetFillColor(12)
        err.SetMarkerStyle(1)
        xmin = dataratio.GetXaxis().GetXmin()
        xmax = entry.xmax or dataratio.GetXaxis().GetXmax()
        cl = ROOT.TLine(xmin,1,xmax,1)
        cl.SetLineStyle(3)
        dataratio.GetXaxis().SetLabelSize(0.04)
        dataratio.Draw('PE')
        cl.Draw("SAME")
        err.Draw("SAME,E2")
        c.Update()

    c.Update()
    c.Print('%s/%s.gif' % (PLOTDIR,cname))
    c.Print('%s/%s.gif+1' % (PLOTDIR,pname))












