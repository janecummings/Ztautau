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


def charged_asymmetry(mu):
    rf = ROOT.TFile('~/panalysis/output/prodr/%s_ttbar/%s_regions.root' % (mu.name,mu.name))    
    entry = upsilon_20
    k_entry = upsilon_CR
    verbose = False
    region = 'TTBAR_OS'
    # scale factors
    k_entry = upsilon_CR

    groups = mu.MC + [mu.signal]
    kW_OS = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, groups = groups)
    kW_SS = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose, groups = groups)
    rQCD = hist.calc_rQCD(rf,mu,k_entry,'QCD',kW_OS = kW_OS, kW_SS = kW_SS, verbose = verbose, groups = groups)

    # Canvas
    if entry.blind: ratio = False
    ratio = False
    cname = '%s_%s_%s' % (mu.name,entry.name,region)
    #pname = '%s_%s_%s' % (mu.name,region.rstrip('_OS'), hand)
    if ratio:
        c = ROOT.TCanvas(cname,cname,1500,1000)
        #Split for ratio plot
        c.Divide(1,2)
        c.cd(1).SetPad(0.0,0.2,1.0,1.0)
        c.cd(2).SetPad(0.0,0.0,1.0,0.2)
        c.cd(1)
    else:
        c = ROOT.TCanvas(cname,cname, 1200,800)
    rb = 1
    # data  
    if 'SR_OS' in region: dregion = 'SR_OS'
    if 'SR_SS' in region: dregion = 'SR_SS'
    data = rf.Get('h_%s_%s_%s' % (entry.name,mu.data.name,dregion)).Clone()
    data.SetMarkerStyle(20)
    data.SetMarkerSize(2)
    data.SetLabelSize(0.04,"y")
    data.Rebin(rb)
    # Blind data
    data.Reset()
    if entry.xmax:
        data.SetAxisRange(entry.xmin or 0,entry.xmax)
    # stat. unc. on stack
    stacksum = data.Clone("tmp")
    stacksum.Reset()

    # Legend
    x0 = 0.73; y0 = 0.55; step = 0.045
    
    nl = 4 + len(mu.MC) - 1*entry.blind 
    l = ROOT.TLegend(x0,y0,x0+.16,y0+nl*step)
    l.SetTextSize(0.0345)
    l.SetBorderSize(0)
    l.SetFillColor(0)

    # Data
    # Stack
    st = ROOT.THStack(cname,"%s; %s; " % (cname,entry.xtitle))

    # QCD ESTIMATE
    #if region == 'SR_SS': rQCD = ErrorFloat(1.,0.)
    qcd = hist.get_QCDW(rf,mu,entry,'SR_SS', rqcd = rQCD , verbose = False, groups = groups)
    qcd.SetLineWidth(0)
    qcd.SetLineColor(3)
    qcd.SetFillColor(3)
    qcd.SetMarkerStyle(0)
    qcd.SetMarkerColor(3)
    qcd.Rebin(rb)


    if entry.xmax:
        qcd.SetAxisRange(entry.xmin or 0,entry.xmax)
    
    st.Add(qcd)
    stacksum.Add(qcd)


    if verbose:
        total_est = ErrorFloat(*hist.IntegralError(qcd))
        print 'Data: %.0f +/- %.0f' % (hist.IntegralError(data)[0],hist.IntegralError(data)[1])
        print 'NQCD: %.0f +/- %.0f' % (hist.IntegralError(qcd)[0],hist.IntegralError(qcd)[1])

    mc_leg = []

    n_total = ErrorFloat(*hist.IntegralError(hist.get_OS(rf,entry,mu.signal,region)))

    for group in groups:
        entry.add_group(group)
        if group.name == 'Wlnu': 
            h = hist.get_dataW(rf, mu, entry, 'WCR_OS', verbose = False, groups = groups)
        else:
            h = hist.get_OS(rf,entry,group,region)
        h.Rebin(rb)

        if entry.xmax:
            h.SetAxisRange(entry.xmin or 0,entry.xmax)

        h.SetLineColor(entry[group.name]['fill'])
        h.SetLineWidth(0)
        h.SetMarkerStyle(0)
        h.SetMarkerColor(entry[group.name]['fill'])
        h.SetFillColor(entry[group.name]['fill'])
        h.SetLabelSize(0.04,"y")

        st.Add(h)
        stacksum.Add(h)
        mc_leg.append( (h,group.legend,hist.IntegralError(h)[0],hist.IntegralError(h)[1]))

        if verbose:
            print '%s: %.0f +/- %.0f' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1])
            total_est += ErrorFloat(*hist.IntegralError(h))

    if verbose:
        print 'Total Estimated: %.0f +/- %.0f' % (total_est.val,total_est.err)

    for li in reversed(mc_leg):
        l.AddEntry(li[0], '%s'% (li[1]))

    l.AddEntry(qcd,'Multijet')#: %.0f +/- %.0f' % (hist.IntegralError(qcd)[0],hist.IntegralError(qcd)[1]))
    data.SetXTitle(entry.xtitle)
    data.GetXaxis().SetTitleOffset(1.25)

    stacksum.SetFillStyle(3354)
    stacksum.SetLineColor(1)
    stacksum.SetFillColor(12)
    stacksum.SetMarkerStyle(1)

    maxi = max(data.GetMaximum(),st.GetMaximum())
    maxi = 6250
    mini = 0.
    data.SetMinimum(1.2*mini)
    data.SetMaximum(maxi)

    bw = data.GetXaxis().GetBinWidth(1)
    data.GetYaxis().SetTitle('Events / %s %s' % (bw, entry.units))

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
        stacksum.Draw("E2,SAME")

    ROOT.gPad.RedrawAxis()
    l.Draw()

    atlas = ROOT.ATLASLabel(0.32,0.86,"   Work in Progress")
    td = ROOT.TPaveText(0.54,0.68,0.67,0.85,"NDC")
    td.SetTextAlign(13)
    td.SetBorderSize(0)
    td.SetFillColor(0)
    td.SetTextSize(0.04)

    td.AddText('#sqrt{s} = 8 TeV')
    td.AddText('#int Ldt = 20.3 fb^{-1}')
    td.Draw()
    # channel
    tm = ROOT.TPaveText(0.19,0.83,0.22,0.90,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.05)
    if mu.name == 'mu':
        tm.AddText('#mu#tau_{had}')
    elif mu.name == 'el':
        tm.AddText('e#tau_{had}')
    tm.Draw()

    c.Update()
    c.Print('../plots/ttbar/%s.eps' % (cname))


def kinematic_plots(mu,entry,ratio=True, verbose = False):
    rf = ROOT.TFile('~/panalysis/output/prod4/%s_fit/%s_regions.root' % (mu.name,mu.name))
    region = 'SR_OS'

    # scale factors
    k_entry = upsilon_CR
    kW_OS = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose)
    kW_SS = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose)
    rQCD = hist.calc_rQCD(rf,mu,k_entry,'QCD',kW_OS = kW_OS, kW_SS = kW_SS, verbose = verbose)
    # left and right handed
    print rQCD
    #sys.exit()
    rb = entry.rebin

    # Canvas
    if entry.blind: ratio = False
    #ratio = False
    cname = '%s_%s_%s' % (mu.name,entry.name,region)
    pname = '%s_%s' % (mu.name,region.rstrip('_OS'))
    if ratio:
        c = ROOT.TCanvas(cname,cname,1500,1000)
        #Split for ratio plot
        c.Divide(1,2)
        c.cd(1).SetPad(0.0,0.2,1.0,1.0)
        c.cd(2).SetPad(0.0,0.0,1.0,0.2)
        c.cd(1)
    else:
        c = ROOT.TCanvas(cname,cname, 1200,800)
    # data  
    dregion = region
    data = rf.Get('h_%s_%s_%s' % (entry.name,mu.data.name,dregion)).Clone()
    data.SetMarkerStyle(20)
    data.SetMarkerSize(2)
    data.SetLabelSize(0.04,"y")
    data.Rebin(rb)
    if entry.xmax:
        data.SetAxisRange(entry.xmin or 0,entry.xmax)
    # stat. unc. on stack
    stacksum = data.Clone("tmp")
    stacksum.Reset()

    # Legend
    x0 = 0.78; y0 = 0.5; step = 0.05
  
    nl = 4 + len(mu.MC) 
    l = ROOT.TLegend(x0,y0,x0+.16,y0+nl*step)
    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)

    if not entry.blind:
        l.AddEntry(data,'Data 2012', "P")# (%.0f)' % data.Integral(0,data.GetNbinsX()+1))
    # Stack
    st = ROOT.THStack(cname,"%s; %s; " % (cname,entry.xtitle))

    # QCD ESTIMATE
    #if region == 'SR_SS': rQCD = ErrorFloat(1.,0.)
    #qcd = hist.get_QCD(rf,mu,entry,'TTBAR_SS', kW = kW_SS, rqcd = rQCD)
    qcd = hist.get_QCD(rf,mu,entry,'SR_SS', kW = kW_SS, rqcd = rQCD)
    qcd.SetLineWidth(0)
    qcd.SetLineColor(3)
    qcd.SetFillColor(3)
    qcd.SetMarkerStyle(0)
    qcd.SetMarkerColor(3)
    qcd.Rebin(rb)

    if entry.xmax:
        qcd.SetAxisRange(entry.xmin or 0,entry.xmax)
    
    st.Add(qcd)
    stacksum.Add(qcd)

    if verbose:
        total_est = ErrorFloat(*hist.IntegralError(qcd))
        print 'Data: %.0f +/- %.0f' % (hist.IntegralError(data)[0],hist.IntegralError(data)[1])
        print 'NQCD: %.0f +/- %.0f' % (hist.IntegralError(qcd)[0],hist.IntegralError(qcd)[1])

    mc_leg = []

    n_total = ErrorFloat(*hist.IntegralError(hist.get_OS(rf,entry,mu.signal,region)))

    groups = mu.MC + [mu.signal]

    for group in groups:
        entry.add_group(group)
        kOS = ErrorFloat(1.0,0.0);kSS = ErrorFloat(1.0,0.0);
        if group.name == 'Wlnu': 
            kOS = kW_OS;
            kSS = kW_SS
        if 'SS' not in region:
            h = hist.get_OS(rf,entry,group,region,kOS=kOS)
        else:
            h = hist.get_SS(rf,entry,group,region,kSS=kSS)
        h.Rebin(rb)

        if entry.xmax:
            h.SetAxisRange(entry.xmin or 0,entry.xmax)


        h.SetLineColor(entry[group.name]['fill'])
        h.SetLineWidth(0)
        h.SetMarkerStyle(0)
        h.SetMarkerColor(entry[group.name]['fill'])
        h.SetFillColor(entry[group.name]['fill'])
        h.SetLabelSize(0.04,"y")

        st.Add(h)
        stacksum.Add(h)
        mc_leg.append( (h,group.legend,hist.IntegralError(h)[0],hist.IntegralError(h)[1]))

        if verbose:
            print '%s: %.0f +/- %.0f' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1])
            total_est += ErrorFloat(*hist.IntegralError(h))
            if syst:
                print '%s: %.0f +/- %.0f' % ('up',hist.IntegralError(h_up)[0],hist.IntegralError(h_up)[1])
                print '%s: %.0f +/- %.0f' % ('down',hist.IntegralError(h_down)[0],hist.IntegralError(h_down)[1])
    if verbose:
        print 'Total Estimated: %.0f +/- %.0f' % (total_est.val,total_est.err)


    for li in reversed(mc_leg):
        l.AddEntry(li[0], '%s'% (li[1]))


    #l.AddEntry(qcd,'QCD: %.0f +/- %.0f' % (hist.IntegralError(qcd)[0],hist.IntegralError(qcd)[1]))
    if 'OS' in region:
        l.AddEntry(qcd,'Multijet')#: %.0f +/- %.0f' % (hist.IntegralError(qcd)[0],hist.IntegralError(qcd)[1]))


    data.SetXTitle(entry.xtitle)
    data.GetXaxis().SetTitleOffset(1.25)

    stacksum.SetFillStyle(3354)
    stacksum.SetLineColor(1)
    stacksum.SetFillColor(12)
    stacksum.SetMarkerStyle(1)


    maxi = max(data.GetMaximum(),st.GetMaximum())
    mini = 0.
    data.SetMinimum(1.2*mini)
    data.SetMaximum(1.4*maxi)

    bw = data.GetXaxis().GetBinWidth(1)

    data.GetYaxis().SetTitle('Events / %s %s' % (bw, entry.units))

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
        stacksum.Draw("E2,SAME")


    data.Draw("SAME,PE")

    ROOT.gPad.RedrawAxis()
    l.Draw()

    atlas = ROOT.ATLASLabel(0.3,0.86,"   Work in Progress")
    #td = ROOT.TPaveText(0.3,0.75,0.7,0.85,"NDC")
    if 'OS' in region:
        td = ROOT.TPaveText(0.52,0.68,0.65,0.85,"NDC")
    else:
        td = ROOT.TPaveText(0.19,0.68,0.32,0.85,"NDC")
    td.SetTextAlign(13)
    td.SetBorderSize(0)
    td.SetFillColor(0)
    td.SetTextSize(0.04)

    #td.AddText('#sqrt{s} = 8 TeV    #int Ldt = 20.3 fb^{-1}')
    #td.AddText('#sqrt{s} = 8 TeV')
    #td.AddText('#int Ldt = 20.3 fb^{-1}')
    td.Draw()
    # channel
    tm = ROOT.TPaveText(0.19,0.83,0.22,0.90,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.05)
    if mu.name == 'mu':
        tm.AddText('#mu#tau_{had}')
    elif mu.name == 'el':
        tm.AddText('e#tau_{had}')
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
        dataratio.SetMaximum(1.6)
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
    c.Print('../plots/ttbar/%s.eps' % (cname))
    pname = '%s_%s' % (mu.name, region)

    c.SaveAs('../plots/ttbar/%s.gif+1'% pname)


def kinematic_plots_SR(mu,entry,ratio=True, verbose = False):
    rf = ROOT.TFile('~/panalysis/output/prod4/%s_fit/%s_regions.root' % (mu.name,mu.name))
    region = 'SR_OS'

    # scale factors
    k_entry = upsilon_CR
    kW_OS = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose)
    kW_SS = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose)
    rQCD = hist.calc_rQCD(rf,mu,k_entry,'QCD',kW_OS = kW_OS, kW_SS = kW_SS, verbose = verbose)
    # left and right handed
    kW_OS_L = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, groups = mu.MC + [mu.LH] )
    kW_OS_R = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, groups = mu.MC + [mu.RH] )
    kW_SS_L = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose, groups = mu.MC + [mu.LH] )
    kW_SS_R = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose, groups = mu.MC + [mu.RH] )    
    rQCD_L = hist.calc_rQCD(rf,mu,k_entry,'QCD',kW_OS = kW_OS_L, kW_SS = kW_SS_L, verbose = verbose, groups = mu.MC + [mu.LH] )
    rQCD_R = hist.calc_rQCD(rf,mu,k_entry,'QCD',kW_OS = kW_OS_R, kW_SS = kW_SS_R, verbose = verbose, groups =mu.MC + [mu.RH] )

    rb = entry.rebin

    # Canvas
    if entry.blind: ratio = False
    ratio = False
    cname = '%s_%s_%s' % (mu.name,entry.name,region)
    pname = '%s_%s' % (mu.name,region.rstrip('_OS'))
    if ratio:
        c = ROOT.TCanvas(cname,cname,1500,1000)
        #Split for ratio plot
        c.Divide(1,2)
        c.cd(1).SetPad(0.0,0.2,1.0,1.0)
        c.cd(2).SetPad(0.0,0.0,1.0,0.2)
        c.cd(1)
    else:
        c = ROOT.TCanvas(cname,cname, 1200,800)

    # Data
    dregion = region
    data = rf.Get('h_%s_%s_%s' % (entry.name,mu.data.name,dregion)).Clone()
    data.SetMarkerStyle(20)
    data.SetMarkerSize(2)
    data.SetLabelSize(0.04,"y")
    data.Rebin(rb)
    if entry.xmax:
        
        data.SetAxisRange(entry.xmin or 0,entry.xmax)

    # stat. unc. on stack
    stacksum = data.Clone("tmp")
    stacksum.Reset()

    # Legend
    x0 = 0.72; y0 = 0.5; step = 0.05
  
    nl = 4 + len(mu.MC) - 1*entry.blind 
    l = ROOT.TLegend(x0,y0,x0+.17,y0+nl*step)
    if entry.leg:
        leg = entry.leg
        l = ROOT.TLegend(leg[0],leg[1],leg[2],leg[3])

    l.SetTextSize(0.035)
    l.SetBorderSize(0)
    l.SetFillColor(0)

    if not entry.blind:
        l.AddEntry(data,'Data 2012', "P")# (%.0f)' % data.Integral(0,data.GetNbinsX()+1))

    # Stack
    st = ROOT.THStack(cname,"%s; %s; " % (cname,entry.xtitle))

    # QCD ESTIMATE
    #if region == 'SR_SS': rQCD = ErrorFloat(1.,0.)
    qcd = hist.get_QCD(rf,mu,entry,'SR_SS', kW = kW_SS, rqcd = rQCD)
    qcd.SetLineWidth(0)
    qcd.SetLineColor(3)
    qcd.SetFillColor(3)
    qcd.SetMarkerStyle(0)
    qcd.SetMarkerColor(3)
    qcd.Rebin(rb)

    # left and right
    lefthand = hist.get_QCD(rf,mu,entry,'SR_SS', kW = kW_SS_L, rqcd = rQCD_L, groups = mu.MC + [mu.LH])
    righthand = hist.get_QCD(rf,mu,entry,'SR_SS', kW = kW_SS_R, rqcd = rQCD_R, groups = mu.MC + [mu.RH])
    lefthand.Rebin(rb)
    righthand.Rebin(rb)
    if not 'SR' in region:
        qcd.Reset()
        lefthand.Reset()
        righthand.Reset()

    if entry.xmax:
        qcd.SetAxisRange(entry.xmin or 0,entry.xmax)
        lefthand.SetAxisRange(entry.xmin or 0 ,entry.xmax)
        righthand.SetAxisRange(entry.xmin or 0, entry.xmax)
    
    st.Add(qcd)
    stacksum.Add(qcd)

    if verbose:
        total_est = ErrorFloat(*hist.IntegralError(qcd))
        print 'Data: %.0f +/- %.0f' % (hist.IntegralError(data)[0],hist.IntegralError(data)[1])
        print 'NQCD: %.0f +/- %.0f' % (hist.IntegralError(qcd)[0],hist.IntegralError(qcd)[1])

    mc_leg = []
    n_total = ErrorFloat(*hist.IntegralError(hist.get_OS(rf,entry,mu.signal,region)))
    groups = mu.MC + [mu.signal]

    for group in groups:
        entry.add_group(group)
        kOS = ErrorFloat(1.0,0.0);kSS = ErrorFloat(1.0,0.0);
        if group.name == 'Wlnu': 
            kOS = kW_OS;
            kSS = kW_SS
        if 'SS' not in region:
            h = hist.get_OS(rf,entry,group,region,kOS=kOS)
        else:
            h = hist.get_SS(rf,entry,group,region,kSS=kSS)
        h.Rebin(rb)

        if entry.xmax:
            h.SetAxisRange(entry.xmin or 0,entry.xmax)


        if group.name != mu.signal.name:    
            lefthand.Add(h)
            righthand.Add(h)

        h.SetLineColor(entry[group.name]['fill'])
        h.SetLineWidth(0)
        h.SetMarkerStyle(0)
        h.SetMarkerColor(entry[group.name]['fill'])
        h.SetFillColor(entry[group.name]['fill'])
        h.SetLabelSize(0.04,"y")

        st.Add(h)
        stacksum.Add(h)
        mc_leg.append( (h,group.legend,hist.IntegralError(h)[0],hist.IntegralError(h)[1]))

        if verbose:
            print '%s: %.0f +/- %.0f' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1])
            total_est += ErrorFloat(*hist.IntegralError(h))
            if syst:
                print '%s: %.0f +/- %.0f' % ('up',hist.IntegralError(h_up)[0],hist.IntegralError(h_up)[1])
                print '%s: %.0f +/- %.0f' % ('down',hist.IntegralError(h_down)[0],hist.IntegralError(h_down)[1])
    if verbose:
        print 'Total Estimated: %.0f +/- %.0f' % (total_est.val,total_est.err)


    ## Add left and right
    hLeft = hist.get_OS(rf,entry,mu.LH,region)
    hRight = hist.get_OS(rf,entry,mu.RH,region)
    for h in [hLeft,hRight]:
        h.Rebin(rb)
        h.SetMarkerStyle(0)
        h.SetLineWidth(1)
        if entry.xmax:
            h.SetAxisRange(entry.xmin or 0,entry.xmax)
    lefthand.Add(hLeft)
    righthand.Add(hRight)

    lefthand.SetLineColor(mu.LH.fill)
    righthand.SetLineColor(mu.RH.fill)


    lefthand.SetLineStyle(2)
    righthand.SetLineStyle(2)

    for li in reversed(mc_leg):
        l.AddEntry(li[0], '%s'% (li[1]))
    #l.AddEntry(qcd,'QCD: %.0f +/- %.0f' % (hist.IntegralError(qcd)[0],hist.IntegralError(qcd)[1]))
    if 'SR' in region:
        l.AddEntry(qcd,'Multijet')#: %.0f +/- %.0f' % (hist.IntegralError(qcd)[0],hist.IntegralError(qcd)[1]))
    l.AddEntry(lefthand, 'Left-Handed', "L")
    l.AddEntry(righthand, 'Right-Handed', "L")

    data.SetXTitle(entry.xtitle)
    data.GetXaxis().SetTitleOffset(1.25)

    stacksum.SetFillStyle(3354)
    stacksum.SetLineColor(1)
    stacksum.SetFillColor(12)
    stacksum.SetMarkerStyle(1)


    maxi = max(data.GetMaximum(),st.GetMaximum())
    mini = 0.
    data.SetMinimum(1.2*mini)
    data.SetMaximum(1.4*maxi)

    bw = data.GetXaxis().GetBinWidth(1)

    data.GetYaxis().SetTitle('Events / %s %s' % (bw, entry.units))



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
        stacksum.Draw("E2,SAME")

    lefthand.Draw("HIST,SAME")
    righthand.Draw("HIST,SAME")
    data.Draw("SAME,PE")

    ROOT.gPad.RedrawAxis()
    l.Draw()

    atlas = ROOT.ATLASLabel(0.3,0.87,"   Work in Progress")
    #td = ROOT.TPaveText(0.3,0.75,0.7,0.85,"NDC")
    if 'SR' in region:
        td = ROOT.TPaveText(0.52,0.68,0.65,0.85,"NDC")
    else:
        td = ROOT.TPaveText(0.19,0.68,0.32,0.85,"NDC")
    td.SetTextAlign(13)
    td.SetBorderSize(0)
    td.SetFillColor(0)
    td.SetTextSize(0.04)

    #td.AddText('#sqrt{s} = 8 TeV    #int Ldt = 20.3 fb^{-1}')
    #td.AddText('#sqrt{s} = 8 TeV')
    #td.AddText('#int Ldt = 20.3 fb^{-1}')
    td.Draw()
    # channel
    tm = ROOT.TPaveText(0.19,0.85,0.22,0.92,"NDC")
    tm.SetBorderSize(0)
    tm.SetFillColor(0)
    tm.SetTextSize(0.05)
    if mu.name == 'mu':
        tm.AddText('#mu#tau_{had}')
    elif mu.name == 'el':
        tm.AddText('e#tau_{had}')
    tm.Draw()

    #atlas.Draw()

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
    c.Print('../plots/ttbar/%s.eps' % (cname))
    #pname = '%s_%s' % (mu.name, region)

if __name__ == '__main__':
    tau_leadTrkPt.blind = True
    kinematic_plots_SR(streams.mu, tau_leadTrkPt) 
