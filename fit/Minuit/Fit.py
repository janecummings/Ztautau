import sys,os
import ROOT
import math
import run
from array import array
import Drawer,hist
from ErrorFloat import ErrorFloat
from template import template,component,fake_data,entry
# Atlas Style
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.LoadMacro("../draw/AtlasStyle.C")
ROOT.gROOT.LoadMacro("../draw/AtlasLabels.C")
ROOT.SetAtlasStyle()
ROOT.gStyle.SetFrameBorderSize(1)


from printing import row
fcol=20; col=20


comment = False

##________________________________________________________________________________________________________________
print '\n------ OUTPUT FILES ------\n'
# get root files
s = run.testFit
stream = s.stream
rf = ROOT.TFile('%s/%s_regions.root' % (s.OUTPUTDIR,s.stream.name))
kCR = ROOT.TFile('%s/%s_CR.root' % (s.OUTPUTDIR,s.stream.name))
#entry = template.entry
for group in stream.groups: entry.add_group(group)

print '*' * 120
print '** Histogram File    =  %s/%s_regions.root' % (s.OUTPUTDIR,s.stream.name)
print '** Control File      =  %s/%s_CR.root' % (s.OUTPUTDIR,s.stream.name)
print '** stream            =  %s' % stream.name
print '** entry name        =  %s' % entry.name
print '*' * 120

# Set number of bins
rb = 10
nbins = 30/rb

# temporary until we have the separate LH and RH samples
p = -0.15
## acceptance normalization
kLH = ErrorFloat(1/(0.5*(1-p)),0.)
kRH = ErrorFloat(1/(0.5*(1+p)),0.)
##________________________________________________________________________________________________________________
print '\n------ FIT METHOD ------\n'
#stat = 1; nPar0 = 2; nPar = nPar0 # simple two parameter fit
#stat = 2; nPar0 = 3; nPar = nPar0 # 2 paramter fit with N_SS; kW_OS and kW_SS with nominal rQCD
#stat = 3; nPar0 =  3; nPar = nPar0 + nbins * (6 * 14 + 2 + 1) # statistical uncertainty on template
#stat = 4; nPar0 =  3; nPar = nPar0 + nbins * (6 * 14 + 2 + 1)  # calculate kW_OS and kW_SS in fit
stat = 5; nPar0 =  2; nPar = nPar0 + nbins * (6 * 14 + 2 + 1)  # calculate kW_OS, kW_SS, nSS in fit

print '*' * 120
print '** stat    =  %i' % stat
print '** nbins   =  %i' % nbins
print '** nPar0   =  %i' % nPar0
print '** nPar    =  %i' % nPar
print '*' * 120
##________________________________________________________________________________________________________________
print '\n------ INITIALIZE ------\n'
# initialize template and components 
print '... Initialize template'
template = template(rf,kCR,stream,fitstat=stat,rb=rb)
# load components / fill and rebin histograms 
print '... Loading components'
template.fill()
# get cross sections and totalEvents for scale factors
print '... Fill cross sections and total events'
template.load()
# scale the left and right signal to the correct cross-section (temp)
print '... Scale LH and RH template components to cross section'

### let's just move these scale factors into the fit equation where we can see them
template.LH.scale(kLH.val)
template.RH.scale(kRH.val)
print '... Populate parameter indices'
template.index(nPar0,nbins)
##________________________________________________________________________________________________________________
print '\n------ NOMINAL SCALE FACTORS ------\n'
# scale factors 
kW_OS = hist.calc_kX(kCR,stream,entry,'WCR','Wlnu',OS='OS',verbose=True)
kW_SS = hist.calc_kX(kCR,stream,entry,'WCR','Wlnu',OS='SS',verbose=True)
rQCD = hist.calc_rQCD(kCR,stream,entry,'QCD',kW_OS = kW_OS, kW_SS = kW_SS,verbose = True)
print '*' * 120
print '** kW_OS    =  %s' % kW_OS
print '** kW_SS    =  %s' % kW_SS
print '** rQCD     =  %s' % rQCD 
print '*' * 120

#if stat == 4:
#    global kW_OS_fit
#    global kW_SS_fit

# make fake data for fit
data = fake_data(rf,stream,kW_OS,kW_SS,rQCD,rb)

##________________________________________________________________________________________________________________
logFact = [0.]
for x in xrange(1,50000):
	logFact.append(logFact[x-1]+math.log(x))

##________________________________________________________________________________________________________________
def fcn(npar,gin,f,par,ierflg):
    # npar is the number of parameters; 
    # f[0] is the minimized function; 
    # par is the list of parameters; 
    # ierflg is how errors are calculated

    f[0] = 0.
    if stat in [1,2]:
        for bin in xrange(nbins): f[0] += logPoisson(bin,par)
    if stat in [3,4,5]:
        #print row('bin',['f from logPoisson','f from logStats'],fcol=fcol,col=col)
        for bin in xrange(nbins):
            #print '*' * 160
            #prow = []
            f[0] += logPoisson(bin,par) 
            #prow.append(f[0])
            f[0] += logStatsTempl(bin,par)
            #prow.append(f[0])
            #print row(bin,prow,fcol=fcol,col=col)

def logPoisson(bin,par):

    #par [0] = norm
    #par [1] = pol
    #par [2] = NSS
    
    if stat == 1: 
        lh_ = 0.5 * (1 - par[1]) * ( template.LH.GTH_OS.GetBinContent(bin + 1) - rQCD.val  * template.LH.GTH_SS.GetBinContent(bin + 1) )
        rh_ = 0.5 * (1 + par[1]) * ( template.RH.GTH_OS.GetBinContent(bin + 1) - rQCD.val * template.RH.GTH_SS.GetBinContent(bin + 1) )
        w_  = kW_OS.val * template.Wlnu.GTH_OS.GetBinContent(bin + 1) - kW_SS.val * rQCD.val * template.Wlnu.GTH_SS.GetBinContent(bin + 1)
        z_ = template.Zll.GTH_OS.GetBinContent(bin + 1) - rQCD.val * template.Zll.GTH_SS.GetBinContent(bin + 1)
        tt_ = template.ttbar.GTH_OS.GetBinContent(bin + 1) - rQCD.val * template.ttbar.GTH_SS.GetBinContent(bin + 1)
        ss_ = rQCD.val * template.SS.GTH_SS.GetBinContent(bin + 1)
 
        sb = par[0] * ( lh_ + rh_ + w_ + z_ + tt_ + ss_)

    if stat == 2:
        lh_ = 0.5 * (1 - par[1]) * ( template.LH.GTH_OS.GetBinContent(bin + 1) - par[2] * template.LH.GTH_SS.GetBinContent(bin + 1) )
        rh_ = 0.5 * (1 + par[1]) * ( template.RH.GTH_OS.GetBinContent(bin + 1) - par[2] * template.RH.GTH_SS.GetBinContent(bin + 1) )
        w_  = kW_OS * template.Wlnu.GTH_OS.GetBinContent(bin + 1) - kW_SS * par[2] * template.Wlnu.GTH_SS.GetBinContent(bin + 1)
        z_ = template.Zll.GTH_OS.GetBinContent(bin + 1) - par[2] * template.Zll.GTH_SS.GetBinContent(bin + 1)
        tt_ = template.ttbar.GTH_OS.GetBinContent(bin + 1) - par[2] * template.ttbar.GTH_SS.GetBinContent(bin + 1)
        ss_ = par[2] * template.SS.GTH_SS.GetBinContent(bin + 1)
        
        sb =  par[0] * ( lh_ + rh_ + w_ + z_ + tt_) + ss_

    if stat == 3:
        lh_ = add_term( template.LH ,   bin , par , norm = 0.5 * (1-par[1]) * kLH.val, kSS = par[2] )
        rh_ = add_term( template.RH ,   bin , par , norm = 0.5 * (1+par[1]) * kRH.val , kSS = par[2] )   
        w_  = add_term( template.Wlnu , bin , par , kOS = kW_OS.val , kSS = kW_SS.val * par[2] )   
        z_  = add_term( template.Zll ,  bin , par , kSS = par[2] )   
        tt_ = add_term( template.ttbar, bin , par , kSS = par[2] ) 
        ss_ = par[2] * par[template.SS.index_SS + bin]

        sb =  par[0] * ( lh_ + rh_ + w_ + z_ + tt_) + ss_

    if stat == 4:
        lh_ = add_term( template.LH ,   bin , par , norm = 0.5 * (1-par[1]) * kLH.val , kSS = par[2] )
        rh_ = add_term( template.RH ,   bin , par , norm = 0.5 * (1+par[1]) * kRH.val , kSS = par[2] )   
        w_  = add_term( template.Wlnu , bin , par , kOS = calc_kW_OS(par) , kSS = calc_kW_SS(par) * par[2] )   
        z_  = add_term( template.Zll ,  bin , par , kSS = par[2] )   
        tt_ = add_term( template.ttbar, bin , par , kSS = par[2] ) 
        ss_ = par[2] * par[template.SS.index_SS + bin]

        sb =  par[0] * ( lh_ + rh_ + w_ + z_ + tt_) + ss_

    if stat == 5:
        print 'HELLO 5'
        rqcd = calc_rQCD(par)
        lh_ = add_term( template.LH ,   bin , par , norm = 0.5 * (1-par[1]) * kLH.val , kSS = rqcd )
        rh_ = add_term( template.RH ,   bin , par , norm = 0.5 * (1+par[1]) * kRH.val , kSS = rqcd )   
        w_  = add_term( template.Wlnu , bin , par , kOS = calc_kW_OS(par) , kSS = calc_kW_SS(par) * rqcd )   
        z_  = add_term( template.Zll ,  bin , par , kSS = rqcd )   
        tt_ = add_term( template.ttbar, bin , par , kSS = rqcd ) 
        ss_ = rqcd * par[template.SS.index_SS + bin]

        sb =  par[0] * ( lh_ + rh_ + w_ + z_ + tt_) + ss_


    if sb <= 0.0: sb = 0.00000000000000001
    try:
        lp = sb - data.GetBinContent(bin + 1) * math.log(sb) + logFact[int(data.GetBinContent(bin + 1))]
    except:
        print 'Out of range in logFact in logPoisson: %i' % data.GetBinContent(bin + 1)
        sys.exit(0)
    print lp

    return lp

def calc_kW_OS(par):
    print 'HERE'
    global kW_OS_fit
    nData = template.WCR_OS
    nW = template.Wlnu.WCR_OS
    nB = 0.
    nB += template.LH.WCR_OS * 0.5 * (1-par[1])
    nB += template.RH.WCR_OS * 0.5 * (1+par[1])
    nB += template.ttbar.WCR_OS
    nB += template.Zll.WCR_OS
    kW_OS_fit = (nData - nB) / nW
    return kW_OS_fit

def calc_kW_SS(par):
    global kW_SS_fit
    nData = template.WCR_SS
    nW = template.Wlnu.WCR_SS
    nB = 0.
    nB += template.LH.WCR_SS * 0.5 * (1-par[1])
    nB += template.RH.WCR_SS * 0.5 * (1+par[1])
    nB += template.ttbar.WCR_SS
    nB += template.Zll.WCR_SS
    kW_SS_fit = (nData - nB) / nW
    return kW_SS_fit  


def calc_rQCD(par):
    global rQCD_fit
    nOS = template.QCD_OS
    nOS -= template.LH.QCD_OS * 0.5 * (1-par[1])
    nOS -= template.RH.QCD_OS * 0.5 * (1+par[1])
    nOS -= template.ttbar.QCD_OS
    nOS -= template.Zll.QCD_OS

    nSS = template.QCD_SS
    nSS -= template.LH.QCD_SS * 0.5 * (1-par[1])
    nSS -= template.RH.QCD_SS * 0.5 * (1+par[1])
    nSS -= template.ttbar.QCD_SS
    nSS -= template.Zll.QCD_SS

    rQCD_fit = nOS/nSS
    return rQCD_fit



def add_term(compo, bin, par, norm = 1.0, kOS = 1.0, kSS = 1.0):
    t_ = 0
    for j,mu in enumerate(compo.mu):
        t_ += norm *  mu * ( kOS * par[compo.index_OS[j] + bin] - kSS * par[compo.index_SS[j] + bin] ) 
        #t_ += norm * (mu_OS * kOS * par[compo.index_OS[j] + bin ] - mu_SS * kSS * compo.TH_SS_UNSCALED[j].GetBinContent(bin+1) ) # OS ONLY
    return t_

def logStatsTempl(bin,par):
    lss = 0

    lss += addLogStats(template.LH,bin,par)
    lss += addLogStats(template.RH,bin,par)
    lss += addLogStats(template.Wlnu,bin,par)
    lss += addLogStats(template.Zll,bin,par)
    lss += addLogStats(template.ttbar,bin,par)
    lss += logStats(par[ template.SS.index_SS + bin ],template.SS.GTH_SS.GetBinContent(bin+1))

    return lss

def logStats(fit,obs):
    ls = 0
    if fit>0:
        try:
            ls = fit - obs * math.log(fit) + logFact[int(obs)]
        except:
            print 'Out of range in logFact in logStats: obs = %i' % obs
            sys.exit()
    return ls

def addLogStats(compo,bin,par):
    ls_ = 0
    for j, (h_OS,h_SS) in enumerate(zip(compo.TH_OS_UNSCALED,compo.TH_SS_UNSCALED)):
        ls_ += logStats( par[compo.index_OS[j] + bin], h_OS.GetBinContent(bin+1) ) 
        ls_ += logStats( par[compo.index_SS[j] + bin], h_SS.GetBinContent(bin+1) ) 
    return ls_


##_______________________________________________________________________________________________________________

# initialize Minuit
print '\n------ INITIALIZE MINUIT ------\n'
m = ROOT.TMinuit(nPar) 
m.SetFCN(fcn)

print '... Set Error'
# set parameters 
erflag = ROOT.Long(0)
arglist = array('d',10*[0.])
arglist[0] = 0.5
m.mnexcm("SET ERR",arglist,1,erflag)
arglist[0] = 0
m.mnexcm("SET PRI",arglist,1,erflag)

print '... Initialize parameters'
m.mnparm(0,"NORM", 1.0,0.1,0.,3.0,erflag)
m.mnparm(1,"POL",  0.0,0.05,-2.,2.,erflag)
if stat > 1 and stat < 5:
    m.mnparm(2,"NSS",1.1,0.01,0.,2.,erflag)
if stat in [3,4,5]:
    step_ = 1.
    min_ = 0.0
    if rb > 5:
        max_ = 50000
    else: max_ = 5000
    tst = 0
    for c in template.components:
        if c.name == 'SS': continue
        for j, (h_OS,h_SS) in enumerate(zip(c.TH_OS_UNSCALED, c.TH_SS_UNSCALED) ):
            for bin in xrange(nbins):
                i_OS = c.index_OS[j] + bin;  i_SS = c.index_SS[j] + bin 
                m.mnparm(i_OS,  "%i_%s" % (i_OS,c.name), h_OS.GetBinContent(bin+1), step_, min_, max_, erflag )
                m.mnparm(i_SS,  "%i_%s" % (i_SS,c.name), h_SS.GetBinContent(bin+1), step_, min_, max_, erflag )
    for bin in xrange(nbins):
        i_SS = c.index_SS + bin
        m.mnparm(i_SS, "%i_%s"%(i_SS,c.name), template.SS.GTH_SS.GetBinContent(bin+1), step_, min_, max_, erflag)

# run MIGRAD algorithm
'... Begin MIGRAD minimization'
arglist[0]=0
m.mnexcm("MIGRAD",arglist,0,erflag)

pol = ROOT.Double(0)
polerr = ROOT.Double(0)
m.GetParameter(1,pol,polerr)
norm = ROOT.Double(0)
normerr = ROOT.Double(0)
m.GetParameter(0,norm,normerr)
if stat > 1:
    nSS = ROOT.Double(0)
    nSSerr = ROOT.Double(0)
    m.GetParameter(2,nSS,nSSerr)

print 'NORM: %.3f +/- %.3f' % (norm,normerr)
print 'POL: %.3f +/- %.3f' % (pol,polerr)
if stat > 1 and stat < 5:
    print 'NSS:  %.3f +/- %.3f' % (nSS,nSSerr)
if stat > 3:
    print '-' * 60
    print 'kW_OS = %.3f' % kW_OS_fit
    print 'kW_SS = %.3f' % kW_SS_fit
if stat == 5:
    print 'rQCD = %.3f' % rQCD_fit


#__________________________________________________________________________________________________________________
def fill_correction(compo):
    if not compo.name == 'SS':
        for j in xrange(len(compo.mu)):
            compo.corr_OS.append([])
            compo.corr_SS.append([])
            for bin in xrange(nbins):
                i_OS = compo.index_OS[j] + bin;  i_SS = c.index_SS[j] + bin 
                compo.corr_OS[j].append(get_p_value(i_OS)[0])
                compo.corr_SS[j].append(get_p_value(i_SS)[0]) 
    else:
        for bin in xrange(nbins):
            i_SS = c.index_SS + bin
            compo.corr_SS.append( get_p_value(i_SS)[0] )

def get_p_value(index):
    val = ROOT.Double(0)
    val_err = ROOT.Double(0)
    m.GetParameter(index,val,val_err)
    return float(val),float(val_err)

for c in template.components:
    fill_correction(c)

#__________________________________________________________________________________________________________________
def draw_fit_result(ratio=True):

    region = 'SR_OS'

    if stat > 1:
        NSS = ErrorFloat(float(nSS),float(nSSerr))
    else: NSS = rQCD
    NORM = float(norm)
    POL = float(pol)
    # Canvas
    cname = 'fit_%s' % stat
    c = ROOT.TCanvas(cname,cname,1500,1000)

    # Split for ratio plot
    c.Divide(1,2)
    c.cd(1).SetPad(0.0,0.2,1.0,1.0)
    c.cd(2).SetPad(0.0,0.0,1.0,0.2)
    c.cd(1)

    # data  
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.2)
    data.SetLabelSize(0.04,"y")
    if entry.xmax:
        data.SetAxisRange(0,entry.xmax)

    
    # Legend
    l = ROOT.TLegend(0.72,0.55,0.90,0.90)

    if entry.name == 'upsilon':
            l = ROOT.TLegend(0.6,0.65,0.90,0.90)
    l.SetTextSize(0.026)
    l.SetBorderSize(0)
    l.SetFillColor(0)

    #l.AddEntry(data,'"Data" (signal MC) (%.0f)' % data.Integral(0,data.GetNbinsX()+1))
    l.AddEntry(data,'"Data" (signal MC)')
    # Stack
    st = ROOT.THStack(cname,cname)
    # Same Sign
    SS = hist.get_SS(rf,entry,stream.data,region,kSS = NSS)
    SS.SetLineWidth(0)
    SS.SetLineColor(3)
    SS.SetFillColor(3)
    SS.SetMarkerStyle(0)
    SS.Rebin(rb)
    if entry.xmax:
        SS.SetAxisRange(0,entry.xmax)
    st.Add(SS)
    if stat == 1:
        SS.Scale(NORM)
    stacksum = SS.Clone('ss')

    groups = stream.MC + [stream.LH,stream.RH]

    for group in groups:

        kOS = ErrorFloat(1.0,0.0); kSS = ErrorFloat(1.0,0.0) 
        if group.name == 'Wlnu': kOS = kW_OS; kSS = kW_SS * NSS
        elif group.name == 'Zlljet': kOS = kOS*(ErrorFloat(1.,0.)-NSS)
        elif group.name == 'Ztautau_LH': 
            kOS = hist.kLH  * ErrorFloat(0.5 * (1-POL),0.)
            kSS = kOS * NSS
        elif group.name == 'Ztautau_RH': 
            kOS = hist.kRH * ErrorFloat(0.5 * (1+POL),0.)
            kSS=kOS * NSS
        elif group.name == 'Zll': kOS = kOS; kSS = kSS * NSS
        elif group.name == 'ttBar': kOS = kOS; kSS = kSS * NSS
        h = hist.get_OS_SS(rf,entry,group,region,kOS=kOS,kSS=kSS)

        h.Rebin(rb)
        if entry.xmax:
            h.SetAxisRange(0,entry.xmax)
        h.SetLineColor(entry[group.name]['fill'])
        h.SetLineWidth(0)
        h.SetMarkerStyle(0)
        h.SetFillColor(entry[group.name]['fill'])
        h.SetLabelSize(0.04,"y")

        h.Scale(NORM)
        stacksum.Add(h)
        st.Add(h)
        #l.AddEntry(h,'%s (OS-SS): %.2f +/- %.2f' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1]))
        l.AddEntry(h,'%s (OS-SS)' % group.name)
        #print '%s: %s +/- %s' % (group.name,hist.IntegralError(h)[0],hist.IntegralError(h)[1])

    #l.AddEntry(SS,'SS: %.2f +/- %.2f' % (hist.IntegralError(SS)[0],hist.IntegralError(SS)[1]))
    l.AddEntry(SS,'SS')
    data.SetXTitle(entry.xtitle)
    data.SetTitle('')

    maxi = max(data.GetMaximum(),st.GetMaximum())
    mini = 0.
    data.SetMinimum(1.2*mini)
    data.SetMaximum(1.5*maxi)

    
    data.Draw()
    st.Draw("HIST,SAME")
    data.Draw("SAME, PE")
    c.RedrawAxis()
    l.Draw()

    t = ROOT.TPaveText(0.2,0.75,0.3,0.9,"NDC")
    t.SetTextAlign(11)
    t.SetBorderSize(0)
    t.SetFillColor(0)
    t.SetTextSize(0.03)
    if comment: t.AddText(comment)
    t.AddText('POL  = %.3f +/- %.3f' % (pol,polerr))
    t.AddText('NORM = %.3f +/- %.3f' % (norm,normerr))
    if stat>1:
        t.AddText('NSS  = %.3f +/- %.3f' % (nSS,nSSerr))

    t.Draw()

    if ratio:
        c.cd(2)

        dataratio = data.Clone('ratio')
        if entry.xmax:
            dataratio.SetAxisRange(0,entry.xmax)
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
    c.Print('%s/%s.gif+1' % (s.plotdir,cname))

def draw_correction(compo,j):
    # Canvas
    cname = 'corr_%s' % compo.name
    c = ROOT.TCanvas(cname,cname,1500,1000)

    # Split for ratio plot
    c.Divide(1,2)
    c.cd(1).SetPad(0.0,0.2,1.0,1.0)
    c.cd(2).SetPad(0.0,0.0,1.0,0.2)
    c.cd(1)

    nominal = compo.TH_OS_UNSCALED[j].Clone()
    fitted  = nominal.Clone('fit')
    fitted.Reset()

    for bin in xrange(fitted.GetNbinsX()):
        fitted.SetBinContent(bin+1,compo.corr_OS[j][bin])

    #fitted.Scale(compo.mu[j])
    fitted.SetLineColor(2)

    nominal.SetTitle('%s_%s' % (compo.name,j))
    nominal.SetMinimum(0)
    nominal.SetMaximum(1.2 * max(nominal.GetMaximum(), fitted.GetMaximum()))
    nominal.Draw()
    fitted.Draw("SAME,L")  

    l = ROOT.TLegend(0.82,0.7,0.90,0.90)
    l.SetTextSize(0.026)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    l.AddEntry(nominal,'nominal')
    l.AddEntry(fitted,'fitted')
    l.Draw()

    c.Print("%s/%s_%s.gif+1" % (s.plotdir,cname,stat))

def draw_sum(compo):
    cname = 'corr_sum_%s' % compo.name
    c = ROOT.TCanvas(cname,cname,1500,1000)

    # Split for ratio plot
    c.Divide(1,2)
    c.cd(1).SetPad(0.0,0.2,1.0,1.0)
    c.cd(2).SetPad(0.0,0.0,1.0,0.2)
    c.cd(1)

    nominal = compo.TH_OS_UNSCALED[0].Clone()
    for j in xrange(1,6):
        nominal.Add(compo.TH_OS_UNSCALED[j])
    fitted  = compo.TH_OS_UNSCALED[0].Clone('fit')
    fitted.Reset()

    for j in xrange(6):
        for bin in xrange(fitted.GetNbinsX()):
            fitted.AddBinContent(bin+1,compo.corr_OS[j][bin])

    #fitted.Scale(compo.mu[j])
    fitted.SetLineColor(2)

    nominal.SetTitle('%s_%s' % (compo.name,j))
    nominal.SetMinimum(0)
    nominal.SetMaximum(1.2 * max(nominal.GetMaximum(), fitted.GetMaximum()))
    nominal.Draw()
    fitted.Draw("SAME,L")  

    l = ROOT.TLegend(0.82,0.7,0.90,0.90)
    l.SetTextSize(0.026)
    l.SetBorderSize(0)
    l.SetFillColor(0)
    l.SetHeader(compo.name)
    l.AddEntry(nominal,'nominal')
    l.AddEntry(fitted,'fitted',"L")
    l.Draw()

    c.Print("%s/%s_%s.gif+1" % (s.plotdir,cname,stat))

draw_fit_result()


if stat>2:
    draw_sum(template.LH)
    draw_sum(template.RH)
    draw_sum(template.Wlnu)        
    draw_sum(template.Zll)
    draw_correction(template.ttbar,0)

