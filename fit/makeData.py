import ROOT
import run,hist, Groups, Drawer
from ErrorFloat import ErrorFloat
import math
import sys

spec = run.el_fit
mu = spec.stream
entry = spec.reg_entries[0]
k_entry = Drawer.upsilon_CR
groups = mu.MC + [mu.signal]
verbose = False

rf = ROOT.TFile('data/%s_NoZpt.root'%mu.name, 'UPDATE')
kW_OS = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, groups = groups)
kW_SS = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose, groups = groups)
rQCD = hist.calc_rQCD(rf,mu,k_entry,'QCD',kW_OS = kW_OS, kW_SS = kW_SS, verbose = verbose, groups = groups)

print kW_OS, kW_SS,rQCD

def addQCDest(region = 'SR_SS'):
    qcd = hist.get_QCDW(rf,mu,entry,'SR_SS', verbose = False, groups = groups)
    qcd.SetName('h_%s_QCD_%s' % (entry.name,region))
    print qcd.GetName(),qcd.Integral(), rQCD.val*qcd.Integral()
    qcd.Write()

def addQCDCR( region = 'QCD_SS', sys = ''):
    if not sys:
        qcd = hist.get_QCD(rf,mu,k_entry,region,kW=kW_SS, groups = groups)
        print qcd.Integral()
        qcd.SetName('h_upsilon_CR_QCD_QCD_SS')
        qcd.Write()

def addSRdata( region = 'SR_OS', pol = ''):
    # first get QCD estimate with data driven W 
    if pol:
        groups = mu.MC + [getattr(mu,pol)]
    else:
        groups = mu.MC + [mu.signal]

    data = hist.get_QCDW(rf,mu,entry,'SR_SS', rqcd= rQCD,verbose = False, groups = groups) 
    for group in groups:
        if group.name == 'Wlnu': 
            h = hist.get_dataW(rf, mu, entry, 'WCR_OS', verbose = False, groups = groups)            
        else:
            h = hist.get_OS(rf,entry,group,region)
        data.Add(h)
    data.SetName('h_upsilon_SRData%s_SR_OS' % pol)
    for bin in xrange(data.GetNbinsX()+1):
        data.SetBinError(bin, math.sqrt(data.GetBinContent(bin)))
    print ErrorFloat(*hist.IntegralError(data))
    #check similar errors
    #d1 = rf.Get('h_visMass_Muons_SR_OS')
    #print ErrorFloat(*hist.IntegralError(d1))
    data.Write()

def addWSR( region = 'OS', sys = ''):

    h = hist.get_dataW(rf, mu, entry, 'WCR_%s' % region, verbose = False, groups = groups, sf = rf, sysf=sys) 
    if not sys:
        h.SetName('h_upsilon_WData_SR_%s' % region )
    else:
        h.SetName('h_upsilon_WData_SR_%s_%s' % (region,sys) )
    #print h.GetName(), h.Integral()
    h.Write()


def addWCR( region = 'WCR_OS', sys = ''):
    if not sys:
        h = rf.Get('h_upsilon_CR_Wlnu_%s' % region)
        if 'OS' in region:
            h.Scale(kW_OS.val)
        elif 'SS' in region:
            h.Scale(kW_SS.val)
        h.SetName('h_upsilon_CR_WData_%s' % region)
        print region, h.Integral()
        h.Write()

    else:
        if 'OS' in region:
            wregion = 'WCR_OS'
        else:
            wregion = 'WCR_SS'
        if 'SR' in region:
            sentry = entry
            #print sentry.name
        else:
            sentry = None

        h = hist.get_nWCR( rf, mu, k_entry, wregion, sys = sys)
        rW = hist.calc_rW(rf,mu,k_entry, region, sys = sys , sentry = sentry)
        h.Scale(rW.val)
        if 'SR' in region:
            h.SetName( 'h_upsilon_WData_%s_%s' % (region,sys)  )
        else:
            h.SetName( 'h_upsilon_CR_WData_%s_%s' % (region,sys)  )
        #print h.GetName(),h.Integral()
        h.Write()



def addWsys( mu, region):
    # nominal W estimate from WCR
    h = hist.get_dataW(rf, mu, entry, 'WCR_%s' % region, verbose = False, groups = groups)
    ## get slope
    if mu.name == 'mu':
        if 'OS' in region:
            p = -0.0269     
        elif 'SS' in region:
            p = 0.0930
    if mu.name == 'el':
        if 'OS' in region:
            p = -0.0494
        elif 'SS' in region:
            p = 0.0626
    hup = h.Clone('up')
    hup.SetName('h_upsilon_WData_SR_%s_Wshape_UP' % (region))
    hdown = h.Clone('down')
    hdown.SetName('h_upsilon_WData_SR_%s_Wshape_DOWN' % (region))
    for bin in xrange(hup.GetNbinsX()+1):
        bc = hup.GetBinCenter(bin)
        hup.SetBinContent(bin, (1 + (bc-0.25) * p) * hup.GetBinContent(bin) )
        hdown.SetBinContent(bin, (1. -  (bc-0.25) * p) * hdown.GetBinContent(bin) )        
        #print hup.GetBinContent(bin)

    hup.Write()
    hdown.Write()
    #print hup.Integral()
    #print hdown.Integral()


####################################################################################################################
## Add Data Estimate
addSRdata()
#addSRdata(pol='LH')
#addSRdata(pol='RH')
## Add Nominal QCD estimate
addQCDest()
## Add SS QCD est
addQCDCR()
## ADD data Driven W in SR
addWSR(region = 'OS')
addWSR(region = 'SS')
## Add Normalized W expected in Control Regions
addWCR( region = 'WCR_OS')
addWCR(region = 'WCR_SS')
addWCR( region = 'QCD_OS')
addWCR( region = 'QCD_SS')

makesys=False
if makesys:
    tausys = ['tau_id_sys','tau_id_stat', 'tau_el', 'tau_fake']
    lepsys = ['el_id','el_iso','mu_id','mu_iso']
    lepES = ['ElES_LowPt','ElES_PS','ElES_R12','ElES_Zee','ElEnResSys', 'MuSys' ]
    JES = ['JER','JES_BJet','JES_Detector1','JES_EtaMethod',
                         'JES_EtaModelling','JES_FlavComp','JES_FlavResp','JES_Modelling1','JES_PUMu','JES_PUNPV','JES_PURho',
                         'JES_Statistical1','JVF','METResSys','METScaleSys']
    PU = ['PU_rescaling']
    eltrig = ['el_trig']
    mutrig = ['mu_trig']
    TES = ['TES_FINAL']

    mu_syst = tausys + lepsys + JES + TES + ['MuSys']
    el_syst = tausys + lepsys + lepES + JES + TES

    addWsys(mu,'OS')
    addWsys(mu,'SS')

    for sys in el_syst:
        for ud in ['_UP','_DOWN']:
            sv = sys+ud
            if 'TES' in sys:
                addWSR( region = 'OS', sys = sv) # this is for W estimate, only TES FINAL
                addWSR( region = 'SS', sys = sv)
            #addWCR
            addWCR( region = 'SR_OS', sys = sv)
            addWCR( region = 'SR_SS', sys = sv)        
            addWCR( region = 'WCR_OS', sys = sv)
            addWCR( region = 'WCR_SS', sys = sv)
            addWCR( region = 'QCD_OS', sys = sv)
            addWCR( region = 'QCD_SS', sys = sv)





