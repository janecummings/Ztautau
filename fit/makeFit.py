import sys
from ROOT import *
import Groups, run, streams, Drawer
import hist as drawhist


fit_ehad = True
fit_muhad = False
statOn = False
qcdOn = False
lumiOff = False
tes = False
norm_sys = False
Wsys = False

name = 'FitTest'
if fit_ehad: fname = 'ehad'
if fit_muhad: fname = 'muhad'
if (fit_ehad and fit_muhad): fname = 'lephad'
name += '_%s' % fname
if statOn:
	name += '_stat'
if not qcdOn:
    name += '_noQCDfit'
if tes:
    name+='_TES_'
fitPol = ''
if fitPol:
    name+='_%s' % fitPol
if Wsys:
    name+='_Wsys'
if norm_sys:
    name+='_norm'
## Compare ZpT weights 
name += '_NoZpt'

s_entry = Drawer.upsilon_20
k_entry = Drawer.upsilon_CR
################################################################
############# Make Measurement #################################
meas = RooStats.HistFactory.Measurement(name, name)
meas.SetOutputFilePrefix('TestResults/%s' % name)
meas.SetBinLow(0)
meas.SetBinHigh(20)
meas.SetExportOnly(True)
################################################################
############# Set Constants and Preprocess #####################
meas.SetPOI("Ptau")
meas.SetLumi( 1.0 )
meas.SetLumiRelErr( 0.028 )
ptau = 0
minPtau = -2
maxPtau = 2
meas.AddPreprocessFunction("Ptau_LH","0.500 * (1-Ptau)","Ptau[%f,%f]" %(minPtau,maxPtau))
meas.AddPreprocessFunction("Ptau_RH","0.500 * (1+Ptau)","Ptau[%f,%f]" %(minPtau,maxPtau))
################################################################
############# Configure ########################################

relError = 0.01
constraintType = 'Poisson'

###################################################################################
###################################################################################
class FitSample():
    def __init__(self,group,name=None,
                 theory = True,
                 normFactor = None,
                 shapeFactor = None,
                 hand = None,
                 region = None,
                 stat = True,
                 overallSys = None):
        if group:
            self.group = group
            self.name = group.name
            self.cname = name or group.name
        else:
            self.group = None
            self.name = name
            self.cname = name

        self.shapeFactor = shapeFactor
        self.normFactor = normFactor
        self.hand = hand
        self.theory = theory
        self.stat = stat
        self.overallSys = overallSys
class NormFactor():
    def __init__(self,name,val,mini,maxi):
        self.name = name
        self.val = val
        self.mini = mini
        self.maxi = maxi
        #self.const = const
        #self.sample = sample
        #self.regions = regions 
    def pstr(self):
        return '%s = (%s,%s,%s) %s' % ( self.name,self.val,self.mini,self.maxi,self.const)
class OverallSys():
    def __init__(self,name,low,high):
        self.name = name
        self.low = low
        self.high = high
class mode():
    def __init__(self, name, data, inputFile, sf = None, osyst = [], hsyst = [], stream = None):
        self.name = name
        self.inputFile = inputFile
        self.inputroot = TFile(inputFile)
        self.data = data
        self.systFile = sf
        self.systroot = TFile(sf)
        self.osyst = osyst
        self.hsyst = hsyst
        self.stream = stream
def addSample(channel, mode, component, region, stat = False):
    entry = 'upsilon'
    if 'SR' not in region: entry += '_CR'
    hist = 'h_%s_%s_%s' % ( entry, component.name, region )
    cname = component.cname + '_' + region + "_" + mode.name
    sample = RooStats.HistFactory.Sample(cname, hist, mode.inputFile)
    if component.shapeFactor: stat=False
    if not component.theory or lumiOff:
        sample.SetNormalizeByTheory(False)
    elif stat and statOn: sample.ActivateStatError()

    if component.normFactor:
        for nF in component.normFactor:
            if 'Ztautau' in component.name: 
                nfname = nF.name
            else:
                nfname = nF.name + '_' + mode.name
            if 'OS' in nF.name and 'SS' in region: continue
            if 'SS' in nF.name and 'OS' in region: continue
            sample.AddNormFactor( nfname, nF.val, nF.mini, nF.maxi)
    if statOn and 'W' in component.name:
        if 'SR' in region:
            rW =  drawhist.calc_rW(mode.inputroot,mode.stream,s_entry,region )
        elif 'QCD' in region: 
            rW =  drawhist.calc_rW(mode.inputroot,mode.stream,k_entry,region )
        else:
            rW = None
        if rW:
            rW = rW/rW
            sample.AddOverallSys( 'alpha_%s_rW_%s' % (mode.name,region) , 1 - rW.err, 1 + rW.err)

    if Wsys and 'W' in component.name and 'SR' in region:
        sample.AddHistoSys('Wshape_%s' % region, '%s_Wshape_DOWN' % hist,  mode.inputFile, '', '%s_Wshape_UP' % hist,  mode.inputFile, ''  )

    if mode.hsyst:
        for hs in mode.hsyst:
            #if 'W' in component.name: continue
            sample.AddHistoSys(hs, '%s_%s_DOWN' % (hist,hs), mode.inputFile, '',  '%s_%s_UP' % (hist,hs),mode.inputFile, '') 


    if mode.osyst:
        for os in mode.osyst:
            
            mode.inputroot.cd()
            print '%s_%s_DOWN' % (hist,os) 
            os_nom = mode.inputroot.Get(hist).Integral()
            os_down = mode.inputroot.Get('%s_%s_DOWN' % (hist,os) ).Integral()
            os_up =mode.inputroot.Get('%s_%s_UP' % (hist,os)).Integral()
            ## Don't include 
            if os_nom - os_down == 0 and os_nom - os_up == 0:
                continue
            if os_nom:
                os_low = os_down/os_nom
                os_high = os_up/os_nom
                sample.AddOverallSys( os, os_low, os_high )
    channel.AddSample(sample)


def addMode( mode, samples ):
    #################################
    # Signal Region
    #################################   
    SR_OS = RooStats.HistFactory.Channel( "SR_OS" + "_" + mode.name )
    SR_OS.SetData( 'h_upsilon_SRData%s_SR_OS' % fitPol, mode.inputFile)
    SR_OS.SetStatErrorConfig( relError, constraintType )

    for component in samples:
        addSample( SR_OS, mode, component, 'SR_OS', stat = statOn )

    qcd_SR_OS = RooStats.HistFactory.Sample( "QCD_SR_OS" + "_" + mode.name, "h_upsilon_QCD_SR_SS", mode.inputFile)
    qcd_SR_OS.SetNormalizeByTheory(False)
    if qcdOn:
        qcd_SR_OS.AddShapeFactor("QCDSHAPE_%s" % mode.name)
    qcd_SR_OS.AddNormFactor( "rQCD_%s" % mode.name, 1, 0., 2 )

    SR_OS.AddSample(qcd_SR_OS)
    meas.AddChannel(SR_OS)

    #################################
    # Same Sign
    #################################   

    SR_SS = RooStats.HistFactory.Channel( "SR_SS" + "_" + mode.name )
    SR_SS.SetData( 'h_upsilon_%s_SR_SS' % mode.data, mode.inputFile)
    SR_SS.SetStatErrorConfig( relError, constraintType )
    for component in samples:
        addSample( SR_SS, mode, component, 'SR_SS', stat = statOn )

    qcd_SR_SS=RooStats.HistFactory.Sample( "QCD_SR_SS_%s" % mode.name, "h_upsilon_QCD_SR_SS", mode.inputFile)
    if qcdOn:
        qcd_SR_SS.AddShapeFactor("QCDSHAPE_%s" % mode.name)
    qcd_SR_SS.SetNormalizeByTheory(False)

    SR_SS.AddSample(qcd_SR_SS)
    meas.AddChannel(SR_SS)

    #################################
    # OS W Control Region
    #################################
    WCR_OS = RooStats.HistFactory.Channel( "WCR_OS_%s" % mode.name )
    WCR_OS.SetData( 'h_upsilon_CR_%s_WCR_OS' % mode.data, mode.inputFile )
    WCR_OS.SetStatErrorConfig( relError, constraintType )
    for component in samples:
        addSample( WCR_OS, mode, component, 'WCR_OS', stat = False )
        #addSample( WCR_OS, mode, component, 'WCR_OS', stat = statOn )
    meas.AddChannel(WCR_OS)

    #################################
    # SS W Control Region
    #################################
    WCR_SS = RooStats.HistFactory.Channel( "WCR_SS_%s" % mode.name )
    WCR_SS.SetData( 'h_upsilon_CR_%s_WCR_SS' % mode.data, mode.inputFile )
    WCR_SS.SetStatErrorConfig( relError, constraintType )
    for component in samples:
    	addSample( WCR_SS, mode, component, 'WCR_SS', stat = False )
        #addSample( WCR_SS, mode, component, 'WCR_SS', stat = statOn )
    meas.AddChannel(WCR_SS)

    #################################
    # OS QCD Control Region
    #################################
    QCD_OS =RooStats.HistFactory.Channel( "QCD_OS_%s" % mode.name )
    QCD_OS.SetData( 'h_upsilon_CR_%s_QCD_OS' % mode.data, mode.inputFile )
    QCD_OS.SetStatErrorConfig( relError, constraintType )
    for component in samples:
        #addSample( QCD_OS, mode, component, 'QCD_OS', stat = statOn )
        addSample( QCD_OS, mode, component, 'QCD_OS', stat = False )

    qcd_QCD_OS =RooStats.HistFactory.Sample( "QCD_QCD_OS_%s" % mode.name, "h_upsilon_CR_QCD_QCD_SS", mode.inputFile)
    qcd_QCD_OS.SetNormalizeByTheory(False)
    qcd_QCD_OS.AddNormFactor( "NSS_%s" % mode.name, 1. , 0., 10 )
    qcd_QCD_OS.AddNormFactor( "rQCD_%s" % mode.name, 1 , 0., 10 )

    QCD_OS.AddSample(qcd_QCD_OS)
    meas.AddChannel(QCD_OS) 

    #################################
    # OS QCD Control Region
    #################################
    QCD_SS =RooStats.HistFactory.Channel( "QCD_SS_%s" % mode.name )
    QCD_SS.SetData( 'h_upsilon_CR_%s_QCD_SS' % mode.data, mode.inputFile )
    QCD_SS.SetStatErrorConfig( relError, constraintType )
    for component in samples:
        #addSample( QCD_SS, mode, component, 'QCD_SS', stat = statOn )
        addSample( QCD_SS, mode, component, 'QCD_SS', stat = False )

    qcd_QCD_SS =RooStats.HistFactory.Sample( "QCD_QCD_SS_%s" % mode.name, "h_upsilon_CR_QCD_QCD_SS", mode.inputFile)
    qcd_QCD_SS.SetNormalizeByTheory(False)
    qcd_QCD_SS.AddNormFactor( "NSS_%s" % mode.name, 1. , 0., 10 )

    QCD_SS.AddSample(qcd_QCD_SS)
    meas.AddChannel(QCD_SS) 


##############################
# Define NormFactors

kW_OS = NormFactor( 'kW_OS', 1, 0, 2 )
kW_SS = NormFactor( 'kW_SS', 1, 0, 2 )
ZtautauNorm = NormFactor("ZtautauNorm", 1.0, 0., 10.)
Ptau_LH = NormFactor("Ptau_LH", 0.5*(1 - ptau),0.5*(1- maxPtau),0.5*(1-minPtau))
Ptau_RH = NormFactor("Ptau_RH", 0.5*(1 + ptau),0.5*(1 + minPtau ),0.5*(1 + maxPtau))

##############################
# Define FitSamples

fit_LH = FitSample(Groups.LeftZtautau, name = 'LH', normFactor = [ZtautauNorm, Ptau_LH])
fit_RH =  FitSample(Groups.RightZtautau, name = 'RH', normFactor = [ZtautauNorm, Ptau_RH])
#fit_LH = FitSample(Groups.LeftZtautau, name = 'LH', normFactor = [Ptau_LH])
#fit_RH =  FitSample(Groups.RightZtautau, name = 'RH', normFactor = [Ptau_RH])
fit_Wlnu = FitSample(None, name = 'WData',  normFactor = [kW_OS,kW_SS], theory = False)
fit_Zlljet = FitSample(Groups.Zlljet)
fit_Zltau = FitSample(Groups.Zlltau)
fit_ttbar = FitSample(Groups.ttBar, name = 'tt')

fit_samples = [fit_LH, fit_RH, fit_Wlnu, fit_Zlljet, fit_Zltau, fit_ttbar]

# Systematics
tausys = ['tau_id_sys','tau_id_stat', 'tau_el', 'tau_fake']
lepsys = ['el_id','el_iso','mu_id','mu_iso']
muMS = ['MuSys']
lepES = ['ElES_LowPt','ElES_PS','ElES_R12','ElES_Zee','ElEnResSys', 'MuSys' ]
JES = ['JER','JES_BJet','JES_Detector1','JES_EtaMethod',
                     'JES_EtaModelling','JES_FlavComp','JES_FlavResp','JES_Modelling1','JES_PUMu','JES_PUNPV','JES_PURho',
                     'JES_Statistical1','JVF','METResSys','METScaleSys']
eltrig = ['el_trig']
mutrig = ['mu_trig']


if norm_sys:
    overall_mu = mutrig + tausys + lepsys + JES + muMS#['mu_trig']
    overall_el =  tausys + lepsys + lepES + JES
else:
    overall_el = []
    overall_mu = []
if tes:
    shape = ['TES_FINAL']
else:
    shape = []

#####################
# Build Model

el = mode( 'ehad', 'Egamma', 'data/el_NoZpt.root', sf = 'data/el_syst.root', osyst = overall_el, hsyst = shape, stream = streams.el)
mu = mode( 'muhad', 'Muons','data/mu_regions.root', sf = 'data/mu_regions.root', osyst = overall_mu, hsyst = shape, stream = streams.mu)

if fit_ehad: addMode( el , fit_samples)
if fit_muhad: addMode( mu , fit_samples)

meas.CollectHistograms()
meas.PrintTree()
meas.PrintXML("xml/%s" % (name))

RooStats.HistFactory.MakeModelAndMeasurementFast( meas )
