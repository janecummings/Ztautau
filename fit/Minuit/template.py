import ROOT
import Groups
import hist
from ErrorFloat import ErrorFloat
import Drawer

entry = Drawer.upsilon_fit

p = -0.15
int_lumi = 20276.9
## acceptance normalization
kLH = ErrorFloat(1/(0.5*(1-p)),0.)
kRH = ErrorFloat(1/(0.5*(1+p)),0.)

class template:
    def __init__(self,rf,kCR,stream,fitstat=1,rb=5):
        self.rf = rf
        self.kCR = kCR
        self.stat = fitstat
        self.stream = stream
        self.LH = None
        self.RH = None
        self.rb = rb
        self.LH = component('Ztautau_LH', stream.LH, fitstat)
        self.RH = component('Ztautau_RH', stream.RH, fitstat)

        self.Wlnu   = component('Wlnu', stream.Wlnu, fitstat)
        self.Zll    = component('Zll', stream.Zll, fitstat)
        self.ttbar  = component('ttbar', stream.ttbar, fitstat)

        self.SS     = component('SS', stream.data, 0)

        self.components = [self.LH,self.RH,
                           self.Wlnu,self.Zll,self.ttbar,
                           self.SS]
        self.WCR_OS = None
        self.WCR_SS = None
        self.QCD_OS = None
        self.QCD_SS = None

    def fill(self):
        if self.stat in [4,5]:
            h_OS = hist.get_OS(self.kCR,entry,self.stream.data,'WCR_OS')
            h_SS = hist.get_SS(self.kCR,entry,self.stream.data,'WCR_SS') 
            self.WCR_OS = h_OS.Integral(0,h_OS.GetNbinsX()+1)
            self.WCR_SS = h_SS.Integral(0,h_SS.GetNbinsX()+1)
            for c in self.components:
                c.fill(self.rf,self.kCR)
                c.rebin(self.rb)
            if self.stat == 5:
                h_OS = hist.get_OS(self.kCR,entry,self.stream.data,'QCD_OS')
                h_SS = hist.get_SS(self.kCR,entry,self.stream.data,'QCD_SS') 
                self.QCD_OS = h_OS.Integral(0,h_OS.GetNbinsX()+1)
                self.QCD_SS = h_SS.Integral(0,h_SS.GetNbinsX()+1)
        else:
            for c in self.components:
                c.fill(self.rf,self.kCR)
                c.rebin(self.rb)

    def load(self):
        for c in self.components:
            c.load()
            for (xsec,total) in zip(c.cross_section,c.totalEvents):
                c.mu.append( xsec * int_lumi / total)

    def index(self,offset,nbins):
        m = offset
        for c in self.components:
            if c.name == 'SS': continue
            c.index_OS = []
            c.index_SS = []
            nsamples = len(c.mu)
            for n in xrange(nsamples):
                #c.index_OS.append( m + n*nbins) # OS ONLY
                c.index_OS.append( m + 2*n*nbins)
                c.index_SS.append( m + (2*n+1)*nbins)
            #m += nsamples*nbins #OS ONLY
            m += nsamples*2*nbins 
        self.SS.index_SS = m

class component:
    def __init__(self,name,group,fitstat,region='SR_OS',nsamples = 1):
        
        self.name = name
        self.group = group
        self.region = region
        self.stat=fitstat
        # scaled group histogram from output
        self.GTH_OS = None
        self.GTH_SS = None
        # list of unscaled sample hist from output
        self.TH_OS_UNSCALED = []
        self.TH_SS_UNSCALED = []
        # list of sample cross-sections and totalEvents; filled in load
        self.cross_section = []
        self.totalEvents = []
        # mu scale factors in fit
        self.mu = []
        # WCR hist
        self.WCR_OS = None
        self.WCR_SS = None
        # scale factors  -- maybe i dont need this 
        self.kOS = ErrorFloat(1.,0.)
        self.kSS = ErrorFloat(1.,0.)
        # indices
        self.index_OS = None
        self.index_SS = None
        # hold corrections
        self.corr_OS = []
        self.corr_SS = []

    def load(self):
        # load totalEvents and cross sections
        if self.stat in [0,1,2]: return 0
        for sample in self.group:
            sample.load()
            self.cross_section.append(sample.xsec)
            self.totalEvents.append(sample.totalEvents)

            #del sample.chain

    def fill(self,rf,kCR):
        # fill ther relevant template components for the fit
        if self.stat == 0:
            self.GTH_SS = hist.get_SS(rf,entry,self.group,self.region)
        if self.stat in [1,2,3,4,5]:
            self.GTH_OS = hist.get_OS(rf,entry,self.group,self.region) # no scaling
            self.GTH_SS = hist.get_SS(rf,entry,self.group,self.region) # no scaling
        if self.stat in [3,4,5]:
            #for (sample, xsec, total) in zip(self.group,self.cross_section,self.totalEvents):
            for sample in self.group:
                region = self.region + '_%s_unscaled' % sample.name
                self.TH_OS_UNSCALED.append(hist.get_OS(rf,entry,self.group,region))
                self.TH_SS_UNSCALED.append(hist.get_SS(rf,entry,self.group,region))
            #    self.mu_OS.append(1.) ##
        if self.stat in [4,5]:
            h_OS = hist.get_OS(kCR,entry,self.group,'WCR_OS')
            h_SS = hist.get_SS(kCR,entry,self.group,'WCR_SS')
            self.WCR_OS = h_OS.Integral(0,h_OS.GetNbinsX()+1)
            self.WCR_SS = h_SS.Integral(0,h_SS.GetNbinsX()+1)
        if self.stat == 5:
            h_OS = hist.get_OS(kCR,entry,self.group,'QCD_OS')
            h_SS = hist.get_SS(kCR,entry,self.group,'QCD_SS')
            self.QCD_OS = h_OS.Integral(0,h_OS.GetNbinsX()+1)
            self.QCD_SS = h_SS.Integral(0,h_SS.GetNbinsX()+1)

    def rebin(self,rb):
        if self.stat == 0:
            self.GTH_SS.Rebin(rb)
        if self.stat in [1,2,3]:
            for h in [self.GTH_OS,self.GTH_SS]:
                h.Rebin(rb)
        if self.stat in [3,4,5]:
            for h in self.TH_OS_UNSCALED + self.TH_SS_UNSCALED:
                if h: h.Rebin(rb)

    def scale(self,kX):
        if self.stat == 0:
            self.GTH_SS.Scale(kX)
        else:
            for h in [self.GTH_OS,self.GTH_SS]:
                h.Scale(kX)
        if self.stat in [4,5]:
            self.WCR_OS = self.WCR_OS * kX
            self.WCR_SS = self.WCR_SS * kX
        return 0

def fake_data(rf,stream,kW_OS,kW_SS,rQCD,rb):
    region = 'SR_OS'
    fake_data = hist.get_SS(rf,entry,stream.data,region,kSS=rQCD).Clone('data')
    fake_data.Rebin(rb)
    groups = [stream.LH,stream.RH] + stream.MC
    for group in groups:
        print group.name
        kOS = ErrorFloat(1.0,0.0); kSS = ErrorFloat(1.0,0.0) 
        if group.name == 'Wlnu': kOS = kW_OS; kSS = kW_SS * rQCD
        elif group.name == 'ttBar': kOS = kOS; kSS = kSS * rQCD
        elif group.name == 'Zll': kOS = kOS; kSS = kSS * rQCD
        elif group.name == 'Ztautau_LH': kOS = kOS; kSS = kOS * rQCD 
        elif group.name == 'Ztautau_RH': kOS = kOS; kSS = kOS * rQCD 
        h = hist.get_OS_SS(rf,entry,group,region,kOS=kOS,kSS=kSS)
        h.Rebin(rb)
        fake_data.Add(h)

    return fake_data





