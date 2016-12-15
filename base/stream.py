"""
Region Definitions
"""
from Sample import Cut,Region
from Groups import *
from ErrorFloat import ErrorFloat

slt_ = Cut('') # placeholder for stream class
slt_el = Cut('evtsel_is_el',name='evtsel_is_el')
slt_mu = Cut('evtsel_is_mu',name='evtsel_is_mu')

class stream:
    def __init__(self,
                 name='',
                 data=None,
                 signal=None,
                 LH = None,
                 RH = None,
                 Wlnu = None,
                 ttbar = None,
                 Zltau = None,
                 Zlljet = None,
                 Zll = None,
                 slt_ = None,
                 weight = None
                 ):
        self.name = name
        self.weight = weight
        self.data = data
        self.signal = signal
        self.LH = LH
        self.RH = RH
        self.Wlnu = Wlnu
        self.ttbar = ttbar 
        self.Zltau = Zltau
        self.Zlljet = Zlljet
        self.Zll = Zll
        self.groups = []
        for group in [data,signal,LH,RH,Wlnu,ttbar,Zltau,Zlljet]:
            if group is None: continue
            if isinstance(group,list): self.groups.extend(group)
            else: self.groups.append(group)
        if not slt_:
            if name == 'el': slt_ = slt_el
            elif name == 'mu': slt_ = slt_mu
        self.slt_ = slt_

    def get_selection(self,cuts = [], f=0, name=''):
        cuts = [self.slt_ if cut is slt_ else cut for cut in cuts] 
        select = []
        for i,cut in enumerate(cuts):
            if f>i: continue
            if name: r_name = '%i_%s_%s' % (i,name,cut.name)
            else: r_name = '%i_%s' % (i,cut.name)
            self.make_single_region(cuts=cuts[:i+1], name = r_name)
            s = Region(r_name,cutset = cuts[:i+1])
            select.append(s)
        return select

    # make selection for cuts 
    def make_selection(self,cuts = [], f=0, name = 'select'):
        cuts = [self.slt_ if cut is slt_ else cut for cut in cuts]
        setattr(self,name,self.get_selection(cuts = cuts, f=f)) 

    # make OS and SS regions 
    def make_region(self,cuts = [], S_ = [OS_,SS_], name = 'SR'):
        cuts = [self.slt_ if cut is slt_ else cut for cut in cuts] 
        region = [Region('%s_%s' % (name,s_.name),cutset = cuts + [s_]) for s_ in S_]
        setattr(self,name,region)
    # make single region
    def make_single_region(self,cuts = [], name = 'region'):
        cuts = [self.slt_ if cut is slt_ else cut for cut in cuts]
        region = [Region('%s' % (name),cutset = cuts)]
        setattr(self,name,region)


##########################
# basic stream definitions 
##########################
class elstream(stream):
    def __init__(self):
        self.name = 'el'
        self.data = Egamma
        self.signal = Ztautau
        self.LH = LeftZtautau
        self.RH = RightZtautau
        self.ttbar = ttBar
        self.Zll = Zll
        self.Zltau = Zlltau
        self.Zlljet = Zlljet
        self.MC = [ttBar,Zlltau,Zlljet,Wlnu]
        self.Wlnu = Wlnu
        self.groups = [self.data,self.signal,self.LH,self.RH] + self.MC
        self.slt_ = slt_el
        self.weight = None


class eltheory(stream):
    def __init__(self):
        self.name = 'el'
        self.data = Egamma
        self.signal = Ztautau
        self.LH = LeftZtautau
        self.RH = RightZtautau
        self.ttbar = ttBar
        self.Zll = Zll
        self.Zltau = Zlltau
        self.Zlljet = Zlljet
        self.MC = [LeftAlpgen, RightAlpgen, LeftPythia, RightPythia, LeftPowheg, RightPowheg, LeftHerwig, RightHerwig]
        self.Pythia = [LeftPythia,RightPythia]
        self.Alpgen = [LeftAlpgen_RW,RightAlpgen_RW]        
        self.Wlnu = Wlnu
        self.groups = self.MC
        self.slt_ = slt_el
        self.weight = None

class mutheory(stream):
    def __init__(self):
        self.name = 'mu'
        self.data = Muons
        self.signal = Ztautau
        self.LH = LeftZtautau
        self.RH = RightZtautau
        self.ttbar = ttBar
        self.Zll = Zll
        self.Zltau = Zlltau
        self.Zlljet = Zlljet
        self.MC = [LeftAlpgen, RightAlpgen, LeftPythia, RightPythia, LeftPowheg, RightPowheg, LeftHerwig, RightHerwig]
        self.Pythia = [LeftPythia,RightPythia]
        self.Alpgen = [LeftAlpgen_RW,RightAlpgen_RW]
        self.Wlnu = Wlnu
        self.groups = self.MC
        self.slt_ = slt_mu
        self.weight = None


class mu_modes(stream):
    def __init__(self):
        self.name = 'mu'
        self.LH = LeftZtautau
        self.RH = RightZtautau
        self.LH_NP = LHNP 
        self.LH_RHO = LHRHO
        self.LH_PINU = LHPINU 
        self.LH_SPX = LHSPX 

        self.RH_NP = RHNP
        self.RH_RHO = RHRHO
        self.RH_PINU = RHPINU 
        self.RH_SPX = RHSPX 

        self.groups = [self.LH,self.RH,self.LH_NP,self.LH_SPX, self.LH_PINU, self.LH_RHO, self.RH_NP,self.RH_SPX, self.RH_PINU, self.RH_RHO]
        self.slt_ = slt_mu
        self.weight = None

class el_modes(stream):
    def __init__(self):
        self.name = 'el'
        self.LH = LeftZtautau
        self.RH = RightZtautau
        self.LH_NP = LHNP 
        self.LH_RHO = LHRHO
        self.LH_PINU = LHPINU 
        self.LH_SPX = LHSPX 

        self.RH_NP = RHNP
        self.RH_RHO = RHRHO
        self.RH_PINU = RHPINU 
        self.RH_SPX = RHSPX 

        self.groups = [self.LH,self.RH,self.LH_NP,self.LH_SPX, self.LH_PINU, self.LH_RHO, self.RH_NP,self.RH_SPX, self.RH_PINU, self.RH_RHO]
        self.slt_ = slt_el
        self.weight = None

class mustream(stream):
    def __init__(self):
        self.name = 'mu'
        self.data = Muons
        self.signal = Ztautau
        self.LH = LeftZtautau
        self.RH = RightZtautau
        self.ttbar = ttBar
        self.Zll = Zll
        self.Zltau = Zlltau
        self.Zlljet = Zlljet
        self.MC = [ttBar,Zlltau,Zlljet,Wlnu]
        self.Wlnu = Wlnu
        self.groups = [self.data,self.signal,self.LH,self.RH] + self.MC
        self.slt_ = slt_mu
        self.weight = None



class p8elstream(stream):
    def __init__(self):
        self.name = 'el'
        self.data = Egamma
        self.signal = P8Ztautau
        self.LH = LeftP8Ztautau
        self.RH = RightP8Ztautau
        self.ttbar = ttBar
        self.Zll = Zll
        self.MC = [ttBar,Zll,Wlnu]
        self.Wlnu = Wlnu
        self.groups = [self.data,self.signal,self.LH,self.RH] + self.MC
        self.slt_ = slt_el
        self.weight = None

class p8mustream(stream):
    def __init__(self):
        self.name = 'mu'
        self.data = Muons
        self.signal = P8Ztautau
        self.LH = LeftP8Ztautau
        self.RH = RightP8Ztautau
        self.ttbar = ttBar
        self.Zll = Zll
        self.MC = [ttBar,Zll,Wlnu]
        self.Wlnu = Wlnu
        self.groups = [self.data,self.signal,self.LH,self.RH] + self.MC
        self.slt_ = slt_mu
        self.weight = None

class Zmustream(stream):
    def __init__(self):
        self.name = 'mu'
        self.data = Muons
        self.signal = Ztautau_On
        self.LH = LeftZtautau_On
        self.RH = RightZtautau_On
        self.ttbar = ttBar
        self.Zll = Zll
        self.MC = [ttBar,Zll,Wlnu,Ztautau_Off]
        self.Zoff = Ztautau_Off
        self.Zlo = Ztautau_lo
        self.Zhi = Ztautau_hi
        self.Wlnu = Wlnu
        self.groups = [self.data,self.signal,self.LH,self.RH] + self.MC
        self.slt_ = slt_mu
        self.weight = None


class P8Zmustream(stream):
    def __init__(self):
        self.name = 'mu'
        self.data = Muons
        self.signal = P8Ztautau_On
        self.LH = LeftP8Ztautau_On
        self.RH = RightP8Ztautau_On
        self.ttbar = ttBar
        self.Zll = Zll
        self.MC = [ttBar,Zll,Wlnu,P8Ztautau_Off]
        self.Zoff = P8Ztautau_Off

        self.Wlnu = Wlnu
        self.groups = [self.data,self.signal,self.LH,self.RH] + self.MC
        self.slt_ = slt_mu
        self.weight = None



