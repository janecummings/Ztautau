from Sample import Cut,CutSet
from stream import *
from Groups import *


################################################################################################################
# STANDARD SELECTION   #########################################################################################
################################################################################################################
all_ = Cut('',name='all')
# Triggers
slt_ = Cut('') 
slt_el = Cut('evtsel_is_el',name='evtsel_is_el')
slt_mu = Cut('evtsel_is_mu',name='evtsel_is_mu')
# bch cleaning
bch_ = Cut('!evtsel_is_badTightBch', name='bch')
# Lepton Isolation
isoLep_ = Cut('evtsel_is_isoLep',name='isoLep')
# Tau Identification
is_tau_ = Cut('evtsel_is_tau',name='is_tau')
# Dilepton Veto
dilepVeto_ = Cut('evtsel_is_dilepVeto',name='dilepVeto')
# Lepton vetos
lepveto_ = Cut('evtsel_is_conf_lep_veto_medium',name = 'lepVeto') 
# Transverse mass < 30 GeV
transMass30_ = Cut('evtsel_transverseMass', '<', 30,name='transMass30')
# Sum DeltaPhi 
dPhiSum_ = Cut('evtsel_dPhiSum','<',3.5,name='dPhi')
# Single Prong Taus
sp_ = Cut('evtsel_tau_numTrack','==',1,name='SP')
# Visible mass 
vismass_lo_ = Cut('evtsel_ditau_visibleMass','>=',40)
vismass_hi_ = Cut('evtsel_ditau_visibleMass','<=',85)
vismass_ = CutSet([vismass_lo_,vismass_hi_], name = 'visMass')
# Upsilon < 1.5
ups_ = Cut('(2*tau_leadTrkPt/(1000.0 * evtsel_tau_et) - 1.0)','<',1.5,name="upscut") 
# Opposite Sign and Same sign Requirements
OS_ = Cut('evtsel_is_oppositeSign',name='OS')
SS_ = Cut('!evtsel_is_oppositeSign',name = 'SS')


################################################################################################################
# Control regions   ############################################################################################
################################################################################################################
# Z Control Region
dilepCR_ = Cut('evtsel_is_dilepCR')
# W Control Region 
antiTM70_ = Cut('evtsel_transverseMass','>',70)
antidPhiSum_ = Cut('evtsel_dPhiSum','>=',3.5)
# WT Control Region 
leadJet_ = Cut('evtsel_jet_leading_pt','>','30')
notMV1_ = Cut('!evtsel_MV1')
# Top CR 
MV1_ = Cut('evtsel_MV1')
antiTM50_ = Cut('evtsel_transverseMass','>',50)
# QCD Control Region
notIso_ = Cut('!evtsel_is_isoLep')


################################################################################################################
# Studies           ############################################################################################
################################################################################################################

### Standard selection alternatives
tight_ = Cut('evtsel_tau_is_Tight',name='Tight')
looseVeto_ = Cut('evtsel_is_conf_lep_veto')
muVeto_ = Cut('!evtsel_tau_is_muVeto',name = 'muVeto')
transMass40_ = Cut('evtsel_transverseMass', '<', 40,name='transMass40')
# Region Definitions
tight_signal = [all_,slt_, isoLep_,tight_,dilepVeto_,looseVeto_,sp_,transMass30_,dPhiSum_]
tight_wcr = [slt_, tight_,isoLep_,dilepVeto_,antiTM70_,antidPhiSum_,sp_,looseVeto_]
qcd = [slt_,is_tau_,notIso_,transMass30_,dPhiSum_, dilepVeto_,sp_,lepveto_]

### Real taus in ttbar
met_ = Cut('evtsel_MET','>',30)
ht_ = Cut('evtsel_sumPt','>',170)
jetpt_ = CutSet([Cut('evtsel_jet_leading_pt','>',25), Cut('evtsel_jet_subleading_pt','>',25)])
jeteta_ = CutSet([Cut('abs(evtsel_jet_leading_eta)','<',2.5), Cut('abs(evtsel_jet_subleading_eta)','<',2.5)])
taupt_ = Cut('evtsel_tau_et','<',100,name='tauEt100')

### Cut on lead track pt
leadPt_ = Cut('tau_leadTrkPt','>',4000.0,name='leadPt4Gev')

### Cutting on number of pi0s
nopi0_ = Cut('evtsel_tau_pi0_n','==',0,name='NoPi0')
pi0_ = Cut('evtsel_tau_pi0_n','==',1,name='Pi0')
Npi0_ = Cut('evtsel_tau_pi0_n','>',0,name='NoPi0')
Api0_ = Cut('evtsel_tau_pi0_n','==',1,name='Pi0')

### Modeling in regions of pile-up
mu0_ = Cut('evtsel_averageIntPerXing','<',10,name='mu0')
mu10_ = CutSet([Cut('evtsel_averageIntPerXing','>=',10),
                Cut('evtsel_averageIntPerXing','<',20)],name='mu10')
mu20_ = CutSet([Cut('evtsel_averageIntPerXing','>=',20),
                Cut('evtsel_averageIntPerXing','<',30)],name='mu20')
mu30_ = Cut('evtsel_averageIntPerXing','>=',30,name='mu30')
mu25_ = Cut('evtsel_averageIntPerXing','<',25)

################################################################################################################
# Base Regions      ############################################################################################
################################################################################################################

signal = [all_,bch_,slt_,isoLep_,is_tau_,dilepVeto_,lepveto_,sp_,transMass30_,dPhiSum_, vismass_, ups_]
signal_thesis = [all_, dilepVeto_, slt_, isoLep_, is_tau_, lepveto_, sp_, OS_, transMass30_, dPhiSum_, ups_ ]
wcr = [slt_,bch_, is_tau_,isoLep_,dilepVeto_,antiTM70_,antidPhiSum_,sp_,lepveto_,vismass_,ups_]
qcd = [slt_,bch_,is_tau_,notIso_,transMass30_,dPhiSum_, dilepVeto_,sp_,lepveto_,vismass_,ups_]

signal_val = [all_,OS_,slt_, isoLep_,is_tau_,dilepVeto_,lepveto_,sp_,transMass30_,dPhiSum_]
zlcr = [slt_,is_tau_,isoLep_,transMass30_,dilepCR_,sp_,lepveto_]
wtcr = [slt_, is_tau_,isoLep_,dilepVeto_,antiTM70_,leadJet_,notMV1_,sp_,lepveto_]
tcr = [slt_, is_tau_,isoLep_,dilepVeto_,leadJet_,MV1_,antiTM50_,sp_,lepveto_]
wcr_noTM = [slt_, is_tau_,isoLep_,dilepVeto_,antidPhiSum_,sp_,lepveto_]


REGIONS = {'SR': signal, 'ZLCR': zlcr, 'WTCR': wtcr, 'TCR': tcr, 'WCR': wcr, 'QCD': qcd}
def make_all_regions(stream,cut=[]):
    for region in REGIONS:
        cuts = REGIONS[region]+cut
        stream.make_region(cuts=cuts,name=region)
def make_regions(stream,cut=[],regions=[]):
	for region in REGIONS:
		if region in regions:
			cuts = REGIONS[region]+cut
			stream.make_region(cuts=cuts,name=region)

################################################################################################################
# Selections        ############################################################################################
################################################################################################################
############### Standard ############################################
el = elstream(); mu = mustream()
mutheory = mutheory();eltheory = eltheory()
el_modes = el_modes(); mu_modes = mu_modes()

# fill all regions 
for st in [el,mu, el_modes, mu_modes, mutheory, eltheory]:
    make_all_regions(st)
    st.make_selection(cuts = signal + [OS_])
    st.make_selection(cuts = signal[0:2] + [OS_] + signal[2:], name = 'noteflow')
    st.make_selection(cuts = signal_val, name = 'val' )
    st.make_selection(cuts = signal_thesis, name = 'cutflow')

################# Studies ##########################################
# Tight Tau ID
mu_tight = mustream()
el_tight = elstream()
for st in [mu_tight,el_tight]:
    st.make_region(cuts = tight_signal, name = 'SR')
    st.make_region(cuts = tight_wcr, name = 'WCR')
    st.weight = 'evtsel_weight_tight'

# Real taus in ttbar 
mu_ttbar = mustream()
mu_ttbar.LH = None; mu_ttbar.RH = None
mu_ttbar.signal = ttBar_real
mu_ttbar.MC = [Zll,Wlnu,Ztautau,ttBar_fake]

ttbar_cuts =[slt_,bch_,is_tau_,isoLep_,dilepVeto_,lepveto_,MV1_,met_,ht_,jetpt_,jeteta_,taupt_,sp_,ups_]
ttbar_qcd_cuts = [slt_,bch_,is_tau_,notIso_,dilepVeto_,lepveto_,MV1_,met_,ht_,jetpt_,jeteta_,taupt_,sp_, ups_]
mu_ttbar.make_region(cuts=ttbar_cuts,name='TTBAR')
mu_ttbar.make_region(cuts=ttbar_qcd_cuts,name='QCD')
mu_ttbar.make_region(cuts=wcr,name='WCR')


el_ttbar = elstream()
el_ttbar.LH = None; el_ttbar.RH = None
el_ttbar.signal = ttBar_real
el_ttbar.MC = [Zll,Wlnu,Ztautau,ttBar_fake]

el_ttbar.make_region(cuts=ttbar_cuts,name='TTBAR')
el_ttbar.make_region(cuts=ttbar_qcd_cuts,name='QCD')
el_ttbar.make_region(cuts=wcr,name='WCR')


