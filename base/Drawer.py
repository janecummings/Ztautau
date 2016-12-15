import streams

class plot(dict):
    def __init__(self,name,variable,nbins,binlo,binhi):
        self.name = name
        self.variable = variable
        self.nbins = nbins
        self.binlo = binlo
        self.binhi = binhi
        self.groups = []
        self.rebin = 1
        self.xtitle = name
        self.xmax = None
        self.xmin = None
        self.leg = []
        self.blind = False
        self.units = ''
    def add_group(self,group,fill='',style=''):
        self.groups.append(group)
        self[group.name] = {}
        if fill:
            self[group.name]['fill'] = fill
        else: self[group.name]['fill'] = group.fill

class plot2d(dict):
    def __init__(self,x_name,x_variable,
                 y_name,y_variable,
                 nxbins,xbinlo,xbinhi,
                 nybins,ybinlo,ybinhi,
                 ):
        self.name = x_name + '_' + y_name
        self.xvariable = x_variable
        self.yvariable = y_variable
        self.nxbins = nxbins
        self.xbinlo = xbinlo
        self.xbinhi = xbinhi
        self.nybins = nybins
        self.ybinlo = ybinlo
        self.ybinhi = ybinhi
        self.groups = []
        self.rebin = 1.0
        self.xtitle = x_name
        self.ytitle = y_name
        self.rebinX = 2
        self.rebinY = 2

    def add_group(self,group,fill='',style=''):
        self.groups.append(group)
        self[group.name] = {}
        if fill:
            self[group.name]['fill'] = fill
        else: self[group.name]['fill'] = group.fill

def join1d(plotx,ploty):
    plot2 = plot2d(plotx.name,plotx.variable,
                   ploty.name,ploty.variable,
                   plotx.nbins,plotx.binlo,plotx.binhi,
                   ploty.nbins,ploty.binlo,ploty.binhi)
    plot2.xtitle = plotx.xtitle
    plot2.ytitle = ploty.xtitle 

    return plot2


########################################################################################################################
## ENTRIES #############################################################################################################
########################################################################################################################

### Tau Kinematic Variables ###
tau_et = plot('evtsel_tau_et','evtsel_tau_et',100,0,100)
tau_et.title = 'evtsel_tau_et'
tau_et.xtitle = 'p_{T}^{#tau} [GeV]'
tau_et.units = 'GeV'

tau_leadTrkPt = plot('tau_leadTrkPt','tau_leadTrkPt/1000.0',100,0,100)
tau_leadTrkPt.title = 'tau_leadTrkPt'
tau_leadTrkPt.xtitle = 'lead track p_{T} [GeV]'
tau_leadTrkPt.units = 'GeV'

tau_eta = plot('tau_eta','evtsel_tau_eta',120,-3,3)
tau_eta.xtitle = '#eta^{#tau}'

upsilon = plot('upsilon','2*tau_leadTrkPt/(1000.0 * evtsel_tau_et) - 1.0',200,-1.5,2.5)
upsilon.title = 'upsilon'
upsilon.xtitle = '#Upsilon'

upsilon_100 = plot('upilon100','2*tau_leadTrkPt/(1000.0 * evtsel_tau_et) - 1.0',100,-1.0,3)
upsilon_150 = plot('upsilon', '2*tau_leadTrkPt/(1000.0 * evtsel_tau_et) - 1.0',150,-1.0,2)
upsilon_true = plot('upsi','tau_true_upsilon',120,-1.2,1.2)
true_reco = join1d(upsilon_true,upsilon_150)

### Test Binning for Upsilon Fit ###
upsilon_40 = plot('upsilon','2*tau_leadTrkPt/(1000.0 * evtsel_tau_et) - 1.0',40,-1.0,1.5)
upsilon_24 = plot('upsilon','2*tau_leadTrkPt/(1000.0 * evtsel_tau_et) - 1.0',24,-1.0,2.0)
upsilon_28 = plot('upsilon','2*tau_leadTrkPt/(1000.0 * evtsel_tau_et) - 1.0',28,-1.5,2.0)
upsilon_20 =  plot('upsilon','2*tau_leadTrkPt/(1000.0 * evtsel_tau_et) - 1.0',20,-1.0,1.5)

tau_leadTrkPt_20 = plot('tau_leadTrkPt','tau_leadTrkPt/1000.0',40,0,80)
tau_leadTrkPt_20.xtitle = 'Lead Track p_{T} [GeV]'

upsilon_24.xtitle = '#Upsilon'
upsilon_20.xtitle = '#Upsilon'
upsilon_28.xtitle = '#Upsilon'
upsilon_40.xtitle = '#Upsilon'

# Fit development
upsilon_simp = plot('upsilon','2*tau_leadTrkPt/(1000.0 * evtsel_tau_et) - 1.0',10,-1.0,1.0)
# Single bin in control region for fitting scale factors
upsilon_CR = plot('upsilon_CR','2*tau_leadTrkPt/(1000.0 * evtsel_tau_et) - 1.0',1,-1.0,200001)


### Event Kinematic Variables ###
visMass = plot('visMass','evtsel_ditau_visibleMass',200,0,200)
visMass.title = 'visibleMass'
visMass.xtitle = 'm_{vis}(l,#tau)[GeV]'
visMass.units = 'GeV'

## Binning same as Z pT reweighting histogram
Zpt = plot('Zpt','evtsel_higgs_pt',32, -3, 93)
Zpt.xtitle = 'p_T^Z'
Zpt.units = 'GeV'

met = plot('met','evtsel_MET',80,0,80)
met.xtitle = 'E_{T}^{miss} [GeV]'
met.units = 'GeV'

transMass = plot('transMass','evtsel_transverseMass',200,0,200)
transMass.title = 'transverseMass'
transMass.xtitle = 'm_{T} [GeV]'
transMass.units = 'GeV'

deltaPhi = plot('deltaPhi','evtsel_dPhiSum',70,0,7)
deltaPhi.xtitle = '#Delta#phi(l,E_{T}^{miss}) + #Delta#phi(#tau_{had},E_{T}^{miss})'

lep_pt = plot('lep_pt','evtsel_lep_pt',80,0,80)
lep_pt.xtitle = 'p_{T}^{l} [GeV]'
lep_pt.units = 'GeV'

lep_eta = plot('lep_eta','evtsel_lep_eta',120,-3,3)
lep_eta.xtitle = '#eta^{l}'

cosDphi = plot('cosDphi','evtsel_ditau_cosDPhi',80,-2,2)
cosDphi.xtitle = 'cosDphi'

dRTauLep = plot('dRTauLep','evtsel_dRTauLep',80,0,8)
sumPt = plot('sumPt','evtsel_sumPt',300,0,600)

### True Event Kinematics 
mtt = plot('mtt', 'evtsel_true_mtt', 250,0,250)
mtt.title = 'Zmass'
mtt.xtitle = 'm_{#tau#tau} [GeV]'
mtt.units = 'GeV'

truZpt = plot('truZpt','true_Zpt', 50,0,150)
truZpt.xtitle = 'p_{T}^{Z}'
truZpt.units = 'GeV'


######################
# Profiles
######################

## Polarization Observable Correlations
ups_et = join1d(upsilon,tau_et)
et_ups = join1d(tau_et,upsilon)
ups_ldpt = join1d(upsilon_40,tau_leadTrkPt_20)
ldpt_ups = join1d(tau_leadTrkPt,upsilon)
et_ldpt = join1d(tau_et,tau_leadTrkPt)
ldpt_et = join1d(tau_leadTrkPt,tau_et)

ups_met = join1d(upsilon_40, met)

#######################
# W Correction Studies 
##########################
pTratio = plot('pTratio','evtsel_pTratio',40,0,10)
pTratio.xtitle = 'p_{T}(#tau_{had})/p_{T}(l)'

lep_dEta = plot('lep_dEta','abs(evtsel_lep_eta - evtsel_tau_eta)',40,0,10)
lep_dEta.xtitle = '#Delta#eta(#tau_{had},l)'

#########################
# Tau Identification 
##########################
BDTJetScore = plot('BDTJetScore','evtsel_tau_JetBDTScore',80,0.2,1.0)
BDTJetScore.xtitle = 'BDT Jet Score'
BDTJetScore.leg = [0.18,0.55,0.45,0.9]

BDTEleScore = plot('BDTEleScore','tau_BDTEleScore',20,0,1.) 

emfrac = plot('emfrac','evtsel_tau_EMFrac',40,0,2)
emfrac.xtitle = 'f_{EM}'

jet_emfrac = plot('jet_emfrac','evtsel_tau_jet_emfrac',40,0,2)

ftrk = plot('ftrk','evtsel_tau_FTrk',30,0,3)
ftrk.xtitle = 'f_{trk}'

centfrac = plot('centfrac','evtsel_tau_CentFrac',30,0,1.5)
centfrac.xtitle = 'f_{cent}'
centfrac.leg = [0.18,0.55,0.45,0.9]

## Tau ID correlations 
bdt_upsilon = join1d(BDTJetScore,upsilon)
emfrac_upsilon = join1d(emfrac,upsilon)
centfrac_upsilon = join1d(centfrac,upsilon)
leadTrkPt_upsilon = join1d(tau_leadTrkPt,upsilon)

#####################
#Tracking variables
#####################

nBLHits = plot('nBLHits','tau_track_nBLHits',6,0,6)
nPixHits = plot('nPixHits','tau_track_nPixHits',6,0,6)


ldTrkEta = plot("ldTrkEta","tau_leadTrack_eta",120,-3,3)
ldTrkPhi = plot("ldTrkPhi","tau_leadTrack_phi",64,-3.2,3.2)
tau_phi = plot("tau_phi","evtsel_tau_phi",64,-3.2,3.2)
cluster_eta = plot("cluster_eta","evtsel_tau_cluster_eta",120,-3,3)
cluster_phi = plot("cluster_phi","evtsel_tau_cluster_phi",64,-3.2,3.2)
dEta = plot("dEta","evtsel_tau_cluster_eta - tau_leadTrack_eta",100,-0.5,0.5)
dPhi = plot("dPhi","evtsel_tau_cluster_phi - tau_leadTrack_phi",100,-0.5,0.5)

avg_int = plot("avg_int",'evtsel_averageIntPerXing',40,0,40)
avg_int.xtitle = 'evtsel_averageIntPerXing'
n_vert = plot('n_vert','evtsel_vertices',35,0,35)
n_vert.xtitle = 'evtsel_vertices'

## Angular Correlations
dEta_dPhi = join1d(dEta,dPhi)
ldPt_dEta = join1d(tau_leadTrkPt,dEta)
ldPt_dPhi = join1d(tau_leadTrkPt,dPhi)


###################
# Event weights 
###################
evtsel_weight = plot('evtsel_weight','evtsel_weight',100,0,10)
evtsel_weight_PU = plot('evtsel_weight_PU','evtsel_weight_PU',100,0,10)
evtsel_weight_lep = plot('evtsel_weight_lep','evtsel_weight_lep',100,0,10)
evtsel_weight_lep_reco = plot('evtsel_weight_lep_reco','evtsel_weight_lep_reco',100,0,10) 
evtsel_weight_lep_id = plot('evtsel_weight_lep_id','evtsel_weight_lep_id',100,0,10)   
evtsel_weight_lep_trigger = plot('evtsel_weight_lep_trigger','evtsel_weight_lep_trigger',100,0,10)  
evtsel_weight_lep_iso = plot('evtsel_weight_lep_iso','evtsel_weight_lep_iso',100,0,10)  
evtsel_weight_SLT = plot('evtsel_weight_SLT','evtsel_weight_SLT',100,0,10)      
evtsel_weight_tau_el_overlap = plot('evtsel_weight_tau_el_overlap','evtsel_weight_tau_el_overlap',100,0,10)    
evtsel_weight_tau_id = plot('evtsel_weight_tau_id','evtsel_weight_tau_id',100,0,10)   


###################
# Truth variables 
###################
trueTau_upsilon = plot('evtsel_trueTau_upsilon', 'evtsel_trueTau_upsilon', 40, -1, 1)
trueTau_upsilon.xtitle = '#Upsilon'

trueTau_fEM = plot('evtsel_trueTau_fEM', 'evtsel_trueTau_fEM', 40, 0, 2)

trueTau_decayMode = plot('evtsel_trueTau_decayMode','evtsel_trueTau_decayMode', 17, 0, 17)
trueTau_decayMode.xtitle = 'Decay Mode'

truethad_pt = plot('evtsel_truethad_pt','evtsel_truethad_pt', 80, 0, 80)
truethad_pt.xtitle = 'p_{T}^{#tau}'

truethad_eta = plot('evtsel_truethad_eta','evtsel_truethad_eta', 40, -4, 4)
truethad_eta.xtitle = '#eta^{#tau}'

truethad_phi = plot('evtsel_truethad_phi', 'evtsel_truethad_phi', 32, -3.2, 3.2) 
truethad_phi.xtitle = '#phi^{#tau}'

truethad_EMfrac = plot('evtsel_truethad_EMfrac', 'evtsel_truethad_EMfrac', 100, 0, 1)

trueZ_m = plot('evtsel_trueZ_m', 'evtsel_trueZ_m', 100, 60, 160)
trueZ_pt = plot('evtsel_trueZ_pt', 'evtsel_trueZ_pt', 100, 0, 100)

true_mtt = plot('evtsel_true_mtt', 'evtsel_true_mtt', 100, 60, 160)
true_mtt.xtitle = 'm_{#tau#tau}'



