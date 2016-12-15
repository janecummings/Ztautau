from Sample import Sample,TruthSample,Group,Cut
import ROOT
###############################################################################################################################################
########## SAMPLES ############################################################################################################################
###############################################################################################################################################

# Signal Aplgen + Pythia6 Z->tautau w. Perugia2011C tune (CTEQ6L1 PDF)
AlpPythZtautauNp0 = Sample(name="ZtautauNp0", DSID = 147121, version = 'l4', xsec = 718.87*1.18*1.0, xstr = "718.87*1.18*1.0")
AlpPythZtautauNp1 = Sample(name="ZtautauNp1", DSID = 147122, version = '', xsec = 175.75*1.18*1.0, xstr = "175.75*1.18*1.0")
AlpPythZtautauNp2 = Sample(name="ZtautauNp2", DSID = 147123, version = '', xsec = 58.856*1.18*1.0, xstr = "58.856*1.18*1.0")
AlpPythZtautauNp3 = Sample(name="ZtautauNp3", DSID = 147124, version = '', xsec = 15.667*1.18*1.0, xstr = "15.667*1.18*1.0")
AlpPythZtautauNp4 = Sample(name="ZtautauNp4", DSID = 147125, version = '', xsec = 4.0121*1.18*1.0, xstr = "4.0121*1.18*1.0")
AlpPythZtautauNp5incl = Sample(name="ZtautauNp5incl", DSID = 147126, version = '', xsec = 1.256*1.18*1.0, xstr = "1.256*1.18*1.0")

# Tau Spinner Polarization in Alpgen Np Samples
p0 = -0.144452519763
p1 = -0.138878336093
p2 = -0.120975830456
p3 = -0.109910391528
p4 = -0.106513371532
p5 = -0.099161344576
# Averaged polarization
Alpol = -0.1413

# Left-Handed Signal
LeftZtautauNp0 = Sample(name="LeftZtautauNp0", DSID = 147121, version = 'LEFT', xsec = 718.87*1.18*1.0, xpol = 0.5*(1-p0))
LeftZtautauNp1 = Sample(name="LeftZtautauNp1", DSID = 147122, version = 'LEFT', xsec = 175.75*1.18*1.0, xpol = 0.5*(1-p1))
LeftZtautauNp2 = Sample(name="LeftZtautauNp2", DSID = 147123, version = 'LEFT', xsec = 58.856*1.18*1.0, xpol = 0.5*(1-p2))
LeftZtautauNp3 = Sample(name="LeftZtautauNp3", DSID = 147124, version = 'LEFT', xsec = 15.667*1.18*1.0, xpol = 0.5*(1-p3))
LeftZtautauNp4 = Sample(name="LeftZtautauNp4", DSID = 147125, version = 'LEFT', xsec = 4.0121*1.18*1.0, xpol = 0.5*(1-p4))
LeftZtautauNp5incl = Sample(name="LeftZtautauNp5incl", DSID = 147126, version = 'LEFT', xsec = 1.256*1.18*1.0, xpol = 0.5*(1-p5))

# Right-Handed Signal
RightZtautauNp0 = Sample(name="RightZtautauNp0", DSID = 147121, version = 'RIGHT', xsec = 718.87*1.18, xpol =  0.5*(1+p0))
RightZtautauNp1 = Sample(name="RightZtautauNp1", DSID = 147122, version = 'RIGHT', xsec = 175.75*1.18, xpol =  0.5*(1+p1))
RightZtautauNp2 = Sample(name="RightZtautauNp2", DSID = 147123, version = 'RIGHT', xsec = 58.856*1.18, xpol =  0.5*(1+p2))
RightZtautauNp3 = Sample(name="RightZtautauNp3", DSID = 147124, version = 'RIGHT', xsec = 15.667*1.18, xpol =  0.5*(1+p3))
RightZtautauNp4 = Sample(name="RightZtautauNp4", DSID = 147125, version = 'RIGHT', xsec = 4.0121*1.18, xpol =  0.5*(1+p4))
RightZtautauNp5incl = Sample(name="RightZtautauNp5incl", DSID = 147126, version = 'RIGHT', xsec = 1.256*1.18, xpol =  0.5*(1+p5))

# Pythia8 Z->tautau production with the AU2 CTEQ6L1 tune
LeftP8 = Sample(name="LeftZtautauP8", DSID = 147818, version = 'LEFT', xsec = 974.42*1.18*1)
RightP8 = Sample(name="LeftZtautauP8", DSID = 147818, version = 'RIGHT', xsec = 974.42*1.18*1)
ZtautauP8 = Sample(name="ZtautauP8", DSID = 147818, xsec = 974.42*1.18*1)

# Pythia + Pythia6 ttbar
Pythia_ttbar = Sample(name="Pythia_ttbar", DSID = 117050, xsec = 252.89*1*0.543, xstr = "252.89*1*0.543")

# Alpgen + Pythia6 W + jets
AlpgenPythiaWenuNp0 = Sample(name="AlpgenPythiaWenuNp0", DSID = 147025, xsec = 8127.3*1.15*1, xstr = "8127.3*1.15*1")#8136.8*1.15*1)
AlpgenPythiaWenuNp1 = Sample(name="AlpgenPythiaWenuNp1", DSID = 147026, xsec = 1792.9*1.15*1, xstr = "1792.9*1.15*1")#1791.5*1.15*1)
AlpgenPythiaWenuNp2 = Sample(name="AlpgenPythiaWenuNp2", DSID = 147027, xsec = 542.3*1.15*1, xstr = "542.3*1.15*1")#541.6*1.15*1)
AlpgenPythiaWenuNp3 = Sample(name="AlpgenPythiaWenuNp3", DSID = 147028, xsec = 147.65*1.15*1, xstr = "147.65*1.15*1")#146.65*1.15*1)
AlpgenPythiaWenuNp4 = Sample(name="AlpgenPythiaWenuNp4", DSID = 147029, xsec = 37.751*1.15*1, xstr = "37.751*1.15*1")#37.295*1.15*1)
AlpgenPythiaWenuNp5incl = Sample(name="AlpgenPythiaWenuNp5incl", DSID = 147030, xsec = 11.369*1.15*1, xstr = "11.369*1.15*1")#11.916*1.15*1)

AlpgenPythiaWmunuNp0 = Sample(name="AlpgenPythiaWmunuNp0", DSID = 147033, xsec = 8127.1*1.15*1, xstr = "8127.1*1.15*1")#8133.4*1.15*1)
AlpgenPythiaWmunuNp1 = Sample(name="AlpgenPythiaWmunuNp1", DSID = 147034, xsec = 1792.8*1.15*1, xstr = "1792.8*1.15*1")#1792.7*1.15*1)
AlpgenPythiaWmunuNp2 = Sample(name="AlpgenPythiaWmunuNp2", DSID = 147035, xsec = 542.42*1.15*1, xstr = "542.42*1.15*1")#541.27*1.15*1)
AlpgenPythiaWmunuNp3 = Sample(name="AlpgenPythiaWmunuNp3", DSID = 147036, xsec = 147.66*1.15*1, xstr = "147.66*1.15*1")#146.49*1.15*1)
AlpgenPythiaWmunuNp4 = Sample(name="AlpgenPythiaWmunuNp4", DSID = 147037, xsec = 37.76*1.15*1, xstr = "37.76*1.15*1")#37.334*1.15*1)
AlpgenPythiaWmunuNp5incl = Sample(name="AlpgenPythiaWmunuNp5incl", DSID = 147038, xsec = 11.934*1.15*1, xstr = "11.934*1.15*1")#11.414*1.15*1)

AlpgenPythiaWtaunuNp0 = Sample(name="AlpgenPythiaWtaunuNp0", DSID = 147041, xsec = 8127.1*1.15*1, xstr = "8127.1*1.15*1")#8135.7*1.15*1)
AlpgenPythiaWtaunuNp1 = Sample(name="AlpgenPythiaWtaunuNp1", DSID = 147042, xsec = 1792.5*1.15*1, xstr = "1792.5*1.15*1")#1793.7*1.15*1)
AlpgenPythiaWtaunuNp2 = Sample(name="AlpgenPythiaWtaunuNp2", DSID = 147043, xsec = 542.21*1.15*1, xstr = "542.21*1.15*1")#541.24*1.15*1)
AlpgenPythiaWtaunuNp3 = Sample(name="AlpgenPythiaWtaunuNp3", DSID = 147044, xsec = 147.64*1.15*1, xstr = "147.64*1.15*1")#146.48*1.15*1)
AlpgenPythiaWtaunuNp4 = Sample(name="AlpgenPythiaWtaunuNp4", DSID = 147045, xsec = 37.738*1.15*1, xstr = "37.738*1.15*1")#37.264*1.15*1)
AlpgenPythiaWtaunuNp5incl = Sample(name="AlpgenPythiaWtaunuNp5incl", DSID = 147046, xsec = 11.905*1.15*1, xstr = "11.905*1.15*1")#11.537*1.15*1)

# Alpgen + AlpgenPythia6 Z + jets 
AlpgenPythiaZeeNp0 = Sample(name="AlpgenPythiaZeeNp0", DSID = 147105, xsec = 718.97*1.18*1, xstr = "718.97*1.18*1")
AlpgenPythiaZeeNp1 = Sample(name="AlpgenPythiaZeeNp1", DSID = 147106, xsec = 175.7*1.18*1, xstr = "175.7*1.18*1")
AlpgenPythiaZeeNp2 = Sample(name="AlpgenPythiaZeeNp2", DSID = 147107, xsec = 58.875*1.18*1, xstr = "58.875*1.18*1")
AlpgenPythiaZeeNp3 = Sample(name="AlpgenPythiaZeeNp3", DSID = 147108, xsec = 15.636*1.18*1, xstr = "15.636*1.18*1")
AlpgenPythiaZeeNp4 = Sample(name="AlpgenPythiaZeeNp4", DSID = 147109, xsec = 4.0116*1.18*1, xstr = "4.0116*1.18*1")
AlpgenPythiaZeeNp5incl = Sample(name="AlpgenPythiaZeeNp5incl", DSID = 147110, xsec = 1.2592*1.18*1, xstr = "1.2592*1.18*1")

AlpgenPythiaZmumuNp0 = Sample(name="AlpgenPythiaZmumuNp0", DSID = 147113, xsec = 719.16*1.18*1, xstr = "719.16*1.18*1")
AlpgenPythiaZmumuNp1 = Sample(name="AlpgenPythiaZmumuNp1", DSID = 147114, xsec = 175.74*1.18*1, xstr = "175.74*1.18*1")
AlpgenPythiaZmumuNp2 = Sample(name="AlpgenPythiaZmumuNp2", DSID = 147115, xsec = 58.882*1.18*1, xstr = "58.882*1.18*1")
AlpgenPythiaZmumuNp3 = Sample(name="AlpgenPythiaZmumuNp3", DSID = 147116, xsec = 15.673*1.18*1, xstr = "15.673*1.18*1")
AlpgenPythiaZmumuNp4 = Sample(name="AlpgenPythiaZmumuNp4", DSID = 147117, xsec = 4.0057*1.18*1, xstr = "4.0057*1.18*1")
AlpgenPythiaZmumuNp5incl = Sample(name="AlpgenPythiaZmumuNp5incl", DSID = 147118, xsec = 1.2544*1.18*1, xstr = "1.2544*1.18*1")

####################################################################################################
############## DATA SAMPLES ########################################################################
####################################################################################################

# Muon
Muons_PeriodA = Sample(name="Muons_PeriodA", DSID = "periodA", stream = "Muons")
Muons_PeriodB = Sample(name="Muons_PeriodB", DSID = "periodB", stream = "Muons")
Muons_PeriodC = Sample(name="Muons_PeriodC", DSID = "periodC", stream = "Muons")
Muons_PeriodD = Sample(name="Muons_PeriodD", DSID = "periodD", stream = "Muons")
Muons_PeriodE = Sample(name="Muons_PeriodE", DSID = "periodE", stream = "Muons")
Muons_PeriodG = Sample(name="Muons_PeriodG", DSID = "periodG", stream = "Muons")
Muons_PeriodH = Sample(name="Muons_PeriodH", DSID = "periodH", stream = "Muons")
Muons_PeriodI = Sample(name="Muons_PeriodI", DSID = "periodI", stream = "Muons")
Muons_PeriodJ = Sample(name="Muons_PeriodJ", DSID = "periodJ", stream = "Muons")
Muons_PeriodL = Sample(name="Muons_PeriodL", DSID = "periodL", stream = "Muons")

# Electron 
Egamma_PeriodA = Sample(name="Egamma_PeriodA", DSID = "periodA", stream = "Egamma")
Egamma_PeriodB = Sample(name="Egamma_PeriodB", DSID = "periodB", stream = "Egamma")
Egamma_PeriodC = Sample(name="Egamma_PeriodC", DSID = "periodC", stream = "Egamma")
Egamma_PeriodD = Sample(name="Egamma_PeriodD", DSID = "periodD", stream = "Egamma")
Egamma_PeriodE = Sample(name="Egamma_PeriodE", DSID = "periodE", stream = "Egamma")
Egamma_PeriodG = Sample(name="Egamma_PeriodG", DSID = "periodG", stream = "Egamma")
Egamma_PeriodH = Sample(name="Egamma_PeriodH", DSID = "periodH", stream = "Egamma")
Egamma_PeriodI = Sample(name="Egamma_PeriodI", DSID = "periodI", stream = "Egamma")
Egamma_PeriodJ = Sample(name="Egamma_PeriodJ", DSID = "periodJ", stream = "Egamma")
Egamma_PeriodL = Sample(name="Egamma_PeriodL", DSID = "periodL", stream = "Egamma")


####################################################################################################
############## OTHER SAMPLES #######################################################################
####################################################################################################

# POWHEG+Pythia8 Z->tautau production without lepton filter and AU2 CT10 tune
Pythia_Ztautau = Sample(name="Pythia_Ztautau", DSID = 147808, xsec = 1109.9*1*1)

# McAtNlo+fHerwig/Jimmy ttbar production with TTbarWToLeptonFilter - CT10 PDF, AUET2 tune
ttbar = Sample(name="ttbar", DSID = 105200, xsec = 252.89*1*0.543)
ttbar_allhad = Sample(name="ttbar_allhad", DSID = 105204, xsec = 252.89*1*0.457)

# ALPGEN+Jimmy Z(->ll) process with AUET2 tune and pdf CTEQ6L1
JimmyZmumuNp1 = Sample(name="JimmyZmumuNp1", DSID = 107661, xsec = 154.77*1.23*1)
JimmyZmumuNp3 = Sample(name="JimmyZmumuNp3", DSID = 107663, xsec = 48.912*1.23*1)
JimmyZmumuNp4 = Sample(name="JimmyZmumuNp4", DSID = 107664, xsec = 3.7838*1.23*1)
JimmyZmumuNp5 = Sample(name="JimmyZmumuNp5", DSID = 107665, xsec = 1.1148*1.23*1)

# ALPGEN+Jimmy W(->Lnu) process with AUET2 tune and pdf CTEQ6L1
JimmyWenuNp0 = Sample(name="JimmyWenuNp0", DSID = 107680, xsec = 8037.1*1.19*1)
JimmyWenuNp1 = Sample(name="JimmyWenuNp1", DSID = 107681, xsec = 1579.2*1.19*1)
JimmyWenuNp2 = Sample(name="JimmyWenuNp2", DSID = 107682, xsec = 477.2*1.19*1)
JimmyWenuNp3 = Sample(name="JimmyWenuNp3", DSID = 107683, xsec = 133.93*1.19*1)
JimmyWenuNp4 = Sample(name="JimmyWenuNp4", DSID = 107684, xsec = 35.622*1.19*1)
JimmyWenuNp5 = Sample(name="JimmyWenuNp5", DSID = 107685, xsec = 10.553*1.19*1)

JimmyWmunuNp0 = Sample(name="JimmyWmunuNp0", DSID = 107690, xsec = 8040*1.19*1)
JimmyWmunuNp1 = Sample(name="JimmyWmunuNp1", DSID = 107691, xsec = 1580.3*1.19*1)
JimmyWmunuNp2 = Sample(name="JimmyWmunuNp2", DSID = 107692, xsec = 477.5*1.19*1)
JimmyWmunuNp3 = Sample(name="JimmyWmunuNp3", DSID = 107693, xsec = 133.94*1.19*1)
JimmyWmunuNp4 = Sample(name="JimmyWmunuNp4", DSID = 107694, xsec = 35.636*1.19*1)

JimmyWtaunuNp0 = Sample(name="JimmyWtaunuNp0", DSID = 107700, xsec = 8035.8*1.19*1)
JimmyWtaunuNp1 = Sample(name="JimmyWtaunuNp1", DSID = 107701, xsec = 1579.8*1.19*1)
JimmyWtaunuNp2 = Sample(name="JimmyWtaunuNp2", DSID = 107702, xsec = 477.55*1.19*1)
JimmyWtaunuNp3 = Sample(name="JimmyWtaunuNp3", DSID = 107703, xsec = 133.79*1.19*0.12038)
JimmyWtaunuNp4 = Sample(name="JimmyWtaunuNp4", DSID = 107704, xsec = 35.583*1.19*0.2512)
JimmyWtaunuNp5 = Sample(name="JimmyWtaunuNp5", DSID = 107705, xsec = 10.54*1.19*0.44231)

# JetTauEtMiss Stream
JetTauEtmiss_PeriodA = Sample(name="JetTauEtmiss_PeriodA", DSID = "periodA", stream = "JetTauEtmiss")
JetTauEtmiss_PeriodB = Sample(name="JetTauEtmiss_PeriodB", DSID = "periodB", stream = "JetTauEtmiss")
JetTauEtmiss_PeriodC = Sample(name="JetTauEtmiss_PeriodC", DSID = "periodC", stream = "JetTauEtmiss")
JetTauEtmiss_PeriodD = Sample(name="JetTauEtmiss_PeriodD", DSID = "periodD", stream = "JetTauEtmiss")
JetTauEtmiss_PeriodE = Sample(name="JetTauEtmiss_PeriodE", DSID = "periodE", stream = "JetTauEtmiss")
JetTauEtmiss_PeriodG = Sample(name="JetTauEtmiss_PeriodG", DSID = "periodG", stream = "JetTauEtmiss")
JetTauEtmiss_PeriodH = Sample(name="JetTauEtmiss_PeriodH", DSID = "periodH", stream = "JetTauEtmiss")
JetTauEtmiss_PeriodI = Sample(name="JetTauEtmiss_PeriodI", DSID = "periodI", stream = "JetTauEtmiss")
JetTauEtmiss_PeriodJ = Sample(name="JetTauEtmiss_PeriodJ", DSID = "periodJ", stream = "JetTauEtmiss")
JetTauEtmiss_PeriodL = Sample(name="JetTauEtmiss_PeriodL", DSID = "periodL", stream = "JetTauEtmiss")


####################################################################################################
############## TRUTH SAMPLES #######################################################################
####################################################################################################
PATH_TO_TRUE = '/data/cummings/common/truth'
TrueNp0 = TruthSample(name='Np0')
TrueNp1 = TruthSample(name='Np1')
TrueNp2 = TruthSample(name='Np2')
TrueNp3 = TruthSample(name='Np3')
TrueNp4 = TruthSample(name='Np4')
TrueNp5 = TruthSample(name='Np5')

# Nominal Alpgen Files
LeftAlpgenNp0 = Sample(name="LeftAlpgenNp0", DSID = 147121, version = 'LEFT', xsec = 718.87*1.18*1.0, xpol = 0.5*(1-p0), theory = True)
LeftAlpgenNp1 = Sample(name="LeftAlpgenNp1", DSID = 147122, version = 'LEFT', xsec = 175.75*1.18*1.0, xpol = 0.5*(1-p1), theory = True)
LeftAlpgenNp2 = Sample(name="LeftAlpgenNp2", DSID = 147123, version = 'LEFT', xsec = 58.856*1.18*1.0, xpol = 0.5*(1-p2), theory = True)
LeftAlpgenNp3 = Sample(name="LeftAlpgenNp3", DSID = 147124, version = 'LEFT', xsec = 15.667*1.18*1.0, xpol = 0.5*(1-p3), theory = True)
LeftAlpgenNp4 = Sample(name="LeftAlpgenNp4", DSID = 147125, version = 'LEFT', xsec = 4.0121*1.18*1.0, xpol = 0.5*(1-p4), theory = True)
LeftAlpgenNp5incl = Sample(name="LeftAlpgenNp5incl", DSID = 147126, version = 'LEFT', xsec = 1.256*1.18*1.0, xpol = 0.5*(1-p5), theory = True)

RightAlpgenNp0 = Sample(name="RightAlpgenNp0", DSID = 147121, version = 'RIGHT', xsec = 718.87*1.18, xpol =  0.5*(1+p0), theory = True)
RightAlpgenNp1 = Sample(name="RightAlpgenNp1", DSID = 147122, version = 'RIGHT', xsec = 175.75*1.18, xpol =  0.5*(1+p1), theory = True)
RightAlpgenNp2 = Sample(name="RightAlpgenNp2", DSID = 147123, version = 'RIGHT', xsec = 58.856*1.18, xpol =  0.5*(1+p2), theory = True)
RightAlpgenNp3 = Sample(name="RightAlpgenNp3", DSID = 147124, version = 'RIGHT', xsec = 15.667*1.18, xpol =  0.5*(1+p3), theory = True)
RightAlpgenNp4 = Sample(name="RightAlpgenNp4", DSID = 147125, version = 'RIGHT', xsec = 4.0121*1.18, xpol =  0.5*(1+p4), theory = True)
RightAlpgenNp5incl = Sample(name="RightAlpgenNp5incl", DSID = 147126, version = 'RIGHT', xsec = 1.256*1.18, xpol =  0.5*(1+p5), theory = True)

LeftAlpgen = Group("AlpgenPythia_L",[LeftAlpgenNp0,LeftAlpgenNp1,LeftAlpgenNp2,LeftAlpgenNp3,LeftAlpgenNp4,LeftAlpgenNp5incl]) 
RightAlpgen = Group("AlpgenPythia_R",[RightAlpgenNp0,RightAlpgenNp1,RightAlpgenNp2,RightAlpgenNp3,RightAlpgenNp4,RightAlpgenNp5incl]) 
LeftAlpgen.weight = 'evtsel_weight_ZPt_OS'
RightAlpgen.weight = 'evtsel_weight_ZPt_OS'
LeftAlpgen.xpol = 0.5*(1-Alpol)
RightAlpgen.xpol = 0.5*(1+Alpol)


LeftAlpgen_RW = LeftAlpgen.clone("AlpgenPythia_L_Pythia8_trueta"  )
RightAlpgen_RW = RightAlpgen.clone("AlpgenPythia_R_Pythia8_trueta"  )
LeftAlpgen_RW.weight = 'evtsel_weight_ZPt_OS*weight_Pythia8_trupt'
RightAlpgen_RW.weight = 'evtsel_weight_ZPt_OS*weight_Pythia8_trupt'

# Pythia8 Tau Spinner Polarization
p8pol = -0.1290

LeftPythia8 = Sample(name="Pythia8_L", DSID = 147818, version = 'LEFT', xsec = 974.42*1.18*1, xpol = 0.5*(1-p8pol), theory = True)
RightPythia8 = Sample(name="Pythia8_R", DSID = 147818, version = 'RIGHT', xsec = 974.42*1.18*1, xpol = 0.5*(1+p8pol), theory = True)
LeftPythia = Group( 'Pythia8_L', [LeftPythia8])
RightPythia = Group( 'Pythia8_R', [RightPythia8])
LeftPythia.xpol = 0.5*(1-p8pol)
RightPythia.xpol = 0.5*(1+p8pol)

# Powheg Truth Samples

LeftPowheg0 = Sample(name="Powheg_L", DSID = 147808, version = 'LEFT', xsec = 974.42*1.18*1, theory = True)
RightPowheg0 = Sample(name="Powheg_R", DSID = 147808, version = 'RIGHT', xsec = 974.42*1.18*1, theory = True)
LeftPowheg = Group( 'Powheg_L', [LeftPowheg0])
RightPowheg = Group( 'Powheg_R', [RightPowheg0])

# Herwig Truth Samples

LeftHerwigNp0 = Sample(name="LeftHerwigNp0", DSID = 107670, version = 'LEFT', xsec = 718.87*1.18*1.0, xpol = 0.5*(1-p0), theory = True)
LeftHerwigNp1 = Sample(name="LeftHerwigNp1", DSID = 107671, version = 'LEFT', xsec = 175.75*1.18*1.0, xpol = 0.5*(1-p1), theory = True)
LeftHerwigNp2 = Sample(name="LeftHerwigNp2", DSID = 107672, version = 'LEFT', xsec = 58.856*1.18*1.0, xpol = 0.5*(1-p2), theory = True)
LeftHerwigNp3 = Sample(name="LeftHerwigNp3", DSID = 107673, version = 'LEFT', xsec = 15.667*1.18*1.0, xpol = 0.5*(1-p3), theory = True)
LeftHerwigNp4 = Sample(name="LeftHerwigNp4", DSID = 107674, version = 'LEFT', xsec = 4.0121*1.18*1.0, xpol = 0.5*(1-p4), theory = True)
LeftHerwigNp5incl = Sample(name="LeftHerwigNp5incl", DSID = 107675, version = 'LEFT', xsec = 1.256*1.18*1.0, xpol = 0.5*(1-p5), theory = True)

RightHerwigNp0 = Sample(name="RightHerwigNp0", DSID = 107670, version = 'RIGHT', xsec = 718.87*1.18, xpol =  0.5*(1+p0), theory = True)
RightHerwigNp1 = Sample(name="RightHerwigNp1", DSID = 107671, version = 'RIGHT', xsec = 175.75*1.18, xpol =  0.5*(1+p1), theory = True)
RightHerwigNp2 = Sample(name="RightHerwigNp2", DSID = 107672, version = 'RIGHT', xsec = 58.856*1.18, xpol =  0.5*(1+p2), theory = True)
RightHerwigNp3 = Sample(name="RightHerwigNp3", DSID = 107673, version = 'RIGHT', xsec = 15.667*1.18, xpol =  0.5*(1+p3), theory = True)
RightHerwigNp4 = Sample(name="RightHerwigNp4", DSID = 107674, version = 'RIGHT', xsec = 4.0121*1.18, xpol =  0.5*(1+p4), theory = True)
RightHerwigNp5incl = Sample(name="RightHerwigNp5incl", DSID = 107675, version = 'RIGHT', xsec = 1.256*1.18, xpol =  0.5*(1+p5), theory = True)

LeftHerwig = Group("Herwig_L",[LeftHerwigNp0,LeftHerwigNp1,LeftHerwigNp2,LeftHerwigNp3,LeftHerwigNp4,LeftHerwigNp5incl]) 
RightHerwig = Group("Herwig_R",[RightHerwigNp0,RightHerwigNp1,RightHerwigNp2,RightHerwigNp3,RightHerwigNp4,RightHerwigNp5incl]) 


###############################################################################################################################################
########## ANALYSIS GROUPS ####################################################################################################################
###############################################################################################################################################

# Signal Ztautau_______________________________________________________________________________________________________________________________
Ztautau = Group("Ztautau",[AlpPythZtautauNp0,AlpPythZtautauNp1,AlpPythZtautauNp2,AlpPythZtautauNp3,AlpPythZtautauNp4,AlpPythZtautauNp5incl])
Ztautau.fill = ROOT.kAzure-4; 
Ztautau.legend = 'Z#rightarrow#tau#tau'
Ztautau.weight = 'evtsel_weight_ZPt_OS' # includes Z pT reweighting 
Ztautau.latex = '\Ztautau'
## Left-Handed
LeftZtautau = Group("LH_Ztautau",[LeftZtautauNp0,LeftZtautauNp1,LeftZtautauNp2,LeftZtautauNp3,LeftZtautauNp4,LeftZtautauNp5incl]) 
LeftZtautau.fill=ROOT.kRed+1
LeftZtautau.legend = 'LH Z#rightarrow#tau#tau'
LeftZtautau.weight = 'evtsel_weight_ZPt_OS' # includes Z pT reweighting 
LeftZtautau.latex = 'Left-Handed'
## Decay Mode Decomposition 
LHNP = LeftZtautau.clone('LH_NP', string = '(!tau_SP)')
LHRHO = LeftZtautau.clone('LH_RHO', string = '(tau_RHO)')
LHPINU = LeftZtautau.clone('LH_PINU', string = '(tau_PINU)')
LHSPX = LeftZtautau.clone('LH_SPX', string = '(tau_SPX)')
## Right-Handed
RightZtautau = Group("RH_Ztautau",[RightZtautauNp0,RightZtautauNp1,RightZtautauNp2,RightZtautauNp3,RightZtautauNp4,RightZtautauNp5incl]) 
RightZtautau.fill=ROOT.kBlue+1
RightZtautau.legend = 'RH Z#rightarrow#tau#tau'
RightZtautau.weight = 'evtsel_weight_ZPt_OS'
RightZtautau.latex = 'Right-Handed'
## Decay Mode Decomposition
RHNP = RightZtautau.clone('RH_NP', string = '(!tau_SP)')
RHRHO = RightZtautau.clone('RH_RHO', string = '(tau_RHO)')
RHPINU = RightZtautau.clone('RH_PINU', string = '(tau_PINU)')
RHSPX = RightZtautau.clone('RH_SPX', string = '(tau_SPX)')

## Pythia 8 Signal
LeftP8Ztautau = Group("LH_Ztautau", [LeftP8], fill =2)
LeftP8Ztautau.legend = 'LH Z#tau#tau'

RightP8Ztautau = Group("RH_Ztautau", [RightP8], fill = 4)
RightP8Ztautau.legend = 'RH Z#tau#tau'

P8Ztautau = Group("Ztautau", [ZtautauP8], fill = ROOT.kAzure-4)
P8Ztautau.legend = 'Z#tau#tau'

# Wlnu + jets__________________________________________________________________________________________________________________________________
AlpgenPythiaWmunu = Group("Wmunu",[AlpgenPythiaWmunuNp0,AlpgenPythiaWmunuNp1,AlpgenPythiaWmunuNp2,AlpgenPythiaWmunuNp3,AlpgenPythiaWmunuNp4,AlpgenPythiaWmunuNp5incl])
AlpgenPythiaWenu = Group("Wenu",[AlpgenPythiaWenuNp0,AlpgenPythiaWenuNp1,AlpgenPythiaWenuNp2,AlpgenPythiaWenuNp3,AlpgenPythiaWenuNp4,AlpgenPythiaWenuNp5incl])
AlpgenPythiaWtaunu = Group("Wtaunu",[AlpgenPythiaWtaunuNp0,AlpgenPythiaWtaunuNp1,AlpgenPythiaWtaunuNp2,AlpgenPythiaWtaunuNp3,AlpgenPythiaWtaunuNp4,AlpgenPythiaWtaunuNp5incl])
Wlnu = AlpgenPythiaWenu + AlpgenPythiaWmunu + AlpgenPythiaWtaunu
Wlnu.name = 'Wlnu'; Wlnu.legend = 'W+jets'; Wlnu.fill = 5
Wlnu.latex = '$W$+jets'

# Zll + jets___________________________________________________________________________________________________________________________________
AlpgenPythiaZee = Group("Zee",[AlpgenPythiaZeeNp0,AlpgenPythiaZeeNp1,AlpgenPythiaZeeNp2,AlpgenPythiaZeeNp3,AlpgenPythiaZeeNp4,AlpgenPythiaZeeNp5incl])
AlpgenPythiaZmumu = Group("Zmumu",[AlpgenPythiaZmumuNp0,AlpgenPythiaZmumuNp1,AlpgenPythiaZmumuNp2,AlpgenPythiaZmumuNp3,AlpgenPythiaZmumuNp4,AlpgenPythiaZmumuNp5incl])
Zll = AlpgenPythiaZee + AlpgenPythiaZmumu
Zll.name = 'Zll'; Zll.legend = 'Z#rightarrow ll'
Zll.fill = ROOT.kOrange+7
Zll.weight = 'evtsel_weight_ZPt_OS'

# Z+jets background w. lepton faking tau
Zlltau = Zll.clone('Zltau',string = '(evtsel_tau_is_mu || evtsel_tau_is_el_parentEW)')
Zlltau.legend = 'Z+jets(e/#mu as #tau)'
Zlltau.latex = '\Zll($\ell\\to\\tau$)'
Zlltau.fill = ROOT.kPink+10
Zlltau.weight = 'evtsel_weight_ZPt_OS'
# Z+jets background w. jet faking tau
Zlljet = Zll.clone('Zlljet',string = '!(evtsel_tau_is_mu || evtsel_tau_is_el_parentEW)')
Zlljet.legend = 'Z+jets(jet as #tau)'
Zlljet.latex = '\Zll+jet'
Zlljet.fill = ROOT.kMagenta+3
Zlljet.weight = 'evtsel_weight_ZPt_OS'

# ttbar________________________________________________________________________________________________________________________________________
ttBar = Group("ttBar",[Pythia_ttbar])
ttBar.fill = ROOT.kCyan-4
ttBar.legend = 't#bar{t}'
ttBar.latex = '\\ttbar'

# data_________________________________________________________________________________________________________________________________________
Muons = Group("Muons", [Muons_PeriodA,Muons_PeriodB,Muons_PeriodC,Muons_PeriodD,Muons_PeriodE,Muons_PeriodG,
                        Muons_PeriodH,Muons_PeriodI,Muons_PeriodJ,Muons_PeriodL], isData = True)
Egamma = Group("Egamma", [Egamma_PeriodA,Egamma_PeriodB,Egamma_PeriodC,Egamma_PeriodD,Egamma_PeriodE,Egamma_PeriodG,
                          Egamma_PeriodH,Egamma_PeriodI,Egamma_PeriodJ,Egamma_PeriodL], isData = True)
Muons.latex = 'Data'
Egamma.latex = 'Data'
Data = Muons + Egamma 
Data.name = 'Data'

###############################################################################################################################################
########## STUDY GROUPS #######################################################################################################################
###############################################################################################################################################

# studies for tau fakes in W + jets____________________________________________________________________________________________________________
Wlnu_bJet  = Wlnu.clone('Wlnu_bjet',string = '(evtsel_tau_is_bjet)', fill = 6)
Wlnu_cJet  = Wlnu.clone('Wlnu_cjet',string = '(evtsel_tau_is_cjet)', fill = 6)
Wlnu_Jet   = Wlnu.clone('Wlnu_jet',string = '(evtsel_tau_is_jet)', fill = 6) 
Wlnu_Quark = Wlnu.clone('Wlnu_LightJet',string = '(evtsel_tau_is_lightjet)', fill = 6) 
Wlnu_Gluon = Wlnu.clone('Wlnu_GluJet',string = '(evtsel_tau_is_glujet)', fill = 6) 
Wlnu_real  = Wlnu.clone('Wlnu_realTau',string = '(evtsel_tau_is_real)', fill = 6) 
Wlnu_El    = Wlnu.clone('Wlnu_el',string = '(evtsel_tau_is_el)', fill = 6) 
Wlnu_Mu    = Wlnu.clone('Wlnu_mu',string = '(evtsel_tau_is_mu)', fill = 6) 

# studies for real taus in ttbar control region________________________________________________________________________________________________
ttBar_real = ttBar.clone('ttbar_real',string = '(evtsel_tau_is_real)')
ttBar_real.cut = Cut('evtsel_tau_is_real')
ttBar_real.legend = 't#bar{t} (true tau)'
ttBar_fake = ttBar.clone('ttbar_fake',string = '(!evtsel_tau_is_real)')
ttBar_fake.cut = Cut('!evtsel_tau_is_real')
ttBar_fake.legend = 't#bar{t} (fake tau)'
ttBar_fake.fill = ROOT.kCyan+2

# Np Sample comparison studies________________________________________________________________________________________________________________

Np0 = Group("Np0", [AlpPythZtautauNp0], fill = ROOT.kAzure-4)
Np1 = Group("Np1", [AlpPythZtautauNp1], fill = ROOT.kAzure-4)
Np2 = Group("Np2", [AlpPythZtautauNp2], fill = ROOT.kAzure-4)
Np3 = Group("Np3", [AlpPythZtautauNp3], fill = ROOT.kAzure-4)
Np4 = Group("Np4", [AlpPythZtautauNp4], fill = ROOT.kAzure-4)
Np5 = Group("Np5", [AlpPythZtautauNp5incl], fill = ROOT.kAzure-4)
Np0.legend = 'Z#tau#tau Np0'

Left_Np0 = Group("Left_Np0", [LeftZtautauNp0], fill = 2)
Left_Np1 = Group("Left_Np1", [LeftZtautauNp1], fill = 49)
Left_Np2 = Group("Left_Np2", [LeftZtautauNp2], fill = 48)
Left_Np3 = Group("Left_Np3", [LeftZtautauNp3], fill = 47)
Left_Np4 = Group("Left_Np4", [LeftZtautauNp4], fill = 46)
Left_Np5 = Group("Left_Np5", [LeftZtautauNp5incl], fill = 45)

Right_Np0 = Group("Right_Np0", [RightZtautauNp0], fill = 4)
Right_Np1 = Group("Right_Np1", [RightZtautauNp1], fill = 39)
Right_Np2 = Group("Right_Np2", [RightZtautauNp2], fill = 38)
Right_Np3 = Group("Right_Np3", [RightZtautauNp3], fill = 37)
Right_Np4 = Group("Right_Np4", [RightZtautauNp4], fill = 36)
Right_Np5 = Group("Right_Np5", [RightZtautauNp5incl], fill = 35)

# Z mass peak Studies_______________________________________________________________________________________________________________________

Ztautau_Off = Ztautau.clone('Ztautau_Off', string = '(evtsel_true_mtt < 80 || evtsel_true_mtt > 100)', fill = ROOT.kViolet-5)
Ztautau_Off.legend = 'Z#tau#tau (off Z peak)'
Ztautau_lo = Ztautau.clone('Ztautau_lo', string = '(evtsel_true_mtt < 80)', fill = ROOT.kViolet-3)
Ztautau_lo.legend = 'Z#tau#tau (m_{#tau#tau} < 80)'
Ztautau_hi = Ztautau.clone('Ztautau_hi', string = '(evtsel_true_mtt > 100)', fill = ROOT.kViolet-7)
Ztautau_hi.legend = 'Z#tau#tau (m_{#tau#tau} > 100)'

Ztautau_On = Ztautau.clone('Ztautau_On', string = '(evtsel_true_mtt >= 80 && evtsel_true_mtt <= 100)')
Ztautau_On.legend = 'Z#tau#tau (on Z peak)'
LeftZtautau_On = LeftZtautau.clone('LH_Ztautau_On', string = '(evtsel_true_mtt >= 80 && evtsel_true_mtt <= 100)')
LeftZtautau_On.legend = 'LH Z#tau#tau (on Z peak)'
RightZtautau_On = RightZtautau.clone('RH_Ztautau_On', string = '(evtsel_true_mtt >= 80 && evtsel_true_mtt <= 100)')
RightZtautau_On.legend = 'RH Z#tau#tau (on Z peak)'

# Pythia8 Signal
P8Ztautau_Off = P8Ztautau.clone('P8Ztautau_Off', string = '(evtsel_true_mtt < 71 || evtsel_true_mtt > 111)', fill = ROOT.kViolet-5)
P8Ztautau_Off.legend = 'Z#tau#tau (off Z peak)'
P8Ztautau_On = P8Ztautau.clone('P8Ztautau_On', string = '(evtsel_true_mtt > 71 && evtsel_true_mtt < 111)')
P8Ztautau_On.legend = 'Z#tau#tau (on Z peak)'
LeftP8Ztautau_On = LeftP8Ztautau.clone('LH_P8Ztautau_On', string = '(evtsel_true_mtt > 71 && evtsel_true_mtt < 111)')
LeftP8Ztautau_On.legend = 'LH Z#tau#tau (on Z peak)'
RightP8Ztautau_On = RightP8Ztautau.clone('RH_P8Ztautau_On', string = '(evtsel_true_mtt > 71 && evtsel_true_mtt < 111)')
RightP8Ztautau_On.legend = 'RH Z#tau#tau (on Z peak)'


