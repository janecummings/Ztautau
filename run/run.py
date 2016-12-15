from spec import spec
from Drawer import *
#import streams,Groups

## Systematic categories
tausys = ['tau_id_sys','tau_id_stat', 'tau_el'] # SF
lepsys = ['el_id','el_iso','mu_iso','mu_trig','el_trig'] # SF
emsys = ['ElES_LowPt','ElES_PS','ElES_R12','ElES_Zee','ElEnResSys' ] # branches
JES = ['JER','JES_Detector1','JES_EtaMethod',
                     'JES_EtaModelling','JES_FlavComp','JES_FlavResp','JES_Modelling1','JES_PUNPV','JES_PURho'] # branch
MET = ['METResSys','METScaleSys','METRESOSOFT','METSCALESOFT'] # branches 
PU = ['PU_rescaling'] 

################################################################
# Standard Event Selection
################################################################
el_fit = spec( streams.el,
            DIR = 'prod4/el_fit',
            regions = ['SR'],
            control = ['WCR','QCD'],
            reg_entries = [upsilon_20],
            cr_entries = [upsilon_CR],
            syst = emsys + JES + MET,
            SFsys = tausys + lepsys + PU
            )


mu_fit = spec( streams.mu,
            DIR = 'prod4/mu_fit',
            regions = ['SR'],
            control = ['WCR','QCD'],
            reg_entries = [upsilon_20],
            cr_entries = [upsilon_CR],
            syst = emsys + JES + MET,
            SFsys = tausys + lepsys + PU
            )

mu_theory  = spec( streams.mutheory,
    DIR = 'prod4/theory',
    regions = ['SR'],
    control = None,
    reg_entries = [upsilon_20],
    #reg_entries = [trueTau_upsilon, trueTau_fEM, trueTau_decayMode, truethad_pt, truethad_EMfrac, truethad_eta, truethad_phi, trueZ_m, trueZ_pt, true_mtt],  
    groups = streams.mutheory.Alpgen
    )

el_theory  = spec( streams.eltheory,
    DIR = 'prod4/theory',
    regions = ['SR'],
    control = None,
    reg_entries = [upsilon_20],
    #reg_entries = [trueTau_upsilon, trueTau_fEM, trueTau_decayMode, truethad_pt, truethad_EMfrac, truethad_eta, truethad_phi, trueZ_m, trueZ_pt, true_mtt],
    groups = streams.eltheory.Alpgen
    )


mu_2d = spec( streams.mu,
  DIR = 'prod4/leadpt',
  regions = ['SR'],
  control = None,
  reg_entries = [ups_met],
  groups = [streams.mu.signal])


el_2d = spec( streams.el,
  DIR = 'prod4/leadpt',
  regions = ['SR'],
  control = None,
  reg_entries = [ups_met],
  groups = [streams.el.signal])

mu_ttbar = spec( streams.mu_ttbar,
            DIR = 'prod4/mu_ttbar',
            regions = ['TTBAR'],
            control = ['QCD'],
            reg_entries = [tau_leadTrkPt, upsilon_20],
            cr_entries = [upsilon_CR],
            groups = [streams.mu_ttbar.data, streams.mu_ttbar.signal] + streams.mu_ttbar.MC)

el_ttbar = spec( streams.el_ttbar,
            DIR = 'prod4/el_ttbar',
            regions = ['TTBAR'],
            control = ['QCD'],
            reg_entries = [upsilon_20,tau_leadTrkPt],
            cr_entries = [upsilon_CR],
            groups = [streams.el_ttbar.data, streams.el_ttbar.signal] + streams.el_ttbar.MC)

## For Validation
mu_note = spec( streams.mu,
            DIR = 'prod4/bch',
            regions = ['WCR'],
            #control = ['WCR','QCD'],
            control = None,
            reg_entries = [upsilon_20],
            cr_entries = [upsilon_CR])

el_note = spec( streams.el,
            DIR = 'prod4/bch',
            regions = ['WCR'],
            #control = ['WCR','QCD'],
            control = None,
            reg_entries = [upsilon_20],
            cr_entries = [upsilon_CR])

