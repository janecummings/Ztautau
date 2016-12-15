import os,sys
import hist
import run
import ROOT 
from ErrorFloat import ErrorFloat
from Drawer import *
import draw
from printing import row


def make_sure_dir(s):
    if not os.path.exists(s.plotdir): os.mkdir(s.plotdir)
    if not os.path.exists(s.plotdir+'/control'): os.mkdir(s.plotdir+'/control')
    if not os.path.exists(s.plotdir+'/modes'): os.mkdir(s.plotdir+'/modes')
    if not os.path.exists(s.plotdir+'/syst'): os.mkdir(s.plotdir+'/syst')


def plot_SR(s,d, plots = None , syst = None, crspec = None):
    stream = s.stream
    make_sure_dir(s)
    #print_scale_factors(s)
    if plots == None: plots = s.reg_entries
    for entry in plots:
        entry.blind = False
        for group in s.groups: entry.add_group(group)
        #entry.rebin = entry.rebin * 10
        if 'Q' in d: draw.compare_QCD_shapes(s,entry)
        if 'S' in d: draw.plot_SS(s,entry,'SR_SS')
        if 'W' in d: draw.plot_stat(s,entry,'WCR_OS', k_entry = upsilon, verbose = False, syst = syst)
        
        if entry.name in ['upsilon','ftrk','tau_leadTrkPt']: 
           entry.blind = True
                
        if 'M' in d: 
            draw.plot_modes(s,entry,'SR_OS', verbose = False, hand = False)
            draw.plot_modes(s,entry,'SR_OS', verbose = False, hand = 'LH')
            draw.plot_modes(s,entry,'SR_OS', verbose = False, hand = 'RH')
        if 'D' in d: 
            #draw.plot_stat_qcd(s,entry,'SR_OS', verbose = False, syst = syst, noHands = False)
            #draw.plot_stat_qcd(s,entry,'SR_OS', verbose = True, syst = syst)
            draw.plot_stat_qcd(s,entry,'SR_OS', verbose = False, syst = syst)

            #draw.plot_stat_qcd(s,entry,'SR_OS_%s_UP' % syst, verbose = True)
            #draw.plot_stat_qcd(s,entry,'SR_OS_%s_DOWN' % syst, verbose = True)
            #entry.blind = False
            #draw.plot_stat_qcd(s,entry,'SR_SS', verbose = True, syst = syst)
            #draw.plot_QCD_sys(s,entry,syst=syst  )
        if 'R' in d: draw.plot_stat(s,entry,'SR_OS', verbose = False, syst = syst, crspec = crspec)
        if 'H' in d: 
            draw.plot_LR(s,entry,'SR_OS', verbose = False, hand = 'LR', norm= True)
            #draw.plot_LR(s,entry,'SR_OS', verbose = False, hand = 'OT', norm= False)
            #draw.plot_LR(s,entry,'SR_OS', verbose = False, hand = 'L', norm= True, comp = run.mu_fit)
            #draw.plot_LR(s,entry,'SR_OS', verbose = False, hand = 'R', norm= True, comp = run.mu_fit)            

            #draw.plot_stat_qcd(s,entry,'SR_OS', verbose = False, noHands=False , syst = syst)
    if 'T' in d: 
        for entry in [transMass,deltaPhi]:
            for group in s.groups: entry.add_group(group)
            draw.plot_stat(s,entry,'6_SP_OS', k_entry = upsilon, verbose = False, syst = syst)

def plot_WCR(s):
    for entry in s.cr_entries:
        for group in s.groups: entry.add_group(group)
        draw.plot_CR(s,entry,'WCR_OS', verbose = False)
        draw.plot_CR(s,entry,'WCR_SS', verbose = False)

def print_scale_factors(s, SF='QW'):

    stream = s.stream
    kCR = ROOT.TFile('%s/%s_regions.root' % (s.OUTPUTDIR,stream.name))
    entry = upsilon
    kW_OS = ErrorFloat(1.,0.); kW_SS = ErrorFloat(1.,0.)
    print '-' *200
    if 'W' in SF:
        kW_OS = hist.calc_kX(kCR,stream,entry,'WCR','Wlnu',OS='OS',verbose=True)
        print '-' *200
        kW_SS = hist.calc_kX(kCR,stream,entry,'WCR','Wlnu',OS='SS',verbose=True)
    if 'Q' in SF:
        hist.calc_rQCD(kCR,stream,entry,'QCD',kW_OS = kW_OS, kW_SS = kW_SS, verbose = True)
    print '-' *200

def print_np(s):
    stream = s.stream
    make_sure_dir(s)
    for entry in s.reg_entries:
        for group in s.groups: entry.add_group(group)
        draw.plot_Np(s,entry,'SR_OS',hand='L')
        draw.plot_Np(s,entry,'SR_OS',hand='R')  



if __name__ == '__main__':

    plot_SR(run.el_zpt,'D')
    #plot_SR( run.mu_p8_plots, 'H'  )
    #for s in [run.mu_tight_plots,run.el_tight_plots]:
    #    plot_SR(s, 'W') #, plots = [met,lep_pt])
    #print_np(run.mu_NpX)
    # print row('Systematic',['UP','DOWN'],fmt='twiki')
    # print row('Systematic',['Left Eff UP','Right Eff UP', 'Left Eff DOWN', 'Right Eff DOWN', 'UP','DOWN'],fmt='twiki')
    # for sys in ['TES_TOTAL','ElES_LowPt','ElES_PS','ElES_R12','ElES_Zee','ElEnResSys','JER','JES_BJet','JES_Detector1','JES_EtaMethod',
    #                       'JES_EtaModelling','JES_FlavComp','JES_FlavResp','JES_Modelling1','JES_PUMu','JES_PUMu','JES_PUNPV','JES_PURho',
    #                       'JES_Statistical1','JVF','METResSys','METScaleSys','MuSys','TES_CLOSURE','TES_FINAL','TES_INSITUINTERPOL',
    #                       'TES_MODELING','TES_SINGLEPARTICLEINTERPOL']:
    #     plot_SR(run.mu_fit, 'D', syst = sys)
        #plot_SR(run.mu_fit, 'D', syst = sys)
    #['el_trig','el_id','el_iso','mu_id','mu_iso','tau_el','tau_fake','tau_id_stat','tau_id_sys']: #'PU_rescaling'

    # for sys in ['mu_trig','el_id','el_iso','mu_id','mu_iso','tau_el','tau_fake','tau_id_stat','tau_id_sys']: #'PU_rescaling'
    #     plot_SR(run.mu_fit, 'D', syst = sys)
    #     #plot_SR(run.el_fit, 'D', syst = sys)

    #print_scale_factors(run.mu_cutflow)


