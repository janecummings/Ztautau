import ROOT
import sys
from printing import row
import streams
import hist
from ErrorFloat import ErrorFloat
import math
import Drawer
import Groups
import pdg_simple

fcol=20;col=22

modelegend = {"0":"All",  "1":"Single-Prong","2":"Three-Prong","3":"$\\tau\\to\\pi\\nu$","4":"$\\tau\\to\\pi\\pi^{0}\\nu$",
'5': "$\\tau\\to\\pi N\\pi^{0}\\nu$"}
cuts = ['Preselected', 'BCH Clean', 'Opposite Sign', 'emu,$\pt>26~\GeV$', 'Isolated emu', 'Tau', 'Dilepton Veto','Tau Lepton Vetoes', 'Single-Prong', 
    '$m_{\mathrm{T}}<30~\GeV$', '$\sum\Delta\phi < 3.5$', 'Visible Mass', '$\Upsilon<1.5$']
pol = -0.141259124591
sm_lh = (1-pol)/2.
sm_rh = (1+pol)/2.

OUTPUT_DIR = '~/panalysis/output/prod4/bch'


def table_Signal_Samples():
    col = 15
    fcol = 15
    fmt = "latex"
    xsec_total = 0
    g = Groups.Ztautau
    gL = Groups.LeftZtautau
    gR = Groups.RightZtautau
    Nl = ErrorFloat(0,0)
    Nr = ErrorFloat(0,0)
    for s in g:
        xsec_total+=s.xsec

    NL = 0
    NR = 0
    pol = 0

    nP = ["Np0","Np1","Np2","Np3","Np4","Np5"]
    for nPi,gi,gLi,gRi in zip(nP,g,gL,gR):
        gi.load(); giN = gi.totalEvents; gi.close()
        gLi.load(); gLiN = gLi.totalEvents; gLi.close()
        gRi.load(); gRiN = gRi.totalEvents; gRi.close()

        nL = ErrorFloat(gLiN,math.sqrt(gLiN))
        nR = ErrorFloat(gRiN,math.sqrt(gRiN))

        pol += (gRiN -gLiN)/(giN)  *gi.xsec/xsec_total

        NL += gLiN*gi.xsec/xsec_total/(gLiN+gRiN)
        NR += gRiN*gi.xsec/xsec_total/(gLiN+gRiN)

        Nl += nL*ErrorFloat(gi.xsec/xsec_total,0)/(nL+nR)
        Nr += nR*ErrorFloat(gi.xsec/xsec_total,0)/(nL+nR)
        
        #print nPi, (gRiN -gLiN)/(giN)
        print row(nPi, [int(giN),int(gLiN),int(gRiN),gi.xstr],fcol=fcol,col=col,fmt=fmt)

    # print Nr,Nl
    # print (Nr-Nl)/(Nr+Nl)
    # print pol
    # print (NR-NL)/(NR+NL)


def table_EW_Samples():
	col = 15; fcol = 15; fmt = "latex"
	groups = [Groups.AlpgenPythiaWenu,Groups.AlpgenPythiaWmunu,Groups.AlpgenPythiaWtaunu,Groups.AlpgenPythiaZee,Groups.AlpgenPythiaZmumu]
	nP = ["Np0","Np1","Np2","Np3","Np4","Np5"]
	for g in groups:
		for nPi,gi in zip(nP,g):
			gi.load()
			print row("\\%s+%sjets"%(g.name,Npi[-1]), [int(gi.totalEvents), gi.xstr], fcol=fcol,col=col,fmt=fmt)
			gi.close()

	tt = Groups.Powheg_ttbar
	tt.load()
	print row("\\ttbar", [tt.totalEvents, tt.xstr],fcol=fcol,col=col,fmt="latex" )
	tt.close()


def cutflow(mu, pdg = True, smnorm = False):
        rf = ROOT.TFile('%s/%s_regions.root' % (OUTPUT_DIR, mu.name))

        eff = False
        do_pol = False
        f = False
        fcol=25;col=22;fmt="latex"
        entry = Drawer.upsilon

        regions = mu.noteflow
        pt = -0.144
        groups = [mu.data, mu.signal, mu.LH, mu.RH,  mu.ttbar, mu.Zlljet, mu.Zltau, mu.Wlnu]

        a = []
        for i,region in enumerate(regions):
            cutgroup = []
            n = []
            mc = 0
            nL = None
            for group in groups:
                if group.isData:
                    h = rf.Get('h_%s_%s_%s' % ('upsilon',group.name,region.name))
                else:
                    h = rf.Get('h_%s_%s_%s' % (entry.name,group.name,region.name))
                #print 'h_%s_%s_%s' % (entry.name,group.name,region.name)
                norm = h.Integral(0,h.GetNbinsX()+1)
                if group.isData:
                    cutgroup.append(int(norm))
                else:
                    v,e = hist.IntegralError(h, of = True)
                    if smnorm:
                        if 'LH' in group.name: 
                            v = v*sm_lh
                            e = e*sm_lh
                        if 'RH' in group.name:
                            v = v*sm_rh
                            e = e*sm_rh                            
                    if pdg:
                        cutgroup.append( '%i +/- %i' % pdg_simple.pdg_round(v,e))
                    else:
                        cutgroup.append( '%.2f +/- %.2f' % (v,e))
                n.append(h.Integral(0,h.GetNbinsX()+1))
                if 'L' in group.name: nL = norm
                elif 'R' in group.name: nR = norm
                if group.name not in ['Muons','Egamma','LH_Ztautau','RH_Ztautau']: 
                    #print group.name
                    mc += h.Integral(0,h.GetNbinsX()+1)
            if i == 0: 
                a = list(n); l=a
                print row('', [group.name for group in groups], fcol =fcol, col = col)
            if nL and do_pol:
                nL = 0.5*(1-pt)*nL
                nR = 0.5*(1+pt)*nR
                print row(cuts[i]+'\t\t',cutgroup + [(nR-nL)/(nR+nL)],fcol=fcol,col=col, fmt = fmt)
            else:
                #print row(region.name , cutgroup , fcol=fcol,col=col, fmt = fmt)
                print row(cuts[i], cutgroup , fcol=fcol,col=col, fmt = fmt)

            #print (nR-nL)/nR
            if f:
                print row(region.name, [''] + [x/mc for x in n[1:]],fcol=fcol,col=col,fmt=fmt)
            if eff:
                print row('',[x/y for (x,y) in zip(n,a)],fcol=fcol,col=col, fmt = fmt)
                print row('',[x/y for (x,y) in zip(n,l)],fcol=fcol,col=col, fmt = fmt)
            l = n



def backgrounds(mu):
    rf = ROOT.TFile('%s/%s_regions.root' % (OUTPUT_DIR,mu.name))
    #rf = ROOT.TFile('~/panalysis/output/%s_Zpt_NoWeight/%s_regions.root' % (mu.name,mu.name))

    eff = True
    f = False
    fcol=15;col=15;fmt="latex"
    entry = Drawer.visMass
    k_entry = Drawer.upsilon_CR
    regions = mu.cutflow
    #regions = mu.SR
    verbose = True
    groups = [mu.Zltau,mu.Zlljet,mu.ttbar,mu.Wlnu]
    kW_OS = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose)

    a = []
    for i,region in enumerate(regions):
        cutgroup = []
        n = []
        mc = 0
        for group in groups:
            h = rf.Get('h_%s_%s_%s' % (entry.name,group.name,region.name))
            if 'W' not in group.name:
                cutgroup.append( '%i +/- %i' % hist.IntegralError(h))
                n.append(h.Integral(0,h.GetNbinsX()+1))
            else:
                nW =  ErrorFloat(*hist.IntegralError(h))*kW_OS
                cutgroup.append('%i +/- %i' % (nW.val,nW.err))
                n.append(nW.val)


            if group.name not in ['Muons','Egamma','LH_Ztautau','RH_Ztautau']: 
                mc += h.Integral(0,h.GetNbinsX()+1)
            if i == 0: a = list(n); l=a
        else: print row(cuts[i]+'\t\t',cutgroup,fcol=fcol,col=col, fmt = fmt)
        if f:
            print row(region.name, [''] + [x/mc for x in n[1:]],fcol=fcol,col=col,fmt=fmt)
        if eff:
            print row('',[x/y for (x,y) in zip(n,a)],fcol=fcol,col=col, fmt = fmt)
            print row('',[x/y for (x,y) in zip(n,l)],fcol=fcol,col=col, fmt = fmt)
        l = n


def print_region(mu,reg, pdg=False):#,regions,groups,entry=Drawer.upsilon_CR,eff=False,twiki=False, fmt = None, f = False):
    a = []
    rf = ROOT.TFile('%s/%s_regions.root' % (OUTPUT_DIR, mu.name))
    regions = getattr(mu,reg)
    groups = [mu.data, mu.signal, mu.LH, mu.RH,  mu.ttbar, mu.Zlljet, mu.Zltau, mu.Wlnu]
    entry = Drawer.upsilon_CR
    fmt =  ''
    if pdg: fmt = 'latex'

    print row('',[group.name for group in groups],fcol=fcol,col=col,fmt=fmt)
    for i,region in enumerate(regions):
        cutgroup = []
        n = []
        mc = 0
        nEW = 0
        for group in groups:
            if 'Wlnu' in group.name:
                wreg = region.name.replace('QCD','WCR')
                wreg = wreg.replace('SR','WCR')
                #print region
                h = hist.get_dataW(rf,mu,entry, wreg, where = reg)
                #print h.Integral()
            else:
                h = rf.Get('h_%s_%s_%s' % (entry.name,group.name,region.name))
            norm = h.Integral(0,h.GetNbinsX()+1)
            if group.isData:
                cutgroup.append(int(norm))
            else:
                #print h.GetName()
                v,e = hist.IntegralError(h, of = True)
                if pdg:
                    cutgroup.append( '%g +/- %g' % pdg_simple.pdg_round(v,e))
                else:
                    cutgroup.append( '%.2f +/- %.2f' % (v,e))
            
            n.append(norm)
            if group.isData: 
                nData = norm
            elif 'LH' in group.name:
                nL = norm
            elif 'RH' in group.name:
                nR = norm
            elif 'Ztautau' in group.name:
                nS = norm
            elif 'Wlnu' not in group.name:
                nEW += norm

            if i == 0: a = list(n); l=a
        #if twiki: print twikirow(region.name,cutgroup)
        print row(region.name,cutgroup,fcol=fcol,col=col, fmt = fmt)

        # print '%s \t\t %s \t\t %s \t\t %s \t\t %s \t\t %s' % ('NS', 'NL', 'NR', 'W EST (NOM)', 'W EST (LH)', 'W EST (RH)')
        # print '%i \t\t %i \t\t %i \t\t %i \t\t %i \t\t %i' % (nS, nL, nR, nData - nS - nEW, nData - nL - nEW, nData - nR - nEW)
        # print nL/(nData - nL - nEW)
        # print nR/(nData - nR - nEW)
        # print nS/(nData - nS - nEW)
        # break


def signal_region(mu):
    #rf = ROOT.TFile('~/panalysis/output/%s_thesis/test/%s_regions.root' % (mu.name,mu.name))
    rf = ROOT.TFile('%s/%s_regions.root' % (OUTPUT_DIR, mu.name))
    region = 'SR_OS'
    k_entry = Drawer.upsilon_CR
    entry = Drawer.upsilon
    verbose = False
    fmt = 'latex'
    kW_OS = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose)
    kW_OS_L = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, groups = mu.MC+[mu.LH])
    kW_OS_R = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, groups = mu.MC+[mu.RH])
    pdg_simple.testpdg(kW_OS.val,kW_OS.err)
    print 'kWOS %g +/ %g' % ( pdg_simple.pdg_round(kW_OS.val,kW_OS.err))

    EST = [ErrorFloat(0),ErrorFloat(0),ErrorFloat(0)]
    BG = [ErrorFloat(0),ErrorFloat(0),ErrorFloat(0)]
    for group in [mu.data,mu.Zlljet,mu.Zltau,mu.ttbar]:
        h = rf.Get('h_%s_%s_%s' % (entry.name,group.name,region))

        nG = ErrorFloat(*hist.IntegralError(h))
        if group.isData:
            print row( group.latex, [int(nG.val),'',''], fmt=fmt)
        else:
            print row( group.latex, ['%i +/- %i'%(nG.val,nG.err),'',''], fmt=fmt, makeint=True)
        if not group.isData:
            for i in xrange(len(EST)):
                EST[i]+=nG
                BG[i]+=nG
                
    print row('', ['SM','Left-Handed','Right-Handed'], fmt = fmt)
    signal = []
    SB = [ErrorFloat(0),ErrorFloat(0),ErrorFloat(0)]
    for i,group in enumerate([mu.signal, mu.LH, mu.RH]):
        h = rf.Get('h_%s_%s_%s' % (entry.name,group.name,region))
        nG = ErrorFloat(*hist.IntegralError(h))
        signal.append('%i +/- %i'%(nG.val,nG.err))
        EST[i]+=nG
        SB[i]=nG
    print row(mu.signal.latex, signal, fmt=fmt, makeint = True)
    wjets = []
    for i,group in enumerate([mu.signal, mu.LH, mu.RH]):
        h =  hist.get_dataW(rf, mu, Drawer.upsilon, 'WCR_OS', verbose = False, groups = mu.MC + [group])
        print '**'
        nG = ErrorFloat(*hist.IntegralError(h))
        wjets.append('%g +/- %g'%( pdg_simple.pdg_round((nG).val,(nG).err)))
        EST[i]+=nG
        BG[i]+=nG

    print row(mu.Wlnu.latex, wjets, fmt=fmt)
        #print row(group.latex, [nG])
    QCD = []
    for i,group in enumerate([mu.signal, mu.LH, mu.RH]):
        kW = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='OS', verbose = verbose, groups = mu.MC+[group])
        kWSS = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='SS', verbose = verbose, groups = mu.MC+[group])
        rqcd = hist.calc_rQCD(rf,mu,k_entry,'QCD',kW_OS = kW, kW_SS = kWSS, verbose = verbose, groups = mu.MC+[group])
        print '%g +/- %g' % pdg_simple.pdg_round(rqcd.val,rqcd.err)
        qcd = hist.get_QCD(rf,mu,entry,'SR_SS', kW = kWSS, rqcd = rqcd, groups = mu.MC +[group])

        nG = ErrorFloat(*hist.IntegralError(qcd))
        QCD.append('%g +/- %g'% pdg_simple.pdg_round(nG.val,nG.err))
        EST[i]+=nG
        BG[i]+=nG        
    print row('Multijet', QCD, fmt=fmt)

    for i in xrange(len(EST)):
        EST[i] = '%i +/- %i'%(EST[i].val,EST[i].err)
        SB[i] = SB[i]/BG[i]

    print row('Total Estimate', EST, fmt = fmt)
    print row('S:B', SB, fmt = fmt)

def stat_error(mu):
    rf = ROOT.TFile('%s/%s_regions.root' % (OUTPUT_DIR,mu.name))
    k_entry = Drawer.upsilon_28
    wSR_MC = rf.Get('h_%s_%s_%s' % (k_entry.name,mu.Wlnu.name,'SR_OS'))
    nMC = ErrorFloat(*hist.IntegralError(wSR_MC))
    # print ErrorFloat( wSR_MC.GetBinContent(10), wSR_MC.GetBinError(10)  )
    kW_OS = hist.calc_kX(rf,mu,k_entry,'WCR','Wlnu',OS='OS')
    # print kW_OS * nMC
    # wSR_MC.Scale(kW_OS.val)
    for bin in xrange(wSR_MC.GetNbinsX()+1):
        bi0 = ErrorFloat(wSR_MC.GetBinContent(bin), wSR_MC.GetBinError(bin))
        bi1 = bi0 * ErrorFloat(kW_OS.val,0)
        bi2 = bi0 * kW_OS
        # print bi0,bi1,bi2
        wSR_MC.SetBinContent(bin,bi2.val)
        wSR_MC.SetBinError(bin,bi2.err)

    nMCkw = ErrorFloat(*hist.IntegralError(wSR_MC))
    print nMCkw
    #print ErrorFloat( wSR_MC.GetBinContent(10), wSR_MC.GetBinError(10)  )

def stat_error_wcr(mu):
    rf = ROOT.TFile('%s/%s_regions.root' % (OUTPUT_DIR,mu.name))
    k_entry = Drawer.upsilon_CR
    hist.get_dataW(rf,mu,k_entry,'WCR_OS')

def modes(mu, hand):
    print '====== %s ========' % mu.name
    rf = ROOT.TFile('%s/%s_regions.root' % (OUTPUT_DIR,mu.name))
    region = 'SR_OS'
    modes = ['NP','PINU','SPX','RHO']
    total = rf.Get('h_upsilon_%s_Ztautau_%s' % (hand,region)).Clone()
    nTotal = total.Integral()
    print 'TOTAL \t %i' % nTotal
    nMode = 0

    for m in modes:
        h = rf.Get('h_upsilon_%s_%s_%s' % (hand,m,region)).Clone()
        nm = h.Integral()
        print '%s %s \t %i %.3f' % (hand, m, nm, nm/nTotal)
        nMode+=nm
    print 'Mode Sum \t %i' % nMode


def rW(mu):
    rf = ROOT.TFile('%s/%s_regions.root' % (OUTPUT_DIR, mu.name))
    entry = Drawer.upsilon_CR
    for region in ['SR_OS','SR_SS','QCD_OS','QCD_SS']:
        if 'SR' in region:
            entry = Drawer.upsilon_CR
            sentry = Drawer.upsilon
        else:
            entry = Drawer.upsilon_CR
            sentry = Drawer.upsilon_CR
        rw = hist.calc_rW(rf,mu,entry,region, sentry = sentry )
        print region, '%g $\pm$ %g' % pdg_simple.pdg_round( rw.val, rw.err)


def sim_modes(sim="Reco"): 
    rf = ROOT.TFile('~/analysis/truth/cxx/thesis/PythiaFull30.root')
    entries = ['upsilon','upsilon20','upsilon40','upsilon60']
    modes = ['0','3','4','5']
    fracs = []
    for i,mode in enumerate(modes):
        fracs.append([])
        for j,entry in enumerate(entries):
            fracs[i].append( ErrorFloat(*hist.IntegralError(rf.Get('%s_%s_%s_LH' % (sim,entry,mode)))))
            fracs[i].append(ErrorFloat(*hist.IntegralError(rf.Get('%s_%s_%s_RH' % (sim,entry,mode)))))

    nTotal = list(fracs[0])
    nOther = None

    for i,m in enumerate(modes):
        if i == 1: nOther = list(fracs[0])
        for j in xrange(len(fracs[i])):
            fracs[i][j] = fracs[i][j]/nTotal[j]
            fracs[i][j] = round(fracs[i][j].val,2)
            if nOther:
                nOther[j]-=fracs[i][j]

        print row( modelegend[m], fracs[i], fmt='latex' , flt = '%.2f')
    print row('Other',nOther, fmt = 'latex', flt = '%.2f')

        

if __name__ == '__main__':
    signal_region(streams.mu)

