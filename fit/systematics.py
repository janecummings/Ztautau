import ROOT
import sys
import Groups, hist, run, Drawer
from printing import row
import math
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.LoadMacro("../AtlasStyle.C")
ROOT.gROOT.LoadMacro("../AtlasLabels.C")
ROOT.SetAtlasStyle()
ROOT.gStyle.SetFrameBorderSize(1)

# 
spec = run.el_syst
output_dir = spec.OUTPUTDIR
mu = spec.stream

tesf = ROOT.TFile('%s/%s_regions.root' % (output_dir,mu.name)) # tes_final shape & nominal
regf = ROOT.TFile('%s/%s_regions.root' % (output_dir,mu.name))
sysf = ROOT.TFile('%s/%s_regions.root' % (output_dir,mu.name))
sff =ROOT.TFile('%s/%s_syst.root' % (output_dir,mu.name))

entry = spec.reg_entries[0]

# Systematic Uncertainty Categories
syst = ['ElES_LowPt','ElES_PS','ElES_R12','ElES_Zee','ElEnResSys','JER','JES_BJet','JES_Detector1','JES_EtaMethod',
                     'JES_EtaModelling','JES_FlavComp','JES_FlavResp','JES_Modelling1','JES_PUMu','JES_PUNPV','JES_PURho',
                     'JES_Statistical1','JVF','METResSys','METScaleSys','MuSys']
sysSF = ['el_id','el_iso','mu_id','mu_iso','tau_el','tau_fake','PU_rescaling','tau_id_stat','tau_id_sys']
tes = ['TES_FINAL']
tausys = ['tau_id_sys','tau_id_stat', 'tau_el', 'tau_fake']
lepsys = ['el_id','el_iso','mu_id','mu_iso']
lepES = ['ElES_LowPt','ElES_PS','ElES_R12','ElES_Zee','ElEnResSys', 'MuSys' ]
JES = ['JER','JES_BJet','JES_Detector1','JES_EtaMethod',
                     'JES_EtaModelling','JES_FlavComp','JES_FlavResp','JES_Modelling1','JES_PUMu','JES_PUNPV','JES_PURho',
                     'JES_Statistical1','JVF','METResSys','METScaleSys']
sysSF = ['TES_FINAL']

groups = mu.MC + [mu.LH,mu.RH]

def KS_test(sys,region):
    rf = tesf
    sf = tesf 

    gsys = {}

    for group in groups:
        rf.cd()
        print group.name
        if group.name == 'Wlnu':
            if 'OS' in region: wregion = 'WCR_OS'
            elif 'SS' in region: wregion = 'WCR_SS'
            nom = hist.get_dataW(rf, mu, entry, wregion, verbose = False)
        else:
            nom = rf.Get('h_upsilon_%s_%s' % (group.name,region))
        #print group.name, nom.Integral()
        #continue
        for i,s in enumerate(sys):
            #if 'PU' in s: continue
            for ud in ['UP','DOWN']:
                #ud = ud.lower()
                if not s+'_'+ud in gsys:
                    gsys[s+'_'+ud] = []
                if not nom.Integral():
                    gsys[s+'_'+ud].append('-')
                    continue
                if group.name == 'Wlnu':
                    sf.cd()
                    shift = hist.get_dataW(rf, mu, entry, wregion, verbose = False, sf = sf, sysf = s+'_'+ud)
                else:
                    sf.cd()
                    shift = sf.Get('h_upsilon_%s_%s_%s_%s' % (group.name,region,s,ud ))
                print shift.Integral()
                gsys[s+'_'+ud].append(nom.KolmogorovTest(shift))


    print row('',[group.name for group in groups], fcol =15, col=15, fmt = 'latex')

    for s in sys:
        #if 'PU' in s: continue
        for ud in ['UP','DOWN']:
            #ud = ud.lower()
            what = s+'_'+ud
            for i,k in enumerate(gsys[what]):
                if k<0.9: 
                    gsys[what][i] = '\\textcolor{red}{%.4f}' % k
                    whatp = '\\textcolor{red}{%s}' % what
                else:
                    whatp = what
            print row(whatp, gsys[what],fcol=15,col=10, fmt = 'latex')



def NormStat(sys,region = 'SR_OS'):
    gsys = {}

    for group in groups:
        rf.cd()
        if group.name == 'Wlnu':
            if 'OS' in region: wregion = 'WCR_OS'
            elif 'SS' in region: wregion = 'WCR_SS'
            nom = hist.get_dataW(rf, mu, entry, wregion, verbose = False)
        else:
            nom = rf.Get('h_upsilon_%s_%s' % (group.name,region))
        #print group.name,nom.Integral()
        #continue
        for i,s in enumerate(syst):
            if 'PU' in s: continue
            for ud in ['UP','DOWN']:
                ud = ud.lower()
                if not s+'_'+ud in gsys:
                    gsys[s+'_'+ud] = []
                if not nom.Integral():
                    gsys[s+'_'+ud].append('-')
                    continue;
                if group.name == 'Wlnu':
                    sf.cd()
                    shift = hist.get_dataW(rf, mu, entry, wregion, verbose = False, sf = sf, sysf = s+'_'+ud)
                else:
                    sf.cd()
                    shift = sf.Get('h_upsilon_%s_%s_%s_%s' % (group.name,region,s,ud ))
                sigmastat = abs(shift.Integral() - nom.Integral())/math.sqrt(nom.Integral())
                gsys[s+'_'+ud].append(sigmastat)


    print row('',[group.name for group in groups], fcol =15, col=15, fmt = 'latex')

    for s in syst:
        if 'PU' in s: continue
        for ud in ['UP','DOWN']:
            ud = ud.lower()
            what = s+'_'+ud
            whatp = what
            imax = -1
            kmax = 0
            for i,k in enumerate(gsys[what]):
                if k>0.01: 
                    gsys[what][i] = '\\textcolor{red}{%.4f}' % k
                    whatp = '\\textcolor{red}{%s}' % what
                    if k>kmax:
                        kmax = k
                        imax = i
                if k <= 0.: gsys[what][i] = '-'
            if imax>=0:
                gsys[what][imax]=gsys[what][imax].replace('red','blue')

            print row(whatp, gsys[what],fcol=15,col=10, fmt = 'latex')


def NormStatFull(s,region):
    rf = regf
    sf = tesf
    sigma_up = []
    sigma_down = []
    intervals = []
    entry = Drawer.upsilon_20
    if 'SR' not in region:
        entry = Drawer.upsilon_CR
    for group in groups:
        rf.cd()
        nom = rf.Get('h_%s_%s_%s' % (entry.name,group.name,region))
        nom = nom.Integral()
        if not nom:
            sigma_up.append('-')
            sigma_down.append('-')
            intervals.append('-')
            continue;
        for ud in ['UP','DOWN']:
            #ud = ud.lower()

            #if group.name == 'Wlnu':
            #   sf.cd()
            #   shift = hist.get_dataW(rf, mu, entry, wregion, verbose = False, sf = sf, sysf = s+'_'+ud)
            #else:
            sf.cd()
            shift = sf.Get('h_%s_%s_%s_%s_%s' % (entry.name,group.name,region,s,ud ))
            sigmastat = abs(shift.Integral() - nom)/math.sqrt(nom)


            if ud == 'UP': 
                sigma_up.append(sigmastat)
                shift_up = shift.Integral()
            else:
                sigma_down.append(sigmastat)
                shift_down = shift.Integral()
        s_up = shift_up/nom
        s_down = shift_down/nom
        if (s_up<1 and s_down<1) or (s_up>1 and s_down>1):
            intervals.append('\\textcolor{red}{%.3f,%.3f}' % (s_down,s_up))
        else:
            intervals.append('%.3f,%.3f' % (s_down,s_up))
    print row(region, ['' for group in groups], fmt = 'latex')
    print row('$\\sigma_up$', sigma_up, fmt='latex'  )
    print row('$\\sigma_down$', sigma_down, fmt='latex'  )
    print row( 'down,up', intervals, fmt = 'latex')
    print '\n\\hline'


def draw_diff( s, region ):
    rf = regf
    sf = tesf
    entry = Drawer.upsilon_20
    for group in groups:
        rf.cd()
        ratio = True
        c = ROOT.TCanvas(s+region+group.name,s+region+group.name,1500,1000)
        #Split for ratio plot
        c.Divide(1,2)
        c.cd(1).SetPad(0.0,0.2,1.0,1.0)
        c.cd(2).SetPad(0.0,0.0,1.0,0.2)
        c.cd(1)
        l = ROOT.TLegend(.75,.75,.9,.9)

        l.SetTextSize(0.035)
        l.SetBorderSize(0)
        l.SetFillColor(0)

        nom = rf.Get('h_%s_%s_%s' % (entry.name,group.name,region)).Clone()
        nom_stat = nom.Clone()
        nom_stat.SetFillStyle(3354)
    
        nom_stat.SetFillColor(12)
        l.AddEntry(nom,group.name)

        if not nom.Integral():
            continue;
        up = sf.Get('h_%s_%s_%s_%s_%s' % (entry.name,group.name,region,s,'UP' ))
        down = sf.Get('h_%s_%s_%s_%s_%s' % (entry.name,group.name,region,s,'DOWN' ))

        up.SetLineStyle(2)
        down.SetLineStyle(4)
        up.SetLineColor(6)
        down.SetLineColor(8)

        nom.SetMaximum(1.4*nom.GetMaximum())
        l.AddEntry(up, '%s Up' % s)
        l.AddEntry(down, '%s Down' % s)

        for h in [nom,up,down,nom_stat]:
            h.SetMarkerSize(0)
            h.SetLineWidth(1)
        for h in [up,down]:
            h.SetLineWidth(2)

        nom.Draw("HIST")
        nom_stat.Draw("SAME,E2")
        up.Draw("SAME,HIST")
        down.Draw("SAME,HIST")
        l.Draw()
        tm = ROOT.TPaveText(0.19,0.83,0.22,0.90,"NDC")
        tm.SetBorderSize(0)
        tm.SetFillColor(0)
        tm.SetTextSize(0.05)
        tm.AddText(region)
        tm.Draw()

        if ratio:
            c.cd(2)
            ROOT.gPad.SetTopMargin(0)
            dataratio = up.Clone('ratio')
            downratio = down.Clone('down')

            dataratio.GetYaxis().SetTitle(' SHIFT / NOM')
            dataratio.GetYaxis().SetTitleOffset(0.25)
            dataratio.GetYaxis().SetTitleSize(0.16)
            dataratio.SetLabelSize(.1,"y")
            dataratio.Divide(nom)
            dataratio.SetMaximum(1.2)
            dataratio.SetMinimum(0.8)

            downratio.Divide(nom)
            
            err = nom_stat.Clone('err')
            err.Divide(err)
            err.SetFillStyle(3354)
            err.SetLineColor(1)
            err.SetFillColor(12)
            err.SetMarkerStyle(1)
            xmin = dataratio.GetXaxis().GetXmin()
            xmax = dataratio.GetXaxis().GetXmax()
            cl = ROOT.TLine(xmin,1,xmax,1)
            cl.SetLineStyle(3)
            dataratio.GetXaxis().SetLabelSize(0.04)
            dataratio.Draw('PE')
            downratio.Draw('PE,SAME')
            cl.Draw("SAME")
            err.Draw("SAME,E2")
            c.Update()

        c.SaveAs('../plots/systematics/%s_lepES.gif+1' % mu.name)

def print_latex_full():
    print '\\documentclass[11pt,oneside,a4paper]{article}'
    print '\\usepackage{rotating}'
    print '\\usepackage[usenames]{color}'
    print '\\begin{document}'

    for sys in lepsys + tausys + lepES + JES:
        print '\\begin{table}'
        print '\\centering'
        print '\\caption{\\bf{%s}}' % sys.replace('_','\_')
        print '\\begin{tabular}{lrrrrrr}'
        print row('',[group.name for group in groups], fmt = 'latex')


        for region in ['SR_OS','SR_SS','WCR_OS','WCR_SS','QCD_OS','QCD_SS']:
            NormStatFull( sys, region)
        print '\\end{tabular}'
        print '\\end{table}'

    print '\\end{document}'

for sys in lepES:
    for region in ['SR_OS','SR_SS']:
        draw_diff(sys, region )