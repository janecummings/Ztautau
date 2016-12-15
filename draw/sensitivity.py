import ROOT
import streams
import Drawer
import hist
from ErrorFloat import ErrorFloat
import math
from printing import row


Recofile = 'TestXX'


def sensitivity(rh,lh,p):
    c = ROOT.TCanvas()

    lh.Scale(1.0/(lh.Integral())) #(0,lh.GetNbinsX()+1)))
    rh.Scale(1.0/(rh.Integral())) #(0,rh.GetNbinsX()+1)))

    summ = lh.Clone()
    diff = lh.Clone()
    multi = lh.Clone()
    ratio = lh.Clone()
    sump = lh.Clone()

    summ.Add(lh,rh,1.0,1.0)
    diff.Add(rh,lh,1.0,-1.0)

    scale = 1/(summ.Integral())#(0,summ.GetNbinsX()+1 ))
    summ.Scale(scale)
    diff.Scale(scale)

    summ.SetLineColor(2)
    diff.SetLineColor(4)
    summ.SetAxisRange(-1,1)
    diff.SetAxisRange(-1,1)
    summ.SetMinimum(1.2*diff.GetMinimum())

    summ.Draw()

    diff.Draw("SAME")

    #diff.Draw()
    multi.Multiply(diff,diff)
    sump.Add(summ,diff,1,p)
    ratio.Divide(multi,sump)

    c.SaveAs('sensTest.gif+1')

    s_square = ErrorFloat(*hist.IntegralError(ratio))
    s = math.sqrt(s_square.val)
    s_err = 0.5 * s_square.err / s_square.val
    return ErrorFloat(round(s,2),round(s_err,2))


def get_sens(mu, p):
    rf = ROOT.TFile('~/panalysis/output/%s_thesis/%s_regions.root' % (mu.name,mu.name))

    entry = Drawer.upsilon_28
    lh = rf.Get('h_%s_%s_SR_OS' % (entry.name, mu.LH.name  )).Clone("lh")
    rh = rf.Get('h_%s_%s_SR_OS' % (entry.name, mu.RH.name  )).Clone("lh")

    return sensitivity(rh,lh,p)

def sens_reco(p,mode, entry = 'upsilon'):
    rf = ROOT.TFile('~/analysis/truth/cxx/thesis/Pythia%s.root' % Recofile)
    lh = rf.Get('Reco_%s_%s_LH' % (entry,mode)).Clone('lh')
    rh = rf.Get('Reco_%s_%s_RH' % (entry,mode)).Clone('rh')

    return sensitivity(rh,lh,p)

def sens_true(p, mode, entry = 'upsilon'):
    rf = ROOT.TFile('~/analysis/truth/cxx/thesis/Pythia%s.root' % Recofile)
    lh = rf.Get('True_%s_%s_LH' % (entry,mode)).Clone('lh')
    rh = rf.Get('True_%s_%s_RH' % (entry,mode)).Clone('rh')
    return sensitivity(rh,lh,p)

def reco():
    mode = '0'
    for entry in ['','20','40','60']:
        entry = 'upsilon'+entry
        print row(entry,[sens_reco(-1,mode,entry = entry), sens_reco(0,mode,entry = entry), sens_reco(1,mode,entry = entry)], fmt = 'latex')
def true():
    mode = '4'
    print row('',[sens_true(-1,mode), sens_true(0,mode), sens_true(1,mode)], fmt = 'latex')
    entry = 'upsilon20'
    print row(entry,[sens_true(-1,mode,entry=entry), sens_true(0,mode,entry=entry), sens_true(1,mode,entry=entry)], fmt = 'latex')
    entry = 'upsilon60'
    print row('upsilon40',[sens_true(-1,mode,entry=entry), sens_true(0,mode,entry=entry), sens_true(1,mode,entry=entry)], fmt = 'latex')
    entry = 'upsilon40'
    print row('upsilon60',[sens_true(-1,mode,entry=entry), sens_true(0,mode,entry=entry), sens_true(1,mode,entry=entry)], fmt = 'latex')
def xfrac():
    mode = '3'
    entry = 'costh'
    print row(entry,[sens_true(-1,mode,entry=entry), sens_true(0,mode,entry=entry), sens_true(1,mode,entry=entry)], fmt = 'latex')
    entry = 'xfrac'
    print row(entry,[sens_true(-1,mode,entry=entry), sens_true(0,mode,entry=entry), sens_true(1,mode,entry=entry)], fmt = 'latex')
    mode = '4'
    entry = 'costh'
    print row(entry+',rho',[sens_true(-1,mode,entry=entry), sens_true(0,mode,entry=entry), sens_true(1,mode,entry=entry)], fmt = 'latex')




