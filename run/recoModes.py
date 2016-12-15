import os,sys,time
import ROOT
import run
import Groups
import trueTau
from Sample import Cut
import streams
from printing import row
from Drawer import *

fcol=30;col=15
int_lumi = 20276.9

#############################
# Modes #####################
#############################
p0n = Cut('mode','==','1p0n', name = '1p0n')
p1n = Cut('mode','==','1p1n', name = '1p1n')
pXn = CutSet([Cut('nP','==',1),Cut('nN','>',1)],name = '1pXn')
p = Cut('nP','==',1, name = '1p')
np = Cut('nP','>',1, name = 'Xp')


def makeModes(stream,group,regions,OUTPUTDIR,name='regions'):

    modes = [p0n,p1n,pXn,np]
    hists = {};clones={}

    if not os.path.exists(OUTPUTDIR): os.mkdir(OUTPUTDIR)
    out = ROOT.TFile('%s/%s_%s.root'%(OUTPUTDIR,stream.name,name),'UPDATE')
    for m in modes:
        hists[m.name]={};clones[m.name]={}
        for region in regions:
            hists[m.name][region.name]={};clones[m.name][region.name]={}
            for entry in entries:
                hname = 'h_%s_%s_%s_%s' % (entry.name,group.name,m.name,region.name)
                th = ROOT.TH1F(hname,hname,entry.nbins,entry.binlo,entry.binhi)
                hists[m.name][region.name][entry.name] = th
                clones[m.name][region.name][entry.name] = {}
                for sample in group:
                    clones[m.name][region.name][entry.name][sample.name] = th.Clone(sample.name)
    n_events = {}
    n_good = {}
    for sample in group:
        sample.load()
        ch = sample.chain
        
        n = 0 
        sr_os = regions[0].name; sr_ss = regions[1].name
        n_events[sample.name] = {sr_os:0,sr_ss:0}
        n_good[sample.name] = {sr_os:0,sr_ss:0}
        #print sample.name
        #print sample.totalEvents
        #print sample.nEvents

        for ientry in xrange(ch.GetEntries()):
            jentry = ch.LoadTree(ientry)
            if jentry < 0: continue
            nb = ch.GetEntry(ientry)
            if nb <= 0: continue
            if group.cut and not group.cut.evaluate(ch): continue ## this cuts on left and right handed samples in first production and real taus in ttbar
            n=n+1            
            

            for region in regions:
                if not region.evaluate(ch): continue
                n_events[sample.name][region.name] +=1
                out.cd()
                tau = trueTau.recoTau(ch)
                tau.dmode()
                if not tau.calculate(): 
                    continue
                n_good[sample.name][region.name]+=1
            
                if stream.weight:
                    weight = getattr(ch,stream.weight)
                else:
                    weight = ch.evtsel_weight 
                for m in modes:
                    if not m.evaluate(tau): continue
                    for entry in entries:
                        dsth = clones[m.name][region.name][entry.name][sample.name]
                        if entry.name in ['upsilon','tau_leadTrkPt']:
                            dsth.Fill(getattr(tau,entry.name),weight)
                        else:
                            dsth.Fill(getattr(ch,entry.variable),weight)
        #logging.info(row(sample.name,[sample.totalEvents,sample.nEvents,n],fcol=fcol,col=col))
        print '--------'
        print n_good
        print n_events

    for m in modes:
        for region in regions:
            for entry in entries:
                for sample in group:
                    scale = sample.xsec * int_lumi * n_events[sample.name][region.name] / (sample.totalEvents * n_good[sample.name][region.name])
                    hists[m.name][region.name][entry.name].Add(clones[m.name][region.name][entry.name][sample.name],scale)
                hists[m.name][region.name][entry.name].Write()


    del clones
    out.Write()
    del hists
    


if __name__ == '__main__':
    
    spec = run.el_tight_plots
    stream = spec.stream
    groups = [stream.LH, stream.RH]

    OUTPUTDIR = spec.OUTPUTDIR
    entries = spec.reg_entries
    regions = stream.SR
    for group in groups:
        makeModes(stream,group,regions, OUTPUTDIR)


