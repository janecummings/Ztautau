import sys
import ROOT
from printing import row
import streams
import run
import hist
from ErrorFloat import ErrorFloat
import Drawer
fcol=20;col=20

def twikirow(cut,entries):
    twow = '| %s |' % cut
    for entry in entries:
        twow += ' %s |' % entry
    return twow


def print_region(rf,regions,groups,entry=Drawer.upsilon,eff=False,twiki=False, fmt = None, f = False):
    a = []
    for i,region in enumerate(regions):
        cutgroup = []
        n = []
        mc = 0
        for group in groups:
            h = rf.Get('h_%s_%s_%s' % (entry.name,group.name,region.name))
            cutgroup.append( '%i +/- %i' % hist.IntegralError(h))
            n.append(h.Integral(0,h.GetNbinsX()+1))
            if group.name not in ['Muons','Egamma','LH_Ztautau','RH_Ztautau']: 
                #print group.name
                mc += h.Integral(0,h.GetNbinsX()+1)
            if i == 0: a = list(n); l=a
        if twiki: print twikirow(region.name,cutgroup)
        else: print row(region.name,cutgroup,fcol=fcol,col=col, fmt = fmt)
        if f:
            print row(region.name, [''] + [x/mc for x in n[1:]],fcol=fcol,col=col,fmt=fmt)
        if eff:
            print row('',[x/y for (x,y) in zip(n,a)],fcol=fcol,col=col, fmt = fmt)
            print row('',[x/y for (x,y) in zip(n,l)],fcol=fcol,col=col, fmt = fmt)
        l = n
def print_regions(s,twiki=False, fmt = None, f = False):
    stream = s.stream
    if twiki: print twikirow('Region/Sample',[group.name for group in stream.groups])
    else: print row('Region/Sample',[group.name for group in stream.groups],fcol=fcol,col=col, fmt = fmt)
    sf = ROOT.TFile('%s/%s_regions.root'%(s.OUTPUTDIR,s.stream.name))
    if not fmt: print '-'*(fcol+col*len(stream.groups))
    print_region(sf,stream.SR,stream.groups, fmt = fmt, f = f)
    if not s.control: return 0
    for regions in s.control:
        if not fmt: print '-'*(fcol+col*len(stream.groups))
        print_region(sf,regions,stream.groups,twiki=twiki, fmt = fmt, f = f )

def print_selection(s,twiki=False, fmt = None):
    stream = s.stream
    if twiki: print twikirow('Region/Sample',[group.name for group in stream.groups])
    else: print row('Region/Sample',[group.name for group in stream.groups],fcol=fcol,col=col, fmt = fmt)
    sf = ROOT.TFile('%s/%s_regions.root'%(s.OUTPUTDIR,s.stream.name))
    if not twiki: print '-'*(fcol+col*len(s.groups))
    print_region(sf,stream.select,s.groups,eff = False,twiki=twiki, fmt = fmt)


def print_all_selection(s,twiki=False, fmt = None):
    stream = s.stream
    if twiki: print twikirow('Region/Sample',[group.name for group in stream.groups])
    else: print row('Region/Sample',[group.name for group in stream.groups],fcol=fcol,col=col, fmt = fmt)
    sf = ROOT.TFile('%s/%s_regions.root'%(s.OUTPUTDIR,s.stream.name))
    if not fmt: print '-'*(fcol+col*len(s.groups))
    print_region(sf,stream.select,s.groups,entry = Drawer.upsilon, eff = False,twiki=twiki, f = False, fmt = fmt)

def print_region_definitions(s):
    for region in s.control + [s.stream.select]:
        for r in region:
            print r.string()
            



if __name__ == '__main__':

    if len(sys.argv)>1:
        s = getattr(run,sys.argv[1])
        #print_regions(s)
        #print_selection(s,twiki=True)
        #print_regions(s, fmt = 'twiki')
        #print_all_selection(s, fmt = 'twiki')
        #print_region_definitions(s)
