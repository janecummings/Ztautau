import sys
import ROOT

f = ROOT.TFile(sys.argv[1])
if len(sys.argv)>2: 
    t = f.Get(sys.argv[2])
    if type(t) is ROOT.TTree:
        t.Print('toponly')
    elif type(t) is ROOT.TDirectoryFile:
        t.ls()
else:
    f.ls()
f.Close()
