import ROOT
import sys

rf = ROOT.TFile( sys.argv[1] )
tau = rf.Get('tau')
keys = tau.GetListOfBranches()
print len(keys)
names = []
for key in keys:
    if not sys.argv[1] in key.GetName(): continue
    if key.GetName() in names: continue
    print key.GetName()
    names.append(key.GetName())

    c = ROOT.TCanvas()
    print tau.Draw('%s>>hist0(100,-2,2)'%(key.GetName()),"tau_mode_RHO")
    h = ROOT.gDirectory.Get("hist0")

    h.Draw()
    c.SaveAs("rootfile.gif+1")
