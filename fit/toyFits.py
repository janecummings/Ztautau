import ROOT
import sys

# name of HistFactory model
name = 'FitTest_ehad_noQCDfit_NoZpt'

# change path to HistFactory Workspace
rf = ROOT.TFile('TestResults/%s_combined_%s_model.root' % (name,name))

w = rf.Get('combined')
mc = w.obj('ModelConfig')
pois = mc.GetParametersOfInterest()
data = w.data("obsData")
poi = pois.first()
pdf = mc.GetPdf()
obs = mc.GetObservables()
nuis = mc.GetNuisanceParameters()
globs = mc.GetGlobalObservables()

nuisAndGlobsAndPOI = ROOT.RooArgSet(pois)
nuisAndGlobsAndPOI.add(nuis)
nuisAndGlobsAndPOI.add(globs)

# Use an initial fit to the expected data to initialize fit parameters for generating toys 
def fit_expected():
    # expected data stored in "obsData"
    result = pdf.fitTo(data,ROOT.RooFit.Extended(),
                        ROOT.RooCmdArg(ROOT.RooFit.Save(ROOT.kTRUE)),
                        ROOT.RooCmdArg(ROOT.RooFit.Minos()), ) 

    result.Print("v")
    print result.edm()
    w.saveSnapshot("nuisAndGlobsAndPOI",nuisAndGlobsAndPOI)

# generate nToys and save FitResult objects in a root file
def fit_toys( nToys ):

    fit_expected()
    
    pc = ROOT.RooStats.ProofConfig(w, 4, "workers=4", False)
    sampler = ROOT.RooStats.ToyMCSampler()
    sampler.SetPdf( pdf ) 
    sampler.SetObservables(mc.GetObservables())
    sampler.SetNuisanceParameters(mc.GetNuisanceParameters())
    sampler.SetGlobalObservables(mc.GetGlobalObservables())
    if not mc.GetPdf().canBeExtended() and (data.numEntries()==1):
         sampler.SetNEventsPerToy(1)
    sampler.SetParametersForTestStat(mc.GetParametersOfInterest())
    sampler.SetProofConfig(pc)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

    f = ROOT.TFile('FitResults.root','UPDATE')
    for seed in xrange(nToys):
        ROOT.RooRandom.randomGenerator().SetSeed(seed)
        w.loadSnapshot("nuisAndGlobsAndPOI")
        toyData = sampler.GenerateToyData(nuisAndGlobsAndPOI)
        result = pdf.fitTo(toyData,ROOT.RooCmdArg(ROOT.RooFit.Extended()),                                  # add extended term in likelihood to fit to Nevents
                                     ROOT.RooCmdArg(ROOT.RooCmdArg(ROOT.RooFit.PrintLevel(-1))),            # turn off messages
                                     ROOT.RooCmdArg(ROOT.RooFit.Minos()),                                   # turn on MINOS (default migrad+hesse)
                                     ROOT.RooCmdArg(ROOT.RooFit.Save(ROOT.kTRUE)) )                         # return FitResult object

        f.cd()
        result.Write( 'fr_%s'% seed )
    f.Write()
    f.Close()

# you can develop code to read the fit results from the file / plot distributions 
# here is an example of printing the initial and final values of the fit in a single fit result
def fit_result( seed ):
    f = ROOT.TFile('FitResults.root')
    fr = f.Get('fr_%i' % seed )

    initPars = fr.floatParsInit()
    initPars.Print("v")
    
    finalPars = fr.floatParsFinal()
    finalPars.Print("v")
    f.Close()

if __name__ == '__main__':

    ###############################################
    fit_expected()
    ###############################################
    nToys = 1
    #fit_toys(nToys)
    #fit_result(0)