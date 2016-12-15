import ROOT
import math

pdg = {16: "nu",111: "pi0",211: "pi",
       311: "K0",321: "K",310: "K0S",130: "K0L"
       }
neutral = [111,310,311,130]
charged = [321,211]

class trueTau():
    def __init__(self,event):
        self.pt = getattr(event,'evtsel_truethad_pt')
        self.eta = getattr(event,'evtsel_truethad_eta')
        self.phi = getattr(event,'evtsel_truethad_phi')
        self.dau_pdgId = list(getattr(event,'evtsel_truethad_dau_pdgId'))
        self.dau_pt = list(getattr(event,'evtsel_truethad_dau_pt'))
        self.dau_eta = list(getattr(event,'evtsel_truethad_dau_eta'))
        self.dau_phi = list(getattr(event,'evtsel_truethad_dau_phi'))
        
        self.mode = ''
        self.nP = -1
        self.nN = -1

    def dmode(self):  
        m = [abs(x) for x in self.dau_pdgId]
        name = ''
        nc = 0
        nn = 0
        for p in [111,211,311,321,310,130]:
            if m.count(p):
                if p in charged: nc+=m.count(p)
                if p in neutral: nn+=m.count(p)
                name = '%ip%in' % (nc,nn)
        self.nP = nc
        self.nN = nn
        self.mode = name
        return name

    def calculate(self):
        m = [abs(x) for x in self.dau_pdgId]
        if m == []: 
            return False
        self.charged = 0.0
        self.neutral = 0.0
        self.neutrino = ROOT.TLorentzVector()
        vistau = ROOT.TLorentzVector()
        for i,p in enumerate(m):
            if p not in [16,111,211,311,321,310,130]: continue
            vt = ROOT.TLorentzVector()
            vt.SetPtEtaPhiM(
                self.dau_pt[i],
                self.dau_eta[i],
                self.dau_phi[i],
                0.177682)
            vistau += vt
            if p == 16: self.neutrino.SetPtEtaPhiM(
                self.dau_pt[i],
                self.dau_eta[i],
                self.dau_phi[i],
                0.0)
            else:
                if p in charged: 
                    self.charged += self.dau_pt[i]
                elif p in neutral: 
                    self.neutral += self.dau_pt[i]
        self.vistau = vistau - self.neutrino
        self.vispt = self.vistau.Pt()
        self.upsilon = (self.charged - self.neutral) / (self.charged + self.neutral)
        return True


class recoTau():
    def __init__(self,event):
         self.dau_pdgId = list(getattr(event,'evtsel_truethad_dau_pdgId'))
         self.tau_leadTrkPt = list(getattr(event,'tau_leadTrkPt'))[0]/1000.
         self.upsilon = -1. + 2 * self.tau_leadTrkPt / getattr(event,'evtsel_tau_et')
         self.mode = ''
         self.nP = -1
         self.nN = -1

    def dmode(self):  
        m = [abs(x) for x in self.dau_pdgId]
        name = ''
        nc = 0
        nn = 0
        for p in [111,211,311,321,310,130]:
            if m.count(p):
                if p in charged: nc+=m.count(p)
                if p in neutral: nn+=m.count(p)
                name = '%ip%in' % (nc,nn)
        self.nP = nc
        self.nN = nn
        self.mode = name
        return name

    def calculate(self):
        m = [abs(x) for x in self.dau_pdgId]
        if m == []: 
            return False
        return True


class trueLep():
    def __init__(self,event):
        self.dau_pdgId = list(getattr(event,'evtsel_truetlep_dau_pdgId'))
        self.dau_pt = list(getattr(event,'evtsel_truetlep_dau_pt'))
        self.dau_eta = list(getattr(event,'evtsel_truetlep_dau_eta'))
        self.dau_phi = list(getattr(event,'evtsel_truetlep_dau_phi'))
        
    def calculate(self):
        m = [abs(x) for x in self.dau_pdgId]
        if m == []: return False
        self.neutrino = ROOT.TLorentzVector()
        self.lep = ROOT.TLorentzVector()
        if not m.count(11) and not m.count(13): return False

        for i,p in enumerate(m):
            if p == 11:
                self.flavor = 'el'
                self.pt = self.dau_pt[i]
                self.phi = self.dau_phi[i]
                self.lep.SetPtEtaPhiM(self.dau_pt[i],
                                      self.dau_eta[i],
                                      self.dau_phi[i],
                                      0.000510999)
            elif p == 13:
                self.flavor = 'mu'
                self.pt = self.dau_pt[i]
                self.phi = self.dau_phi[i]
                self.lep.SetPtEtaPhiM(self.dau_pt[i],
                                      self.dau_eta[i],
                                      self.dau_phi[i],
                                      0.105658367)
            if p in [12,14,16]: 
                nu = ROOT.TLorentzVector()
                nu.SetPtEtaPhiM(
                    self.dau_pt[i],
                    self.dau_eta[i],
                    self.dau_phi[i],
                    0.0)
                self.neutrino+=nu
        return True

class trueEvent():
    def __init__(self, tau_vispt, lep_pt, transverseMass, CosDeltaPhi,leadPt):
        self.tau_vispt = tau_vispt
        self.lep_pt = lep_pt
        self.transverseMass = transverseMass
        self.CosDeltaPhi = CosDeltaPhi
        self.leadPt = leadPt

def transverseMass(tau,lep):
    met = tau.neutrino + lep.neutrino
    cos = 1. - math.cos(lep.lep.DeltaPhi(met))
    return math.sqrt(2.0*lep.lep.Pt()*met.Et()*cos)

def deltaPhi(tau,lep):
    met = tau.neutrino + lep.neutrino
    sumdphi = (tau.phi - met.Phi()) + (lep.phi - met.Phi())
    return sumdphi
    #return math.cos(lep.phi - tau.phi)
    #return math.cos(lep.lep.DeltaPhi(tau.vistau))
    
