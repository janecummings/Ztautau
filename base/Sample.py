import glob,ROOT

PATH_TO_MC = '/data/cummings/common/prod4_mc/'  
PATH_TO_DATA = '/data/cummings/common/prod4_data/' 
PATH_TO_TRUE = '/data/cummings/common/truth/'    
PATH_TO_THRY = '/data/cummings/common/prod4_theory/'    

class Sample:
    def __init__(self,
                 name = '',
                 files=[],
                 xsec = 0.0,
                 xpol = None,
                 int_lumi = 0.0,
                 weight = 'evtsel_weight/evtsel_weight_FF',
                 DSID = None,
                 version = None,
                 stream = None,
                 string = None,
                 xstr = None,
                 theory = False):
        
        self.name = name
        if files == [] and 'period' not in str(DSID):
            version = version or ''
            if not theory:
                files = glob.glob(PATH_TO_MC+"user*%s*%s*/*root*"%(DSID,version))
            else:
                files = glob.glob(PATH_TO_THRY+"user*%s*%s*/*root*"%(DSID,version))

        elif files==[] and DSID and stream:
            if not version: version = 'grp14'
            files = glob.glob(PATH_TO_DATA+"user*%s*%s*PhysCont.%s*/*root*"%(DSID,stream,version))
        self.files = files
        self.xsec = xsec
        self.xpol = xpol
        self.xstr = xstr
        self.isData = bool(stream)
        self.stream = stream
        self.int_lumi = int_lumi
        self.DSID = DSID
        self.version = version
        self.weight = weight
        #self.ldsn = files[0].split('/')[-1].split('
        if self.isData: self.weight = '1.0'

    def load(self,chain='tau'):
        ch = ROOT.TChain(chain)
        self.totalEvents = 0
        self.root_files = []
        for f in self.files:
            try:
                rf = ROOT.TFile(f,"READ")
                self.root_files.append(rf)
                self.totalEvents += rf.Get('TotalEvents').GetBinContent(1)
            except:
                self.totalEvents += 0
            ch.Add(f)
        self.chain = ch
        self.nEvents = ch.GetEntries()

    def close(self):
        for rf in self.root_files:
            rf.Close()

class TruthSample:
    def __init__(self, name = '', PATH = PATH_TO_TRUE,  ls = ''):
        self.name = name
        self.files = glob.glob(PATH_TO_TRUE+"user.cummings*%s*/*root*" %(ls or name))

    def load(self,chain="tau"):
        ch = ROOT.TChain(chain)
        for f in self.files:
            ch.Add(f)
        self.chain = ch
                
class Group(list):
    def __init__(self,name,samples,string=None,weight = None, cut = None, fill = '', kfactor = 1.0 , legend = None, isData = False, latex = None):
        self.name = name
        self.legend = legend or name
        self.latex = latex or name
        self.extend(samples)
        self.string = string
        self.cut = cut
        self.kfactor = kfactor
        self.isData = isData
        if not self.string and self.cut:
            self.string = self.cut.string()

        self.weight = weight
        self.fill = fill
        self.xpol = None
    def __add__(self,other):
        x = Group(self.name,self,string=self.string,cut=self.cut)
        x.extend(other)
        return x
    def clone(self,name,cut = None, string = None, weight = None, fill = '', kfactor = None):
        cut = cut or self.cut
        string = string or self.string
        weight = weight or self.weight
        kfactor = kfactor or self.kfactor 
        fill = fill or self.fill
        return self.__class__(name,self,cut = cut, string = string, weight = weight, fill = fill, kfactor = kfactor)

class Cut(list):
    def __init__(self, variable, comparator = '', value = 0, fn=None, var = None, name = ''):
        self.extend([variable,comparator,value])
        self.fn = fn
        self.var = var
        self.name = name

    def string(self):
        if self[0] == '':
            return '1.0'
        
        if self[1]:
            if self[1] == '==' and self[2] == False:
                return '(!%s)' % self[0]
            else:
                if self.fn: 
                    if self.fn == 'abs': return '(%s(%s)%s%s)' %(self.fn,self[0],self[1],self[2])
                    if self.fn == 'div': return '((%s/%s)%s%s)' %(self[0],self.var,self[1],self[2])
                return '(%s%s%s)' % (self[0], self[1], self[2])

        else:
            return '(%s)' % self[0]
        
    def evaluate(self,event):
        if self[0] == '': return True
        event_value = getattr(event, self[0].lstrip('!'))
        if self.fn: 
            if self.fn == 'abs': event_value = abs(event_value)
            if self.fn == 'div':
                if type(self.var) == str:
                    div_value = getattr(event, self.var)
                else: 
                    div_val = self.var
                    event_value = event_value/div_value
        comparator = self[1]
        compared_to = self[2]
         ## smaller than
        if comparator == '<':
            if not event_value < compared_to: return False
             
            ## smaller than or equal
        elif comparator == '<=':
            if not event_value <= compared_to: return False

            ## greater than
        elif comparator == '>':
            if not event_value > compared_to: return False

            ## greater than or equal
        elif comparator == '>=':
            if not event_value >= compared_to: return False

            ## Checking boolean to be true
        elif comparator == '':
            if self[0][0] == '!':
                if event_value: return False
            else:
                if not event_value: return False
                
            ## equal
        elif comparator == '==':
            if not event_value == compared_to: return False

        else:
            raise ValueError('Comparator %s unknown, add it into the evaluate method of the CutList class.' % comparator)

        return True

class CutSet(list):
    def __init__(self,cuts,logic='and',name = ''):
        self.extend(cuts)
        self.logic = logic
        self.name = name

    def string(self):
        if self.logic == 'or': lo = '||'
        elif self.logic == 'and': lo = '&&'
        return '(%s)'%lo.join([cut.string() for cut in self])

    def evaluate(self,event):
        cut_booleans = [cut.evaluate(event) for cut in self]
        answer = False
        if self.logic == 'and':
            answer = True
            for cut_boolean in cut_booleans:
                answer = answer and cut_boolean
        if self.logic == 'or':
            answer = False
            for cut_boolean in cut_booleans:
                answer = answer or cut_boolean      
        return answer

class Region:
    def __init__(self,name,cutset = [],ln=''):
        self.name = name
        self.cutset = cutset
        if ln == '': self.ln = name
        else: self.ln = ln

    def string(self):
        return ' && '.join([cut.string() for cut in self.cutset])

    def evaluate(self,event):
        cut_booleans = [cut.evaluate(event) for cut in self.cutset]
        #print cut_booleans
        answer = True
        for cut_boolean in cut_booleans:
            answer = answer and cut_boolean
        return answer

