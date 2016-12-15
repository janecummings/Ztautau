from Drawer import *
import streams,Groups

ANALYSIS_DIR = '/afs/cern.ch/user/c/cummings/panalysis'

    #####################
    ### run specs #######
    #####################

class spec:
    def __init__(self,
                 stream,
                 DIR='test',
                 loglady = True,
                 comment= '',
                 regions = [],
                 control = [],
                 groups = [],
                 reg_entries = [],
                 cr_entries = [],
                 unscaled = False, # for the fit,
                 cutflow = None,
                 syst = None,
                 SFsys = None,
                 update = False,
                 filename = None
                 ):
        self.stream = stream
        self.loglady = loglady
        self.DIR = DIR
        self.OUTPUTDIR = '%s/output/%s' % (ANALYSIS_DIR,DIR)
        self.plotdir = '%s/plots/%s' % (ANALYSIS_DIR,DIR)
        self.comment = comment
        self.cutflow = cutflow
        self.syst = syst
        self.SFsys = SFsys
        self.runSys = bool(syst or SFsys)
        self.update = update
        self.filename = filename
        #files (list of regions for each output file to produce)
        if regions is None: self.regions = False
        elif regions == []:
            # default is the full selection 
            # for only the final SR (OS and SS) selection do regions = ['SR']
            self.regions = [stream.selectOS,stream.selectSS] 
        else: self.regions = [getattr(stream,r) for r in regions]
        if control is None: self.control = False
        elif control == []: 
            # default is every control region in (WCR,WTCR,ZLCR,TCR,QCD)
            # for only WCR: control = ['WCR']
            self.control = [stream.ZLCR,
                            stream.WTCR,
                            stream.WCR,
                            stream.TCR,
                            stream.QCD,
                            ]
        else: self.control = [getattr(stream,c) for c in control]
        #format groups  // hmm should i just leave this to stream? 
        if groups == []:
            if stream.groups: self.groups = stream.groups
            else:
                self.groups = []
                for group in [stream.data, stream.signal,stream.LH,stream.RH] + stream.MC:
                    if group: self.groups.append(group)
        else: self.groups = groups#[getattr(stream,g) for g in groups]
        #entries
        if regions and reg_entries == []:
            reg_entries = [tau_et]
        self.reg_entries = reg_entries
        if control and cr_entries == []:
            cr_entries = [tau_et]
        self.cr_entries = cr_entries
        self.unscaled = unscaled

