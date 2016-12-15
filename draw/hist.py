import sys
import ROOT
from ErrorFloat import ErrorFloat
from printing import row
fcol = 25;col=25

# scale LH and RH to Ztautau cross section
kLH = ErrorFloat(1.,0.)#ErrorFloat(x=1/(1-0.425),dx=0) # first production
kRH = ErrorFloat(1.,0.)#ErrorFloat(x=1/0.425,dx=0) # first production

def IntegralError(hist, lo = None, hi = None, of = True):
    error = ROOT.Double(0)
    if not lo: lo = 0
    if not hi:
        if of: 
            hi = hist.GetNbinsX()+1   
        else:
            hi = hist.GetNbinsX() 
    integral = hist.IntegralAndError(lo,hi,error)
    error = float(error)
    return integral,error

def calc_kX(rf,stream,entry,region,skipon,OS='OS',kZ=ErrorFloat(1.0,0.0),kT=ErrorFloat(1.0,0.0),kW=ErrorFloat(1.0,0.0),groups=[],cr='', ptau = None,
            kZL = kLH, kZR = kRH,verbose=False, sys = None, sf = None, ud = None):

    if cr: region = region + '_' + cr
    if not isinstance(skipon,list): skipon = [skipon]
    hname = 'h_%s_%s_%s_%s'% (entry.name,stream.data.name,region,OS)
    #print hname
    data = rf.Get(hname)
    nData = ErrorFloat(*IntegralError(data))

    printK = '(nData ' 
    what = '(%s' % nData
    
    if sys and sf: rf = sf
    nW = ErrorFloat(0.0,0.0)
    if groups == []:
        groups = stream.MC + [stream.signal]
    if ptau: 
        groups = stream.MC + [stream.LH,stream.RH]
        ptau = ErrorFloat(ptau,0)

    for group in groups:
        hname = 'h_%s_%s_%s_%s' % (entry.name,group.name,region,OS)
        if sys: hname += '_%s_%s' % (sys,ud)
        hBG = rf.Get(hname).Clone()
        nBG = ErrorFloat(*IntegralError(hBG))
        if group.name in skipon:
            nW += nBG
            continue
        if ptau and 'LH' in group.name: nBG = nBG*ErrorFloat(0.5*(1-ptau),0)
        if ptau and 'RH' in group.name: nBG = nBG*ErrorFloat(0.5*(1+ptau),0)
        if group.name == 'Zlljet': nBG = nBG*kZ 
        if group.name == 'ttBar': nBG = nBG*kT 
        if group.name == 'Wlnu' : nBG = nBG*kW
        printK += '- N_%s' % group.name
        what += '- (%s)' % nBG
        nData = nData - nBG
    printK += ') / N_%s' % skipon[0]
    what += ') / (%s)' % nW
    # print calculation
    if verbose:
        print 'k_%s_%s = %s = \n %s = %s' % (region,OS,printK,what,nData/nW)

    return nData/nW

def calc_rQCD(rf,stream,entry,region,kW_OS=ErrorFloat(1.,0.),kW_SS=ErrorFloat(1.,0.),groups = [], kZL = kLH, kZR = kRH, verbose=False, sys = None, sf = None, ud = None):
    #print 'h_%s_%s_%s_OS' % (entry.name,stream.data.name,region)
    qcdOS = rf.Get('h_%s_%s_%s_OS' % (entry.name,stream.data.name,region)).Clone()
    qcdSS = rf.Get('h_%s_%s_%s_SS' % (entry.name,stream.data.name,region)).Clone()

    num_name = 'OS Data '
    den_name = 'SS Data '
    num = '%s ' % ErrorFloat(*IntegralError(qcdOS))
    den = '%s ' % ErrorFloat(*IntegralError(qcdSS))

    if sys and not sf: sf = rf

    if groups == []:
        groups = stream.MC + [stream.signal]

    for group in groups:

        num_name += '- OS %s' % group.name
        den_name += '- SS %s' % group.name
        if not sys:
            hOS = rf.Get('h_%s_%s_%s_OS' % (entry.name,group.name,region)).Clone()
            hSS = rf.Get('h_%s_%s_%s_SS' % (entry.name,group.name,region)).Clone()
        else:
            hOS = sf.Get('h_%s_%s_%s_OS_%s_%s' % (entry.name,group.name,region,sys,ud)).Clone()
            hSS = sf.Get('h_%s_%s_%s_SS_%s_%s' % (entry.name,group.name,region,sys,ud)).Clone()

        fOS = ErrorFloat(-1.0,0.)
        fSS = ErrorFloat(-1.0,0.)

        if group.name == 'Ztautau_LH': fOS = fOS * kZL; fSS = fSS * kZL
        if group.name == 'Ztautau_RH': fOS = fOS * kZR; fSS = fSS * kZR

        if group.name == 'Wlnu': fOS = kW_OS*fOS; fSS = kW_SS*fSS
        
        qcdOS.Add(hOS,fOS.val)
        qcdSS.Add(hSS,fSS.val)
        
        num += ' %s' % (ErrorFloat(*(IntegralError(hOS)))*fOS)
        den += ' %s' % (ErrorFloat(*(IntegralError(hSS)))*fSS)

    nOS = ErrorFloat(*IntegralError(qcdOS))
    nSS = ErrorFloat(*IntegralError(qcdSS))

    if verbose:
        print num_name + ' /// ' + den_name
        print num + ' /// ' + den
        print nOS/nSS

    return nOS/nSS

def get_QCD(rf,stream,entry,region, kW=ErrorFloat(1.,0.), rqcd = ErrorFloat(1.,0.), groups = [], kZL = kLH, kZR = kRH, verbose=False, ptau = None, sysf = None):
    #print 'h_%s_%s_%s' % (entry.name,stream.data.name, region)
    qcd = rf.Get('h_%s_%s_%s' % (entry.name,stream.data.name, region)).Clone()
    if type(kW) is not 'ErrorFloat': kW = ErrorFloat(kW,0.)
    if type(rqcd) is not 'ErrorFloat': rqcd = ErrorFloat(rqcd,0)
    if type(kZL) is not 'ErrorFloat': kZL = ErrorFloat(kZL,0)
    if type(kZR) is not 'ErrorFloat': kZR = ErrorFloat(kZR,0)
    if groups == []:
        groups = stream.MC + [stream.signal]
    if ptau:
        groups = stream.MC + [stream.LH,stream.RH]
    if sysf: rf = sysf

    for group in groups:
        h = rf.Get('h_%s_%s_%s' % (entry.name,group.name,region)).Clone()

        fS = -1.
        if group.name == 'Wlnu': fS = fS * kW.val
        if ptau and 'LH' in group.name: fS = 0.5 * (1 - ptau)
        if ptau and 'RH' in group.name: fS = 0.5 * (1 + ptau)
        #print group.name, h.Integral()
        qcd.Add(h,fS)

    qcd.Scale(rqcd.val)
    if verbose:
        print 'QCD Integral: %s' % ErrorFloat(*IntegralError(qcd))
    return qcd

def calc_rW(rf,stream,entry,region, sys = '', sentry = None ):
    # W MC 
    if sentry == None:
        sentry = entry
    #print sentry.name
    if 'OS' in region: wreg ='WCR_OS'
    else: wreg = 'WCR_SS'
    if not sys:
        #print 'h_%s_%s_%s' % (sentry.name,stream.Wlnu.name, region)
        wmc = rf.Get('h_%s_%s_%s' % (entry.name,stream.Wlnu.name, wreg)).Clone()
        srw = rf.Get('h_%s_%s_%s' % (sentry.name,stream.Wlnu.name, region)).Clone()
    else:
        #print 'h_%s_%s_%s_%s' % (sentry.name,stream.Wlnu.name, region,sys)
        wmc = rf.Get('h_%s_%s_%s_%s' % (entry.name,stream.Wlnu.name, wreg, sys)).Clone()
        srw = rf.Get('h_%s_%s_%s_%s' % (sentry.name,stream.Wlnu.name, region,sys)).Clone()

    #print 'h_%s_%s_%s' % (entry.name,stream.Wlnu.name, wreg)
    #v,e = IntegralError(wmc, of = True)
    #print v,e
    nSRW = ErrorFloat(*IntegralError(srw))
    nWMC = ErrorFloat(*IntegralError(wmc))
    #print nSRW, nWMC
    rW = nSRW/nWMC    

    return rW

def get_nWCR( rf, stream, entry, region, groups = [], sys = ''):
    rf.cd()
    wcr = rf.Get('h_%s_%s_%s' % (entry.name,stream.data.name, region)).Clone()
    nData = ErrorFloat(*IntegralError(wcr))
    if groups == []:
        groups = stream.MC + [stream.signal]
    for group in groups:
        if not sys:
            h = rf.Get('h_%s_%s_%s' % (entry.name,group.name,region)).Clone()
        else:
            print 'h_%s_%s_%s_%s' % (entry.name,group.name,region,sys)
            h = rf.Get('h_%s_%s_%s_%s' % (entry.name,group.name,region,sys)).Clone()
        fS = -1.
        # skip over W MC
        if group.name == 'Wlnu': 
            nWMC = ErrorFloat(*IntegralError(h))
            continue
        # subtract MC EW from estimate
        wcr.Add(h,fS)
    nEST = ErrorFloat(*IntegralError(wcr))
    return wcr

def get_dataW(rf,stream,entry,region, groups = [], verbose=False, sf = None, sysf='', where = 'SR'):
    # data in WCR
    rf.cd()
    #print 'h_%s_%s_%s' % (entry.name,stream.data.name, region)
    wcr = rf.Get('h_%s_%s_%s' % (entry.name,stream.data.name, region)).Clone()
    nData = ErrorFloat(*IntegralError(wcr))
    #print wcr
    if groups == []:
        groups = stream.MC + [stream.signal]

    for group in groups:
        #print group.name
        if not sysf:
            rf.cd()
            #print 'h_%s_%s_%s' % (entry.name,group.name,region)
            h = rf.Get('h_%s_%s_%s' % (entry.name,group.name,region)).Clone()
        else:
            sf.cd()
            #print 'h_%s_%s_%s_%s' % (entry.name,group.name,region,sysf)
            h = sf.Get('h_%s_%s_%s_%s' % (entry.name,group.name,region,sysf)).Clone()
        fS = -1.
        # skip over W MC
        if group.name == 'Wlnu': 
            nWMC = ErrorFloat(*IntegralError(h))
            continue
        # subtract MC EW from estimate
        wcr.Add(h,fS)
    nEST = ErrorFloat(*IntegralError(wcr))

    kW = nEST/nWMC
    SR = region.replace('WCR',where)
    if not sysf:
        rf.cd()
        srw = rf.Get('h_%s_%s_%s' % (entry.name,stream.Wlnu.name, SR)).Clone()
    else:
        sf.cd()
        srw = sf.Get('h_%s_%s_%s_%s' % (entry.name,stream.Wlnu.name, SR, sysf)).Clone()
    nSRW = ErrorFloat(*IntegralError(srw))
    rW = nSRW/nWMC
    #print nSRW, nWMC, rW
    if verbose:
        print kW
        print rW
        print ErrorFloat(*IntegralError(wcr))

    for bin in xrange(wcr.GetNbinsX()+1):
        bc = ErrorFloat( wcr.GetBinContent(bin), wcr.GetBinError(bin))
        bc = bc * rW
        wcr.SetBinContent(bin,bc.val)
        wcr.SetBinError(bin,bc.err)
    #print ErrorFloat(wcr.GetBinContent(10),wcr.GetBinError(10))
    #print ErrorFloat(*IntegralError(wcr)) 
    if verbose:
        #print wcr.GetBinContent(1),wcr.GetBinError(1)
        print ErrorFloat(*IntegralError(wcr))    
        #print wk.GetBinContent(1), wk.GetBinError(1)

    return wcr


def get_nW(rf,stream, entry, region, kentry = None, groups = [], sf = None, sysf = '', verbose=False):
    # data in WCR
    if not kentry:
        kentry = entry
    wcr = rf.Get('h_%s_%s_%s' % (kentry.name,stream.data.name, region)).Clone()
    nData = ErrorFloat(*IntegralError(wcr))
    

    if groups == []:
        groups = stream.MC + [stream.signal]

    for group in groups:
        if not sysf:
            rf.cd()
            h = rf.Get('h_%s_%s_%s' % (kentry.name,group.name,region)).Clone()
        else:
            sf.cd()
            h = sf.Get('h_%s_%s_%s_%s' % (kentry.name,group.name,region,sysf)).Clone()
        fS = -1.
        # skip over W MC
        if group.name == 'Wlnu': 
            nWMC = ErrorFloat(*IntegralError(h))
            continue
        # subtract MC EW from estimate
        wcr.Add(h,fS)
    nEST = ErrorFloat(*IntegralError(wcr))

    kW = nEST/nWMC
    SR = region.replace('WCR','SR')
    if not sysf:
        rf.cd()
        srw = rf.Get('h_%s_%s_%s' % (entry.name,stream.Wlnu.name, SR)).Clone()
    else:
        sf.cd()
        srw = sf.Get('h_%s_%s_%s_%s' % (entry.name,stream.Wlnu.name, SR, sysf)).Clone()
    nSRW = ErrorFloat(*IntegralError(srw))
    rW = nSRW/nWMC
    if verbose:
        print kW
        print rW
        print ErrorFloat(*IntegralError(wcr))
    for bin in xrange(wcr.GetNbinsX()+1):
        bc = ErrorFloat( wcr.GetBinContent(bin), wcr.GetBinError(bin))
        bc = bc * rW
        wcr.SetBinContent(bin,bc.val)
        wcr.SetBinError(bin,bc.err)
    #print ErrorFloat(wcr.GetBinContent(10),wcr.GetBinError(10))
    #print ErrorFloat(*IntegralError(wcr)) 
    if verbose:
        #print wcr.GetBinContent(1),wcr.GetBinError(1)
        print ErrorFloat(*IntegralError(wcr))    
        #print wk.GetBinContent(1), wk.GetBinError(1)

    return wcr.Integral()


def get_QCDW(rf,stream,entry,region, rqcd = ErrorFloat(1.,0.), groups = [], kZL = kLH, kZR = kRH, verbose=False, ptau = None, sysf = None):
    #print 'h_%s_%s_%s' % (entry.name,stream.data.name, region)
    qcd = rf.Get('h_%s_%s_%s' % (entry.name,stream.data.name, region)).Clone()
    if groups == []:
        groups = stream.MC + [stream.signal]
    if ptau:
        groups = stream.MC + [stream.LH,stream.RH]
    if sysf: rf = sysf

    for group in groups:
        h = rf.Get('h_%s_%s_%s' % (entry.name,group.name,region)).Clone()
        fS = -1.
        if group.name == 'Wlnu': 
            #ssw = get_dataW(rf,stream,entry,region.replace('SR','WCR'),groups = groups,verbose=verbose)
            ssw = rf.Get('h_%s_%s_%s' % (entry.name,'WData',region)).Clone()
            qcd.Add(ssw,fS)
            continue
        qcd.Add(h,fS)
    qcd.Scale(rqcd.val)
    if verbose:
        print 'QCD Integral: %s' % ErrorFloat(*IntegralError(qcd))

    return qcd

def get_SS_EW(rf,stream,entry,region,kW=ErrorFloat(1.,0.), rqcd = ErrorFloat(1.,0.), groups = [], kZL = kLH, kZR = kRH, verbose=False):
    if type(kW) is not 'ErrorFloat': kW = ErrorFloat(kW,0.)
    ew = rf.Get('h_%s_%s_%s' % (entry.name,stream.data.name,region)).Clone()
    ew.Reset()

    if groups == []:
        groups = stream.MC + [stream.signal]

    for group in groups:

        h = rf.Get('h_%s_%s_%s' % (entry.name,group.name,region)).Clone()
        fS = 1.
        if group.name == 'Wlnu': fS = fS * kW.val
        
        ew.Add(h,fS)

    return ew


def get_OS(rf,entry,group,region,kOS=ErrorFloat(1.0,0.0)):

    if kOS is None: kOS = ErrorFloat(1.,0.) 
    elif not type(kOS) == ErrorFloat: kOS = ErrorFloat(kOS,0.)
    OS = rf.Get('h_%s_%s_%s' % (entry.name,group.name,region)).Clone()
    OS.Scale(kOS.val)
    return OS

def get_SS(rf,entry,group,region,kSS=ErrorFloat(1.,0.)):
    region = region.replace('_OS','_SS')
    SS = get_OS(rf,entry,group,region,kOS=kSS)
    return SS

def get_OS_SS(rf,entry,group,region,kOS=ErrorFloat(1.0,0.0),kSS=ErrorFloat(1.0,0.0),rQCD=None,verbose=False):
    # you can include rQCD in kSS 
    if rQCD: kSS = kSS * rQCD
    OS = get_OS(rf,entry,group,region,kOS=kOS)
    SS = get_SS(rf,entry,group,region,kSS=kSS)
    os_ = '%.0f +/- %.0f' % (IntegralError(OS)[0],IntegralError(OS)[1])
    ss_ = '%.0f +/- %.0f' % (IntegralError(SS)[0],IntegralError(SS)[1])  
    OS.Add(SS,-1.0)
    os_ss_ = '%.0f +/- %.0f' % (IntegralError(OS)[0],IntegralError(OS)[1])
    if verbose:
        print row(group.name,[os_,ss_,os_ss_],fcol=fcol,col=col)
    return OS


def get_sys_error(stack, up, down):
    errors = ROOT.TGraphAsymmErrors(stack)
    for b in xrange(stack.GetNbinsX()):
        b_up = up.GetBinContent(b+1)
        b_down = down.GetBinContent(b+1)
        b_nom = stack.GetBinContent(b+1)

        hi =0; lo =0
        if b_up > b_nom : 
            hi = b_up - b_nom
        if b_down < b_nom:
            lo = b_nom - b_down

        hi = max(0, hi, b_down - b_nom)
        lo = max(0, lo, b_nom - b_up)

        errors.SetPointEYlow(b+1,lo)
        errors.SetPointEYhigh(b+1,hi)


    return errors




