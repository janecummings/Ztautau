import os,sys,time
import ROOT
import streams
from Drawer import *
import logging, logLady
import run
import json
from printing import row

fcol = 20;col = 15
int_lumi = 20243.6

def print_info(s):
    stream = s.stream
    OUTPUTDIR = s.OUTPUTDIR
    ## Stream Info    
    print '\nstream: %s' % stream.name

    if not s.filename:
        outputname = '%s_%s' % (stream.name, 'regions')
    else:
        outputname = s.filename
    ## Output dir
    print 'OUTPUTDIR = %s/%s.root\n' % (OUTPUTDIR,outputname)

    # Overall stream weight
    if stream.weight: print 'WEIGHT: %s' % stream.weight

    print '*'*100     
    
    if s.regions: 
        time.sleep(3)
        print 'Regions:'
        for region in s.regions:
            print '[%s]:' % ', '.join([r.name for r in region])
            for r in region:
                print '%s: %s' % (r.name,r.string())
        print '\nplotting variables: \n'
        for entry in s.reg_entries:
            print '%s \t %s %s %s' % (entry.name, entry.nbins, entry.binlo, entry.binhi)
    if s.control:
        print '*'*100 
        time.sleep(3)
        print 'Control regions:'
        for region in s.control:
            print '[%s]:' % ', '.join([r.name for r in region])
            for r in region:
                print '%s: %s' % (r.name,r.string())
        #print '\nplotting variables: %s\n' % ', '.join([entry.name for entry in s.cr_entries])
        print '\nplotting variables: \n'
        for entry in s.cr_entries:
            print '%s \t %s %s %s' % (entry.name, entry.nbins, entry.binlo, entry.binhi)

    print '*'*100 
    time.sleep(3)

    print 'Groups: '
    for group in s.groups:
        dsid = [str(ds.DSID) for ds in group]
        print '\n%s : %s' % (group.name,', '.join(dsid))
        print '  %s kfactor: %.2f' % (group.name,group.kfactor)
        if group.weight:
            print '  %s weight: %s' % (group.name,group.weight)
        if group.string:
            print '  %s string: %s' % (group.name,group.string)

    print '*'*100
    print 'unscaled = %s' % s.unscaled
    print 'Tree Systematics: ', s.syst
    print 'Branch Systematics: ', s.SFsys


def makeOut(s):

    print_info(s)
    print 'Preparing to run ...'
    time.sleep(10)

    stream = s.stream
    OUTPUTDIR = s.OUTPUTDIR

    if s.loglady:
        logLady.makeLog(OUTPUTDIR)
        if s.regions: 
            logging.info('Regions: ')
            for region in s.regions:
                logging.info(' %s' % ', '.join([r.name for r in region]))
            logging.info('plotting variables: %s' % ', '.join([entry.name for entry in s.reg_entries]))
        if s.control: 
            logging.info('Control regions: ')
            for region in s.control:
                logging.info(' %s' % ', '.join([r.name for r in region]))
            logging.info('CR variables: %s' % ', '.join([entry.name for entry in s.cr_entries]))
        logging.info('Groups: ')
        for group in s.groups:
            dsid = [str(ds.DSID) for ds in group]
            logging.info('%s : %s' % (group.name,', '.join(dsid)))
            if group.weight:
                logging.info('  %s weight: %s' % (group.name,group.weight))
            if group.string:
                logging.info('  %s string: %s' % (group.name,group.string))

    if not os.path.exists(OUTPUTDIR): os.mkdir(OUTPUTDIR)

    ## Run nominal tree, no systematics 
    if not s.runSys:
        if s.regions:
            makeFiles(s,'regions',s.regions,s.reg_entries,s.groups,s.loglady,s.unscaled)
        if s.control:
            makeFiles(s,'regions',s.control,s.cr_entries,s.groups,s.loglady,False)
    ## Tree Systematics 
    if s.syst:
        for sys in s.syst:
            sysUP = 'SystematicsUP/' + sys
            sysDOWN = 'SystematicsDOWN/' + sys
            if s.regions:
                makeFiles(s,'regions',s.regions,s.reg_entries,s.groups,s.loglady,s.unscaled,syst=sysUP)
                makeFiles(s,'regions',s.regions,s.reg_entries,s.groups,s.loglady,s.unscaled,syst=sysDOWN)
            if s.control:
                makeFiles(s,'regions',s.control,s.cr_entries,s.groups,s.loglady,s.unscaled,syst=sysUP)
                makeFiles(s,'regions',s.control,s.cr_entries,s.groups,s.loglady,s.unscaled,syst=sysDOWN)
    ## Scale Factor systematics 
    if s.SFsys:
        for sf in s.SFsys:
            if s.regions:
                makeFiles(s,'regions',s.regions,s.reg_entries,s.groups,s.loglady,s.unscaled,SFsys=sf)
            if s.control:
                makeFiles(s,'regions',s.control,s.cr_entries,s.groups,s.loglady,s.unscaled,SFsys=sf)

def makeFiles(s,name,regions,entries,groups,loglady,unscaled,syst=None, SFsys = None):
    n=0   
    stream = s.stream
    OUTPUTDIR = s.OUTPUTDIR
    if syst:
        sys = syst.split('/')[1]
        up = syst.split('/')[0].lstrip('Systematics')
        #name = sys
    if SFsys:
        if 'PU_' in SFsys:
            sf = 'evtsel_sys_' + SFsys
        else:
            sf = 'evtsel_sys_sf_' + SFsys
        sf_up = sf + '_up'
        sf_down = sf + '_down'
    

    SM = False

    output = '%s/' % OUTPUTDIR
    if not s.filename: 
        output += '%s_%s' % (stream.name,name)
    else: output += s.filename

    out = ROOT.TFile('%s.root' % output, 'UPDATE') #'RECREATE')
    if s.loglady: 
        logging.info('%s ... writing to %s.root' % (logLady.atime(),output))

    for i,entry in enumerate(entries):
        th2 = False
        if isinstance(entry,plot2d): th2 = True
        if loglady: logging.info('... plotting %s' % entry.name)

        for region in regions: 
            for r in region:
                ### turn this on for running systematics in SR only
                #if (SFsys or syst) and 'SS' in r.name: 
                #    continue

                reg = r.string()
                if loglady and i == 0: 
                    logging.info('%s: %s' % (r.name,r.string()))
                ###test region strings###
                #print r.string()
                #continue
                #########################
                gth = None
                th = None

                for group in groups:
                    if (syst or SFsys) and group.isData: continue
                    if syst:
                        gth_name = 'h_%s_%s_%s_%s_%s' % (entry.name, group.name, r.name, sys, up)
                    elif SFsys:
                        SFsys_name_up = SFsys+"_UP"
                        SFsys_name_down = SFsys+"_DOWN"
                        gth_name_up = 'h_%s_%s_%s_%s' % (entry.name, group.name, r.name, SFsys_name_up )
                        gth_name_down = 'h_%s_%s_%s_%s' % (entry.name, group.name, r.name, SFsys_name_down)
                    else:
                        gth_name = 'h_%s_%s_%s' % (entry.name, group.name, r.name)
                    if th2:
                        gth = ROOT.TH2F(gth_name,gth_name,entry.nxbins,entry.xbinlo,entry.xbinhi,entry.nybins,entry.ybinlo,entry.ybinhi)
                    elif SFsys:
                        gth_up = ROOT.TH1F(gth_name_up,gth_name_up,entry.nbins,entry.binlo,entry.binhi)
                        gth_down = ROOT.TH1F(gth_name_down,gth_name_down,entry.nbins,entry.binlo,entry.binhi)
                        gth_up.Sumw2()
                        gth_down.Sumw2()                        
                    else:
                        gth = ROOT.TH1F(gth_name,gth_name,entry.nbins,entry.binlo,entry.binhi)
                        if not group.isData: gth.Sumw2()     

                    for ds in group:
                        if group.string is not None:
                            draw = reg + ' && ' + group.string
                        else: draw = reg

                        if syst: ds.load(chain=syst)
                        else: ds.load()
                        ch = ds.chain

                        out.cd()

                        if th2:
                            th_name = gth_name + '_' + ds.name + '_unscaled'                            
                            th = ROOT.TH2F(th_name,th_name,entry.nxbins,entry.xbinlo,entry.xbinhi,entry.nybins,entry.ybinlo,entry.ybinhi)                        
                        elif SFsys:
                            th_name_up = gth_name_up + '_' + ds.name + '_unscaled'
                            th_name_down = gth_name_down + '_' + ds.name + '_unscaled'
                            th_up = ROOT.TH1F(th_name_up,th_name_up,entry.nbins,entry.binlo,entry.binhi)                            
                            th_down = ROOT.TH1F(th_name_down,th_name_down,entry.nbins,entry.binlo,entry.binhi)                            
                        else:
                            th_name = gth_name + '_' + ds.name + '_unscaled'                            
                            th = ROOT.TH1F(th_name,th_name,entry.nbins,entry.binlo,entry.binhi)
                        if not ds.isData: 
                            if th:
                                th.Sumw2()
                            elif SFsys:
                                th_up.Sumw2()
                                th_down.Sumw2()

                        if not ds.isData and stream.weight:
                            weight = stream.weight
                        else:
                            weight = ds.weight

                        if group.weight is not None:
                            weight = weight + '*' + group.weight

                        if not ds.isData and SFsys: 
                            weight_up = weight + '*' + sf_up
                            weight_down = weight + '*' + sf_down
                            weightdraw_up = '%s * ( %s )' % (weight_up,draw)
                            weightdraw_down = '%s * ( %s )' % (weight_down,draw)
                        else:
                            weightdraw = '%s * ( %s )' % (weight,draw)

                        if th2:
                            #print '(%s):(%s)>>%s'%(entry.xvariable,entry.yvariable,th_name)
                            ch.Draw('(%s):(%s)>>%s'%(entry.yvariable,entry.xvariable,th_name),weightdraw)
                        elif SFsys:
                            ch.Draw('%s>>%s'%(entry.variable,th_name_up),weightdraw_up)  * ds.xsec * int_lumi / ds.totalEvents
                            ch.Draw('%s>>%s'%(entry.variable,th_name_down),weightdraw_down)  * ds.xsec * int_lumi / ds.totalEvents
                        else:
                            ch.Draw('%s>>%s'%(entry.variable,th_name),weightdraw)  * ds.xsec * int_lumi / ds.totalEvents
                            #print '%s %s: %i' % (entry.name, ds.name, th.GetEntries())

                        if unscaled:
                            th.Write()
                        if not ds.isData:
                            if th:
                                th.Scale(group.kfactor * ds.xsec * int_lumi / ds.totalEvents)
                                if SM and ds.xpol:
                                    th.Scale( ds.xpol )
                                #print '%s : %.2f' % (ds.name, th.Integral(0,th.GetNbinsX()+1))
                            elif SFsys:
                                th_up.Scale(group.kfactor * ds.xsec * int_lumi / ds.totalEvents)
                                th_down.Scale(group.kfactor * ds.xsec * int_lumi / ds.totalEvents)
                                print '%s \t %s \t %s: %.2f' % (SFsys, 'UP', ds.name, th_up.Integral() )
                                print '%s \t %s \t %s: %.2f' % (SFsys, 'DOWN', ds.name, th_down.Integral() )                        

                                #print '%s UP: %.2f' % (ds.name, th_up.Integral(0,th_up.GetNbinsX()+1))
                                #print '%s DOWN: %.2f' % (ds.name, th_down.Integral(0,th_up.GetNbinsX()+1))
                        #print th.Integral(), th.GetEntries()
                        
                        if th:
                            gth.Add(th)
                        elif SFsys:
                            gth_up.Add(th_up)
                            gth_down.Add(th_down)
                        for rf in ds.root_files: rf.Close()
                        del ch
                    ## Progress Print 
                    if gth:
                        if not th2:
                            print '%s \t %s : %.2f' % (r.name, group.name,gth.Integral(0,gth.GetNbinsX()+1))
                            #print gth.GetBinContent(1)
                            if not th2 and loglady and i==0: logging.info('%s: %.2f' % (group.name,gth.Integral(0,gth.GetNbinsX()+1)))
                    elif SFsys:
                        print '%s \t %s \t %s: %.2f' % (SFsys, 'UP', group.name, gth_up.Integral() )
                        print '%s \t %s \t %s: %.2f' % (SFsys, 'DOWN', group.name, gth_down.Integral() )                        

                    ## Write output
                    if gth:
                        gth.Write()
                    elif SFsys:
                        gth_up.Write()
                        gth_down.Write()

        out.Write()
        logging.shutdown()


if __name__ == "__main__":

    if len(sys.argv) == 2:
        makeOut(getattr(run,sys.argv[1]))
    elif len(sys.argv) == 3 and sys.argv[2] == 'info':
        print '\n*** printing info ***'
        time.sleep(1)
        print_info(getattr(run,sys.argv[1]))
    else:
        print 'give me something to do'
        



