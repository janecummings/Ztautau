import os,time,logging
from spec import ANALYSIS_DIR


def atime():
    return '%s' % time.strftime("%m_%d_%H_%M",time.localtime())


def makeLog(OUTPUTDIR,name='',comment=''):
    logger = logging.getLogger(OUTPUTDIR)
    logname = name or OUTPUTDIR.split('/')[-1]
    if not os.path.exists('%s/run/log/%s' % (ANALYSIS_DIR,logname)): os.mkdir('%s/run/log/%s' % (ANALYSIS_DIR,logname))
    logname = '%s/run/log/%s/%s.log' % (ANALYSIS_DIR,logname,atime())
    logging.basicConfig(filename = logname,
                        format = ' %(message)s',
                        level = logging.DEBUG)
    logging.info('%s: %s' %(atime(),comment))
    logging.info('OUTPUTDIR = %s' % OUTPUTDIR)
    
