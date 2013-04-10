import os
import sys
from nipype.interfaces.base import CommandLine

def getsmooth(designfsf):
    with open(designfsf, 'r') as searchfile:
        for line in searchfile:
            if 'fmri(smooth)' in line:
                smooth =  line.split()[-1]
                return smooth

if __name__ == '__main__':
    
    if len(sys.argv) <= 2:
        print 'Usage: python feat_pib_prepare.py <4D-PIB-output> <list of subject codes>'
        print 'Note: these all have to have had registration completed on them'
        sys.exit()
    else:
        args = sys.argv[1:]
        
    pib4d = args[0]
    sublist = args[1:]
    
    ######################################################################
    ################# Specify path and file names ########################
    ### Specify path names
    basedir = '/home/jagust/rsfmri_ica/data' # main project directory
    pibdir = 'pib' # name of directory containing pib data
    featdir = '_4d_OldICA_IC0_ecat_2mm_6fwhm_125.ica' # name of feat directory
    
    for subj in sublist:
        print 'Starting on subject %s'%(subj)
        
        #Subject specific paths
        subjpibdir = os.path.join(basedir, subj, pibdir)
        subjfuncdir = os.path.join(basedir, funcdir)
        subjfeatdir = os.path.join(basedir, subj, 
                                    'func', 
                                    ''.join([subj, featdir]))
                                    
        # estimate how much we will need to smooth the pib data by, in the end
        print 'Estimating smoothness...'
        designfsf = os.path.join(subjfeatdir, 'design.fsf')
        func_smoothing = getsmooth(designfsf)
        
        
    
