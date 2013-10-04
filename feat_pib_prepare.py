import os
import sys
from nipype.interfaces.base import CommandLine
from nipype.utils.filemanip import fname_presuffix
from nipype.interfaces import fsl
from glob import glob


def find_single_file(searchstring):
    """ glob for single file using searchstring
    if found returns full file path """
    file = glob(searchstring)
    if len(file) < 1:
        print '%s not found' % searchstring
        return None
    else:
        outfile = file[0]
        return outfile
        
        
def getsmooth(designfsf):
    with open(designfsf, 'r') as searchfile:
        for line in searchfile:
            try:
                if 'fmri(smooth)' in line:
                    smooth =  line.split()[-1]
                    return smooth
            except:
                print 'Could not find fwhm of smoothing applied to func data'
                print 'Make sure design.fsf is located in feat dir'


def get_resolution(img):
    cmd = ' '.join(['fslval', 
                    img, 
                    'pixdim1'])
    cout = CommandLine(cmd).run()
    if not cout.runtime.returncode == 0:
        print cmd
        print cout.runtime.stderr, cout.runtime.stdout
        return None
    else:
        res = cout.runtime.stdout
        return res


def est_smoothing(featdir, trgimage, stdspaceres=None):   
    print 'Estimating smoothness...' 
    designfsf = os.path.join(featdir, 'design.fsf')
    func_smoothing = getsmooth(designfsf) # Get fwhm smoothing of func data
    if not stdspaceres:
        subjstdimg = os.path.join(featdir,'reg', 
                                    'standard.nii.gz') # Standard space img
        stdspaceres = get_resolution(subjstdimg)
    example_func = os.path.join(featdir, 'example_func')
    cmd  = ' '.join(['match_smoothing', 
                    example_func, 
                    func_smoothing, 
                    trgimage, 
                    str(stdspaceres)])
    cout = CommandLine(cmd).run()
    if not cout.runtime.returncode == 0:
        print cmd
        print cout.runtime.stderr, cout.runtime.stdout
    else:
        smoothing = cout.runtime.stdout
        print 'Subject-space images will be smoothed by sigma=%smm'%(smoothing)
        print 'to match the standard-space functional data'
        return smoothing

        
def apply_smooth(img, smoothing):
    outfname = fname_presuffix(img, prefix=u's')
    cmd = ' '.join(['fslmaths', 
                    img, 
                    '-s %s'%(smoothing), 
                    outfname])
    cout = CommandLine(cmd).run()
    if not cout.runtime.returncode == 0:
        print cmd
        print cout.runtime.stderr, cout.runtime.stdout
    else:
        globstr = ''.join([outfname,'*'])
        outfile = find_single_file(globstr)
        return outfile
        
        
def flirt_coreg(infile, ref, outmat, outfile, cost='corratio'):     
    flt = fsl.FLIRT(bins=256, cost_func=cost)
    flt.inputs.in_file = infile
    flt.inputs.reference = ref
    flt.inputs.out_file = outfile
    flt.inputs.out_matrix_file = outmat
    flt.inputs.dof = 6
    cout = flt.run() 
    if not cout.runtime.returncode == 0:
        print invt.cmdline
        print cout.runtime.stderr, cout.runtime.stdout
        return None
    else:
        return cout.outputs

        
def applywarp(infile, ref, warp, premat, outfile):
    aw = fsl.ApplyWarp()
    aw.inputs.in_file = infile
    aw.inputs.ref_file = ref
    aw.inputs.field_file = warp
    aw.inputs.premat = premat
    aw.inputs.out_file = outfile
    cout = aw.run()   
    if not cout.runtime.returncode == 0:
        print invt.cmdline
        print cout.runtime.stderr, 
        return None
    else:
        return cout.outputs   
 
                
def regsummary(overlay, imglist):
    if type(imglist)==list:
        imglist = ' '.join(imglist)
    cmd = ' '.join(['slicesdir',
                    '-p %s'%(overlay),
                    imglist])
    cout = CommandLine(cmd).run()
    if not cout.runtime.returncode == 0:
        print cout.cmdline
        print cout.runtime.stderr, cout.runtime.stdout
        return None  
    else:
        print 'To view registration summary, point browser at:'
        print cout.runtime.stdout.split('\n')[-2:]
        

def concat4d(outfile, filelist):
    if type(filelist)==list:
        filelist = ' '.join(filelist)
    concatcmd = ' '.join(['fslmerge -t', 
                            outfile, 
                            filelist])
    concatout = CommandLine(concatcmd).run()
    if not concatout.runtime.returncode == 0:
        print concatcmd
        print concatout.runtime.stderr, concatout.runtime.stdout
        return None
    else:
        return outfile


def demean4d(infile, outfile):       
    demeancmd = ' '.join(['fslmaths', 
                            infile, 
                            '-Tmean -mul -1 -add %s'%(infile), 
                            outfile])
    demeanout = CommandLine(demeancmd).run()
    if not demeanout.runtime.returncode == 0:
        print demeancmd
        print demeanout.runtime.stderr, demeanout.runtime.stdout
        return None
    else:
        print '4D file saved to %s. Remember to verify subject order!'%(outfile)
        print 'You may want to add additional smoothing to 4D file in order'
        print 'to ameliorate possible effects of mis-registrations between' 
        print 'functional and structural data, and to lessen the effect of'
        print 'the additional confound regressors'
        return outfile
    



        

    
