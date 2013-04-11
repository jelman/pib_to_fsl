import os, sys, re
from glob import glob
import nipype.interfaces.fsl as fsl
from nipype.interfaces.base import CommandLine
import worldmat2flirtmap_pywrapper as w2f
from nipype.utils.filemanip import split_filename

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

def unzip_file(infile):
    """ looks for gz  at end of file,
    unzips and returns unzipped filename"""
    base, ext = os.path.splitext(infile)
    if not ext == '.gz':
        return infile
    else:
        cmd = CommandLine('gunzip %s' % infile)
        cout = cmd.run()
        if not cout.runtime.returncode == 0:
            print 'Failed to unzip %s'%(infile)
            return None
        else:
            return base
       
def get_subid(instr):
    """ pulls the subid Bxx-xxx out of a string
    returns text matching pattern"""
    m = re.search('B[0-9]{2}-[0-9]{3}', instr)
    try:
        return m.group()
    except:
        return None

def mat2flirt(worldmat,src,trg,flirtmat):
    src = unzip_file(src)
    trg = unzip_file(trg)
    mlabcmd = w2f.worldmat2flirtmap()
    mlabcmd.inputs.src = src
    mlabcmd.inputs.trg = trg
    mlabcmd.inputs.worldmat = worldmat
    mlabcmd.inputs.output_file = flirtmat
    mlabout = mlabcmd.run()  
    if not mlabout.runtime.returncode == 0:
        print mlabout.runtime.stderr, mlabout.runtime.stdout
        return None
    else:
        return mlabout.outputs.output_file
        
    
def invert_coreg(mat, invmat):
    invt = fsl.ConvertXFM()
    invt.inputs.in_file = mat
    invt.inputs.invert_xfm = True
    invt.inputs.out_file = invmat
    cout = invt.run()
    if not cout.runtime.returncode == 0:
        print invt.cmdline
        print cout.runtime.stderr, cout.runtime.stdout
        return None
    else:
        return cout.outputs

def flirt_coreg(infile, ref, outmat):     
    flt = fsl.FLIRT(bins=256, cost_func='normcorr')
    flt.inputs.in_file = infile
    flt.inputs.reference = ref
    flt.inputs.out_matrix_file = outmat
    flt.inputs.dof = 6
    cout = flt.run() 
    if not cout.runtime.returncode == 0:
        print invt.cmdline
        print cout.runtime.stderr, cout.runtime.stdout
        return None
    else:
        return cout.outputs

def concat_xfm(infile1, infile2, outfile):
    concat = fsl.ConvertXFM()
    concat.inputs.in_file = infile1
    concat.inputs.in_file2 = infile2
    concat.inputs.concat_xfm = True
    concat.inputs.out_file = outfile
    cout = concat.run()
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

    
if __name__ == '__main__':

    # Get list of subjects from commandline call
    if len(sys.argv) > 1:   #If specified, load file as list of subject id's
        sublist = sys.argv[1:]
    else:    #If no file name given, search for data
        print 'No subjects specified'
        print 'coreg_pib2fsl.py <list of subject codes>'
        
    ######################################################################
    ################# Specify path and file names below########################
    
    ### Path names ###
    basedir = '/home/jagust/rsfmri_ica/data' # main project directory
    pibdir = 'pib' # name of directory containing pib data
    anatdir = 'anat' # name of directory containing FSL structural data
    funcdir = 'func' # name of directory containing FSL functional data
    
    ### Inputs ###
    # Filenames and filename patterns from pib processing
    matbrain_fname = 'brainmask.nii*' # name of mri file used in pib processing
    pib_fname_pattern = 'DVR*.nii*' # filename pattern of pib image
    mat_coregname = 'mri_to_pet.mat' # mat file used to coreg mri->pet
    # Filenames of files used in FSL processing stream
    flirtbrain_fname = 'anat_brain.nii.gz' # structural image from FSL processing
    flirtwarp_fname = os.path.join('*_4d_OldICA_IC0_ecat_2mm_6fwhm_125.ica',
                                'reg',
                                'highres2standard_warp.nii.gz') # filename pattern of
                                                                # flirtbrain -> std warp
    stdbrain = os.path.join(os.getenv('FSLDIR'),
                            'data/standard',
                            'MNI152_T1_2mm_brain.nii.gz') # Standard brain used
                            
    ### Outputs ###                            
    # Specify names of FSL readable matfiles
    flirt_coregname = 'mri_to_pet_flirt.mat' # name for FSL readable copy of mri->pet coreg
    invflirt_coregname = 'pet_to_mri_flirt.mat' # name for pet-> mat coreg
    struc_coregname = 'brainmask_to_anat.mat' # name for matbrain->flirtbrain coreg
    concat_coregname = 'pet_to_anat_flirt.mat' # name for pet->flirtbrain coreg
    # Specify name for pib file after registration to standard space
    stdpib_fname = 'DVR2std.nii.gz' 
    ###########################################################################
    ###########################################################################
    
    startdir = os.getcwd()
    os.chdir(basedir)

    ### Begin looping over all subjects
    ###################################
    for subj in sublist:
        
        # Sets subject specific subdirectories
        subjpibdir = os.path.join(basedir, subj, pibdir)
        subjanatdir = os.path.join(basedir, subj, anatdir)
        subjfuncdir = os.path.join(basedir, subj, funcdir)
        
                
        ### Begin processing
        ######################
       
        # Sets inputs to worldmat2flirt.m
        globstr = os.path.join(subjpibdir, matbrain_fname)
        matbrain =  find_single_file(globstr) 
        matbrain = unzip_file(matbrain) #mri used in pib processing
        globstr = os.path.join(subjpibdir, pib_fname_pattern) 
        pibfile = find_single_file(globstr) 
        pibfile = unzip_file(pibfile) #raw pib image in subject space
        globstr = os.path.join(subjpibdir, 
                            mat_coregname)
        mat_coreg = find_single_file(globstr) #matfile from mri->pet coreg 
        flirt_coreg_out = os.path.join(subjpibdir, flirt_coregname)
        # Convert mri->pet matfile from matlab to flirt format
        flirt_coreg = mat2flirt(mat_coreg,matbrain,pibfile,flirt_coreg_out)

        # Get inverted transform using convert_xfm to map pet-> mri
        invflirt_out_name = os.path.join(subjpibdir, invflirt_coregname)
        invflirt_coreg = invert_coreg(flirt_coreg, invflirt_out)

        # Register brainmask to structural used in FSL processing stream
        sub_flirtbrain_fname = '_'.join([subj, flirtbrain_fname]) # set subj specific filename
        flirtbrain = os.path.join(subjanatdir, sub_flirtbrain_fname)
        mattoflirt_name = os.path.join(subjpibdir, struc_coregname)
        coreg_mat2flirt = flirt_coreg(matbrain, flirtbrain, mattoflirt_name)

        
        # Concatenate pib->brainmask with brainmask->anat
        concat_coreg = os.path.join(subjpibdir, concat_coregname)
        pibtoflirt_mat = concat_xfm(invflirt_coreg, matbrain_to_flirtbrain, concat_coreg)

        # Apply existing 
        flirtwarp = find_single_file(os.path.join(subjfuncdir,flirtwarp_fname))
        stdpib = os.path.join(subjpibdir, stdpib_fname)
        warp_pib2std = applywarp(pibfile, stdbrain, flirtwarp, concat_coreg, stdpib)
        
        # Create summary of pet->std registration for review
        warpedpib.update({subj:warp_pib2stdout_file})
        regsummary(stdbrain, warpedpib.values())

    os.chdir(startdir)
