import os, sys, re
from glob import glob
import scipy.io as sio
import numpy as np
import nibabel as nib
import nipype.interfaces.fsl as fsl
from nipype.interfaces.base import CommandLine
import worldmat2flirtmap_pywrapper as w2f


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


    
if __name__ == '__main__':

    ######################################################################
    ################# Specify path and file names ########################
    ### Specify path names
    basedir = '/home/jagust/rsfmri_ica/data' # main project directory
    pibdir = 'pib' # name of directory containing pib data
    anatdir = 'anat' # name of directory containing FSL structural data
    funcdir = 'func' # name of directory containing FSL functional data
    
    ### Inputs ###
    # Specify filenames and filename patterns from pib processing
    matbrain_fname = 'brainmask.nii*' # name of mri file used in pib processing
    pib_fname_pattern = 'DVR*.nii*' # filename pattern of pib image
    mat_coregname = 'mri_to_pet.mat' # mat file used to coreg mri->pet
    # Specify filenames of files used in FSL processing stream
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
    concat_coregname = 'pet_to_anat_flirt.mat' # name for pet->flirtbrain coreg
    # Specify name of pib file registered to standard space
    stdpib_fname = 'DVR2std.nii.gz' 
    ###########################################################################
    ###########################################################################
    
    # Specify list of subject codes on command line, otherwise searches basedir for data
    if len(sys.argv) > 1:   #If specified, load file as list of subject id's
        sublist = sys.argv[1:]
    else:    #If no file name given, search for data
        print 'No subjects specified, searching for data in ', basedir
        # look for subjects already registered filtered func ica data
        globstr = os.path.join(basedir,'B*')   #subject id 
        # List of paths to subjects dirs
        allsubdirs = glob(globstr)
        # Create list of subject id's
        sublist = []
        for subdir in allsubdirs:
            subid = get_subid(subdir)
            sublist.append(subid)
            if subid is None:
                print subdir, 'does not have a valid subid'
                continue
    
    ### Begin looping over all subject
    ##################################
    for subj in sublist:
        
        # Sets subject specific subdirectories
        subjpibdir = os.path.join(basedir, subj, pibdir)
        subjanatdir = os.path.join(basedir, subj, anatdir)
        subjfuncdir = os.path.join(basedir, subj, funcdir)
        
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
        flirt_coreg = os.path.join(subjpibdir, flirt_coregname)
        
        ###Convert mri->pet matfile from matlab to flirt format
        #Outputs firtmat specified above
        mattoflirt = w2f.worldmat2flirtmap()
        mattoflirt.inputs.src = matbrain
        mattoflirt.inputs.trg = pibfile
        mattoflirt.inputs.worldmat = mat_coreg
        mattoflirt.inputs.output_file = flirt_coreg
        mlabrun = mattoflirt.run()
        
        # Get inverted transform using convert_xfm to map pet-> mri
        invflirt_coreg = os.path.join(subjpibdir, invflirt_coregname)
        
        invt = fsl.ConvertXFM()
        invt.inputs.in_file = flirt_coreg
        invt.inputs.invert_xfm = True
        invt.inputs.out_file = invflirt_coreg
        invt.run()
        
        # Register brainmask to structural used in FSL processing stream
        sub_flirtbrain_fname = '_'.join([subj, flirtbrain_fname]) # set subj specific filename
        flirtbrain = os.path.join(subjanatdir, sub_flirtbrain_fname)
        matbrain_to_flirtbrain = os.path.join(subjpibdir, 'brainmask_to_anat.mat')

        flt = fsl.FLIRT(bins=256, cost_func='normcorr')
        flt.inputs.in_file = matbrain
        flt.inputs.reference = flirtbrain
        flt.inputs.out_matrix_file = matbrain_to_flirtbrain
        flt.inputs.dof = 6
        flt.run() 
        
        # Concatenate pib->brainmask with brainmask->anat
        concat_coreg = os.path.join(subjpibdir, concat_coregname)
        
        concat = fsl.ConvertXFM()
        concat.inputs.in_file = invflirt_coreg
        concat.inputs.in_file2 = matbrain_to_flirtbrain
        concat.inputs.concat_xfm = True
        concat.inputs.out_file = concat_coreg
        concat.run()
        
        # Apply existing 
        flirtwarp = find_single_file(os.path.join(subjfuncdir,flirtwarp_fname))
        stdpib = os.path.join(subjpibdir, stdpib_fname)
        
        aw = fsl.ApplyWarp()
        aw.inputs.in_file = pibfile
        aw.inputs.ref_file = stdbrain
        aw.inputs.field_file = flirtwarp
        aw.inputs.premat = concat_coreg 
        aw.inputs.out_file = stdpib
        aw.run()

        
