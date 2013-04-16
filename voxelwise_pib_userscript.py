import sys, os
sys.path.insert(0, '/home/jagust/jelman/CODE/pib_to_fsl')
import coreg_pib2fsl as coregpib
import feat_pib_prepare as pibprep
from nipype.utils.filemanip import split_filename
from nipype.utils.filemanip import save_json



        
if __name__ == '__main__':

    # Get list of subjects from commandline call
    if len(sys.argv) > 1:   #If specified, load file as list of subject id's
        sublist = sys.argv[1:]
    else:    #If no file name given, search for data
        print 'No subjects specified'
        print 'Usage voxelwise_pib_userscript.py <list of subject codes>'
        
    ######################################################################
    ################# Specify path and file names ########################
    ### Specify path names
    basedir = '/home/jagust/rsfmri_ica/data' # main project directory
    pibdir = 'pib' # name of directory containing pib data
    anatdir = 'anat' # name of directory containing FSL structural data
    funcdir = 'func' # name of directory containing FSL functional data
    featdir = '_4d_OldICA_IC0_ecat_2mm_6fwhm_125.ica'
    
    ### Inputs ###
    # Specify filenames and filename patterns from pib processing
    matbrain_fname = 'brainmask.nii*' # name of mri file used in pib processing
    pib_fname_pattern = 'DVR*.nii*' # filename pattern of pib image
    mat_coregname = 'mri_to_pet.mat' # mat file used to coreg mri->pet
    # Specify filenames of files used in FSL processing stream
    flirtbrain_fname = 'anat_brain.nii.gz' # structural image from FSL processing
    flirtwarp_fname = 'highres2standard_warp.nii.gz'
    stdbrain = os.path.join(os.getenv('FSLDIR'),
                            'data/standard',
                            'MNI152_T1_2mm_brain.nii.gz') # Standard brain used
                            
    ### Outputs ###                            
    # Specify names of FSL readable matfiles
    flirt_coregname = 'mri_to_pet_flirt.mat' # name for FSL readable copy of mri->pet coreg
    invflirt_coregname = 'pet_to_mri_flirt.mat' # name for pet-> mat coreg
    struc_coregname = 'brainmask_to_anat_flirt.mat' # name for matbrain->flirtbrain coreg
    concat_coregname = 'pet_to_anat_flirt.mat' # name for pet->flirtbrain coreg
    # Specify name of pib file registered to standard space
    stdpib_fname = 'DVR2std.nii.gz' 
    pib4d_fname = os.path.join(basedir,
                        'OldICA_IC0_ecat_2mm_6fwhm_125.gica',
                        '4D_pib_dvr.nii.gz')
    ###########################################################################
    ###########################################################################
    
    startdir = os.getcwd()
    os.chdir(basedir)

    warpedpib = {} # create dict to save list of subjects and registered pib files
    
    ### Begin looping over all subject
    ##################################
    for subj in sublist:
        print 'Beginning subject %s'%(subj)
        
        # Sets subject specific subdirectories
        subjpibdir = os.path.join(basedir, subj, pibdir)
        subjanatdir = os.path.join(basedir, subj, anatdir)
        subjfuncdir = os.path.join(basedir, subj, funcdir)
        subjfeatdir = os.path.join(basedir, subj, funcdir, 
                                ''.join([subj,featdir]))
        
        # Estimate and apply smoothing to PIB images
        globstr = os.path.join(subjpibdir, pib_fname_pattern) 
        pibfile = coregpib.find_single_file(globstr) #raw pib image in subject space
        smoothing = pibprep.est_smoothing(subjfeatdir, pibfile)
        smoothedpib = pibprep.apply_smooth(pibfile, smoothing)
        
        # Sets inputs to worldmat2flirt.m
        globstr = os.path.join(subjpibdir, matbrain_fname)
        matbrain =  coregpib.find_single_file(globstr) #mri used in pib processing
        globstr = os.path.join(subjpibdir, 
                            mat_coregname)
        mat_coreg = coregpib.find_single_file(globstr) #matfile from mri->pet coreg 
        flirt_coreg = os.path.join(subjpibdir, flirt_coregname)
        
        # Convert mri->pet matfile from matlab to flirt format
        flirt_coreg_out = coregpib.mat2flirt(mat_coreg,matbrain,pibfile,flirt_coreg)
        
        # Get inverted transform using convert_xfm to map pet-> mri
        invflirt_coreg = os.path.join(subjpibdir, invflirt_coregname)
        invflirt_out = coregpib.invert_coreg(flirt_coreg, invflirt_coreg)

        # Register brainmask to structural used in FSL processing stream
        sub_flirtbrain_fname = '_'.join([subj, flirtbrain_fname]) # set subj specific filename
        flirtbrain = os.path.join(subjanatdir, sub_flirtbrain_fname)
        mattoflirt_coreg = os.path.join(subjpibdir, struc_coregname)
        mattoflirt_coreg_out = coregpib.flirt_coreg(matbrain, flirtbrain, mattoflirt_coreg)

        # Concatenate pib->brainmask with brainmask->anat
        concat_coreg = os.path.join(subjpibdir, concat_coregname)
        concat_coreg_out = coregpib.concat_xfm(invflirt_coreg, mattoflirt_coreg, concat_coreg)

        # Warp pib->std. Uses existing structural->std warp and pib->structural premat
        flirtwarp = os.path.join(subjfeatdir,
                                'reg',
                                flirtwarp_fname)
        stdpib = os.path.join(subjpibdir, stdpib_fname)
        warp_pib2std = coregpib.applywarp(pibfile, stdbrain, flirtwarp, concat_coreg, stdpib)
        
        # update dict with subj and file
        warpedpib.update({subj:warp_pib2std.out_file}) 
        
        
    # Create summary of pet->std registration for review
    pibfilelist = []
    for key in sorted(warpedpib):
        pibfilelist.append(warpedpib[key])
    coregpib.regsummary(stdbrain, pibfilelist)
    
    # Concatenate into 4D file and demean for use as covariate
    pib4dfile = pibprep.concat_demean(pib4d_fname, pibfilelist)
    
    # Save out dict listing subject order and files in 4D file
    pth, fname_base, ext = split_filename(pib4d_fname)
    fname = ''.join([os.path.join(pth,fname_base), '_subjorder.txt'])
    with open(fname, 'w+') as f:
        f.write('\n'.join(pibfilelist))
        
    os.chdir(startdir)
