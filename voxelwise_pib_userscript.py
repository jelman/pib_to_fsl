import sys, os
sys.path.insert(0, '/home/jagust/jelman/CODE/pib_to_fsl')
import feat_pib_prepare as pibprep
from nipype.utils.filemanip import split_filename
from nipype.utils.filemanip import save_json



        
if __name__ == '__main__':

    # Get list of subjects from commandline call
    if len(sys.argv) > 1:   #If specified, load file as list of subject id's
        sublist = sys.argv[1:]
    else:    
        print 'No subjects specified'
        print 'Usage voxelwise_pib_userscript.py <list of subject codes>'
        
    ######################################################################
    ################# Specify path and file names ########################
    ### Specify path names
    basedir = '/home/jagust/rsfmri_ica/data' # main project directory
    pibdir = 'pib' # name of directory containing pib data
    anatdir = 'anat' # name of directory containing FSL structural data
    funcdir = 'func' # name of directory containing FSL functional data
    featdir = '_4d_OldYoungICA.ica'
    
    ### Inputs ###
    # Specify filenames and filename patterns from pib processing
    dvr_suffix = '_dvr.nii' # filename pattern of pib dvr image
    meandvr = '_mean20min.nii' #filename pattern of pib mean20min image
    # Specify filenames of files used in FSL processing stream
    highres_suffix = '_anat_brain.nii.gz' # structural image from FSL processing
    highres2std_fname = 'highres2standard_warp.nii.gz'
    refbrain = os.path.join('/home/jagust/jelman/templates/MNI/data/standard',
                            'MNI152_T1_3mm_brain.nii.gz') # Standard brain used
                            
    ### Outputs ###                            
    dvr2highres_mat = 'dvr2highres.mat' # name for pet-> highres coreg mat
    dvr2highres_fname = 'dvr2highres.nii.gz'
    dvr2std_fname = 'dvr2std.nii.gz'    #dvr image in standard space
    pib4d_fname = os.path.join(
                        basedir,        #name of 4d voxelwise covariate file
                        '4D_pib_dvr.nii.gz')
    pib4d_demean_fname = os.path.join(
                        basedir,        #name of demeaned 4d voxelwise covariate file
                        '4D_pib_dvr_demeaned.nii.gz')
    ###########################################################################
    ###########################################################################
    
    startdir = os.getcwd()
    os.chdir(basedir)

    warpedpib = {} # create dict to save list of subjects and registered pib files
    
    ### Begin looping over all subject
    ##################################
    for subj in sublist:
        print 'Beginning subject %s'%(subj)
        
        # Sets subject specific subdirectories and files
        subjpibdir = os.path.join(basedir, subj, pibdir)
        subjanatdir = os.path.join(basedir, subj, anatdir)
        subjfuncdir = os.path.join(basedir, subj, funcdir)
        subjfeatdir = os.path.join(basedir, subj, funcdir, 
                                ''.join([subj,featdir]))
        subjdvr = os.path.join(subjpibdir, 
                                ''.join([subj, dvr_suffix]))
        subjmeandvr = os.path.join(subjpibdir, 
                                ''.join([subj, meandvr_suffix]))
        subjhighres = os.path.join(subjanatdir, 
                                    ''.join([subj, highres_suffix]))
        
        # Estimate and apply smoothing to PIB images
        smoothing = pibprep.est_smoothing(subjfeatdir, subjdvr, 3)
        smoothedpib = pibprep.apply_smooth(subjdvr, smoothing)
        
        # Coregister dvr->highres and save out matfile
        print 'Running coregistrations...'
        subjoutmat = os.path.join(subjpibdir, dvr2highres_mat)
        subjoutfile = os.path.join(subjpibdir, dvr2highres_fname)
        dvr2highres_coreg = pibprep.flirt_coreg(subjmeandvr, 
                                        subjhighres, 
                                        subjoutmat,
                                        subjoutfile)
        dvr2highres_outmat = dvr2highres_coreg.out_matrix_file
        # Warp pib->std. Uses existing structural->std warp and pib->structural premat
        print 'Warping image to standard space...'
        highres2std_warp = os.path.join(subjfeatdir,
                                'reg',
                                highres2std_fname)
        stdpib = os.path.join(subjpibdir, dvr2std_fname)
        warp_pib2std = pibprep.applywarp(subjdvr, 
                                            refbrain, 
                                            highres2std_warp, 
                                            dvr2highres_outmat, 
                                            stdpib)
        
        # update dict with subj and file
        warpedpib.update({subj:warp_pib2std.out_file}) 
        
        
    # Create summary of pet->std registration for review
    dvrlist = []
    for key in sorted(warpedpib):
        dvrlist.append(warpedpib[key])
    pibprep.regsummary(refbrain, dvrlist)
    
    # Concatenate into 4D file and demean for use as covariate
    pib4dfile = pibprep.concat4d(pib4d_fname, dvrlist)
    pib4dfile_demeaned = pibprep.demean4d(pib4dfile, pib4d_demeaned_fname)
    # Save out dict listing subject order and files in 4D file
    pth, fname_base, ext = split_filename(pib4d_fname)
    fname = ''.join([os.path.join(pth,fname_base), '_subjorder.txt'])
    with open(fname, 'w+') as f:
        f.write('\n'.join(dvrlist))
        
    os.chdir(startdir)
