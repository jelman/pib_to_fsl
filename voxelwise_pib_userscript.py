import sys, os
sys.path.insert(0, '/home/jagust/jelman/CODE/voxelwise_pib')
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
    basedir = '/home/jagust/DST/FSL' # main project directory
    pibdir = os.path.join(basedir,'pib','%(subj)s') # name of directory containing pib data
    anatdir = os.path.join(basedir,'structural','%(subj)s') # name of directory containing FSL structural data
    featdir = os.path.join(basedir,'functional/%(subj)s/run02/Detail.feat')
    
    ### Inputs ###
    ## Specify filenames and filename patterns from pib processing
    dvr_pattern = os.path.join(pibdir,'rawdvr.nii') # filename pattern of pib image
    meandvr_pattern = os.path.join(pibdir,'mean20min.nii.gz') # filename pattern of mean 20 min pib image
    ## Specify filenames of files used in FSL processing stream
    highres_pattern = os.path.join(anatdir,'T1_brain.nii.gz') # structural image from FSL processing
    refbrain = os.path.join('/home/jagust/jelman/templates/MNI/data/standard',
                        'MNI152_T1_2mm_brain.nii.gz') # Standard brain used
    #Specify standard space smoothing resolution. Otherwise, set to None in order to
    #calculate from image.
    stdspaceres = None 
    # Determinr whether to use linear or non-linear registration
    # If linear, set highres2std_pattern as mat file. If non-linear, use warp image
    nonlinear_reg = False
    highres2std_pattern = os.path.join(featdir,
                                        'reg', 'highres2standard.mat')

                    
    ### Outputs ###                            
    dvr2highres_mat = 'dvr2highres.mat' # name for pet-> highres coreg mat
    dvr2highres_fname = 'mean20min2highres.nii.gz'
    dvr2std_mat = 'dvr2std.mat' # name for pet > std coreg mat (not used for non-linear coreg)
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
        subjpibdir = pibdir%{"subj":subj}
        subjfeatdir = featdir%{"subj":subj}
        subjhighres = highres_pattern%{"subj":subj}
        subjdvr = dvr_pattern%{"subj":subj}
        subjmeandvr = meandvr_pattern%{"subj":subj}
        subjhighres = highres_pattern%{"subj":subj}
        subjhighres2std_warp = highres2std_pattern%{"subj":subj}
        subjdvr2std_mat = os.path.join(subjpibdir, dvr2std_mat)
        stdpib = os.path.join(subjpibdir, dvr2std_fname)
        
        # Estimate and apply smoothing to PIB images
        smoothing = pibprep.est_smoothing(subjfeatdir, subjdvr, stdspaceres)
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


        if nonlinear_reg:
            warp_pib2std = pibprep.applywarp(smoothedpib, 
                                                refbrain, 
                                                highres2std_warp, 
                                                dvr2highres_outmat, 
                                                stdpib)
        else:
            dvr2std_outmat = pibprep.concat_xfm(dvr2highres_outmat,
                                        subjhighres2std_warp,
                                        subjdvr2std_mat)
            warp_pib2std = 
        
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
