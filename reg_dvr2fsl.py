import os, sys, re
from glob import glob
import scipy.io as sio
import numpy as np
import nibabel as nib
import nipype

def get_subid(instr):
    """ pulls the subid Bxx-xxx out of a string
    returns text matching pattern"""
    m = re.search('B[0-9]{2}-[0-9]{3}', instr)
    try:
        return m.group()
    except:
        return None

def nifti2scl(img):
    imgmat = img.get_affine() # get image affine
    # not sure if this is always correct with rotations in mat, but seems okay! - Ged
    imgmat_sumsq = np.sum(imgmat[0:3,0:3]**2, axis=0)
    imgmat_rtsumsq = np.sqrt(imgmat_sumsq)
    scl = np.diag(np.append(imgmat_rtsumsq,1))
    if np.linalg.det(imgmat) > 0:
        # neurological, x-axis is flipped, such that [3 2 1 0] and [0 1 2 3]
        # have the same *scaled* coordinates - Ged
        xflip = np.diag((-1, 1, 1, 1))
        xflip[0,3] = img.shape[0] - 1 # reflect about centre
        scl = np.dot(scl, xflip)
    return scl
    
def worldmat2flirtmat(worldmat, src, trg):
    srcimg = nib.load(src)
    trgimg = nib.load(trg)
    worldmat_contents = sio.loadmat(worldmat)  
    # get affine of files
    worldmat_M = worldmat_contents['M'] 
    srcmat = srcimg.get_affine()
    trgmat = trgimg.get_affine()
    
    # get voxel matrices, voxels are one-based in Matlab and zero-based in FSL
    spmvoxmat = np.dot(np.linalg.inv(srcmat), np.dot(worldmat_M,trgmat))
    spmvoxmat = np.around(spmvoxmat, decimals=4)
    addone = np.identity(4)
    addone[:,3] = 1
    fslvoxmat = np.dot(np.linalg.inv(addone), np.dot(spmvoxmat,addone))
    fslvoxmat = np.around(fslvoxmat,decimals=4)

    # get flirt scaled elements
    trgscl = nifti2scl(trgimg)
    trgscl = np.around(trgscl, decimals=4)
    srcscl = nifti2scl(srcimg)
    srcscl = np.around(srcscl, decimals=4)
    inv_trgscl = np.linalg.inv(trgscl)
    flirtmat = np.linalg.inv(np.dot(srcscl, np.dot(fslvoxmat, inv_trgscl)))
    flirtmat = np.around(flirtmat, decimals=4)
    
    return flirtmat, spmvoxmat, fslvoxmat
    
    
    
if __name__ == '__main__':

    ### Specify path names
    basedir = '/home/jagust/rsfmri_ica/data' # main project directory
    pibdir = 'pib' # name of directory containing pib data
    anatdir = 'anat' # name of directory containing FSL structural data
    
    ### Specify filenames used in pib processing
    mri_pibproc = 'brainmask.nii' # name of mri file used in pib processing
    pib_fname_pattern = 'mean20min*.nii*' # filename pattern of pib image
    coreg_mat = 'mri_to_pet.mat' # mat file used to coreg mri->pet
    anat_fname = 'anat_brain.nii.gz' # file name of structural image
                                     #   used in FSL processing
    
    ### Specify filenames of flirt-style mat files
    flirtmatout = 'mri_to_pet_flirt.mat' # FSL readable copy of coreg_mat
    invflirtmatout = 'pet_to_mri_flirt.mat' # FSL readable copy of inverted coreg_mat

    ### Specify list of subject codes on command line, otherwise searches basedir for data
    if len(sys.argv) ==2:   #If specified, load file as list of subject id's
        args = sys.argv[1]
        print 'Using subjects specified in ', args
        with open(args, 'r') as f:
            sublist = f.read().splitlines()
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
                
    ### Loop over all subjects.
    for subj in sublist:
        
        ## Specify subject specific subdirectories
        subjpibdir = os.path.join(basedir, subj, pibdir)
        subjanatdir = os.path.join(basedir, subj, anatdir)
        
        ## Specify inputs to worldmat2flirt.m
        brainmask = os.path.join(subjpibdir, 
                            mri_pibproc) #mri used in pib processing
        pibglob = os.path.join(subjpibdir, 
                            pib_fname_pattern) # filename pattern of pib
        pibfile = glob(pibglob)[0] #raw pib image in subject space
        matfile = os.path.join(subjpibdir, 
                            coreg_mat) #matfile from mri->pet coreg 
        ###Convert mri->pet matfile from matlab to flirt format
        flirtmat, spmvoxmat, fslvoxmat = worldmat2flirtmat(matfile, brainmask, pibfile)
        invflirtmat = np.linalg.inv(flirtmat) #get inverted transform
        # Save out flirtmat and invflirtmat
        np.savetxt(flirtmatout, flirtmat,fmt='%10.4f')
        np.savetxt(invflirtmatout, invflirtmat,fmt='%10.4f')
        
        # Register brainmask to structural used in FSL processing stream
        anat_brain = os.path.join(subjanatdir, anat_fname)
        
        
        """
        TO DO:
        X save out flirtmat (mri -> pet)
        X invert flirtmat to get pet -> mri)
        - test to verify it maps between mri and pet
        
        - flirt brainmask->flirt_brain
        - concat pet->mri &  brainmask->anat_brain
        - apply anat_brain2std warp to pib image,
        use concat mat as initial linear transform
        - adapt feat_gm_prepare to create 4D voxelwise pib covariate
        
