import os, sys, re
from glob import glob


def unzip(infile):
    gunzipfile, gz = os.path.splitext(infile)
    if not 'gz' in gz:
        #when running gunzip on file when
        return infile
    else:
       c3 = CommandLine('gunzip %s'%(infile))
       c3.run()
       return gunzipfile


def get_subid(instr):
    """ pulls the subid Bxx-xxx out of a string
    returns text matching pattern"""
    m = re.search('B[0-9]{2}-[0-9]{3}', instr)
    try:
        return m.group()
    except:
        return None


if __name__ == '__main__':

    ###Specify project directory containing subject dirs.
    basedir = '/home/jagust/rsfmri_ica/data'

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
        
        ## Specify subject subdirectories
        pibdir = os.path.join(basedir, subj, 'pib')
        anatdir = os.path.join(basedir, subj, 'anat')

