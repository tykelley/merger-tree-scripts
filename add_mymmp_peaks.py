#!/usr/bin/pytthon

import sys

def make_pbar(maxi,label=None):
    """
    Returns a progressbar instance, so I don't have to keep copying and pasting that code to make it
    """

    try:
        from progressbar import ProgressBar,Bar,Percentage,ETA
    except ImportError:
        return None


    if label is not None:
        widgets = ['{0} ('.format(label),Percentage(),') ',Bar(left='[',right=']',marker='-'),' ',ETA()]
    else:
        widgets = ['(',Percentage(),') ',Bar(left='[',right=']',marker='-'),' ',ETA()]
    return ProgressBar(widgets=widgets,maxval=maxi)





usage = "usage:  {0} <input HDF5 file>".format(sys.argv[0].split('/')[-1])

if len(sys.argv) < 2:
    print usage
    sys.exit(1337)
    
    
import h5py
from numpy import max,empty
#from mytools import make_pbar

fname = sys.argv[1]

print "\nOpening "+fname

f = h5py.File(fname)

mt = f['RockstarMergerTrees']

snaps = [k for k in f.keys() if k.startswith('Snapshot')]
snaps.sort()
snap = f[snaps[-1]]

cat = snap['HaloCatalog_RockstarMergerTree']

print "Adding 'Vpeak_TotalMass', 'Mpeak_TotalMass', and 'apeak_TotalMass' to {0}".format(cat.name)

tnums = cat['TreeNumber'][:]

Mpeakar = empty(len(tnums))
Mpeakar[:] = -1

Vpeakar = empty(len(tnums))
Vpeakar[:] = -1

apeakar = empty(len(tnums))
apeakar[:] = -1

pbar = make_pbar(tnums.shape[0],'Looping over all trees')
if pbar is not None:
    pbar.start()

for ii in range(tnums.shape[0]):
    tree = mt['Tree_{0}'.format(tnums[ii])]
    
    scale = tree['scale'][:]
    mymmp = tree['TotalMass_mmp'][:]
    if 'orig_mvir' in tree.keys():
        mvir = tree['orig_mvir'][:]
    else:
        print "No 'orig_mvir' dataset found; using 'mvir'" 
        mvir = tree['Mvir'][:]
    
    vmax = tree['vmax'][:]
    
    mainbranch = mymmp == 1
    
    mb_mvir = mvir[mainbranch]
    mb_vmax = vmax[mainbranch]
    mb_scale = scale[mainbranch]
    
    if len(mb_mvir) == 0:
        Mpeakar[ii] = -2
        Vpeakar[ii] = -2
        apeakar[ii] = -2
        if pbar is not None:
            pbar.update(ii)
        continue
    
    mpeak = max(mb_mvir)
    Mpeakar[ii] = mpeak
    
    peak = mb_mvir == mpeak
    
    vpeak = mb_vmax[peak]
    
    if len(vpeak == 1):
        Vpeakar[ii] = vpeak[0]
        apeak = mb_scale[peak]
        apeakar[ii] = apeak[0]
        
    else:
        apeak = mb_scale[peak]
        vpeakmax = max(vpeak)
        
        maxpeak = vpeak == vpeakmax
        
        apeak = apeak[maxpeak][0]
        
        Vpeakar[ii] = vpeakmax
        apeakar[ii] = apeak
        
    if pbar is not None:
        pbar.update(ii)

if pbar is not None:
    pbar.finish()
    
if 'Mpeak_TotalMass' in cat.keys():       
    ompeak = cat['Mpeak_TotalMass'][:]
    if 'Vpeak_TotalMass' in cat.keys(): 
        ovpeak = cat['Vpeak_TotalMass'][:]
    else:
        ovpeak = None
    if 'apeak_TotalMass' in cat.keys():
        oapeak = cat['apeak_TotalMass'][:]
    else:
        oapeak = None
    
    if 'Mpeak_TotalMassBackup' in cat.keys():
        backup = raw_input("Going to delete three (backup) datasets from the file; do you want a backup of the file? (Y/n) ")
        if backup != "y" and backup != "Y":
            print "Not creating a backup."
        else:
            from time import time
            from shutil import copyfile
            outname = fname+".{0}".format(int(time()))
            print "Creating backup as "+outname
            copyfile(fname,outname)
    
    print "Creating TotalMassBackup datasets."
    cat.create_dataset('Mpeak_TotalMassBackup',data=ompeak)
    if ovpeak is not None:  cat.create_dataset('Vpeak_TotalMassBackup',data=ovpeak)
    if oapeak is not None:  cat.create_dataset('apeak_TotalMassBackup',data=oapeak)

    del cat['Mpeak_TotalMass']
    if 'Vpeak_TotalMass' in cat.keys():  del cat['Vpeak_TotalMass']
    if 'apeak_TotalMass' in cat.keys():  del cat['apeak_TotalMass']
        
cat.create_dataset('Mpeak_TotalMass',data=Mpeakar)
cat.create_dataset('Vpeak_TotalMass',data=Vpeakar)
cat.create_dataset('apeak_TotalMass',data=apeakar)


f.close()

print "Done!\n"
