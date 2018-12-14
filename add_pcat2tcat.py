import sys
import h5py
from numpy import empty_like, empty

usage = ("usage:  python {0} <input HDF5 file> [snapshot = last] "
         "[pcat name = 'HaloCatalog_ParallelRockstar'] "
         "[tcat name = 'HaloCatalog_RockstarMergerTree']"
         .format(sys.argv[0].split('/')[-1]))

if len(sys.argv) < 2:
    print(usage)
    sys.exit(1337)

pcatname = 'HaloCatalog_Rockstar'
tcatname = 'HaloCatalog_RockstarMergerTree'
snap = None

fname = sys.argv[1]
if len(sys.argv) > 2:
    snap = int(sys.argv[2])
if len(sys.argv) > 3:
    pcatname = sys.argv[3]
    if not pcatname.startswith('HaloCatalog_'):
        pcatname = 'HaloCatalog_' + pcatname
if len(sys.argv) > 4:
    tcatname = sys.argv[4]
    if not tcatname.startswith('HaloCatalog_'):
        tcatname = 'HaloCatalog_' + tcatname

print("\nOpening "+fname)
f = h5py.File(fname)

if snap is None:
    snaps = [k for k in list(f.keys()) if k.startswith('Snapshot')]
    snaps.sort()
    snap = snaps[-1]
else:
    snap = 'Snapshot{0:05}'.format(snap)
s = f[snap]
pcat = s[pcatname]
tcat = s[tcatname]

print("Adding data from {0} to {1}".format(pcat.name, tcat.name))

tvmax = tcat['Vmax'][:]
pvmax = pcat['Vmax'][:]


pkeys = list(pcat.keys())
tkeys = list(tcat.keys())

if 'Orig_halo_id' in list(tcat.keys()):
    print("Matching 'Orig_halo_id' with 'ID'")
    orig_id = tcat['Orig_halo_id'][:]
elif 'Orighaloid' in list(tcat.keys()):
    print("Matching 'Orighaloid' with 'ID'")
    orig_id = tcat['Orighaloid'][:]
else:
    print("No original halo ID dataset found; exiting.")
    sys.exit(1337)
pid = pcat['ID'][:]

assert pid.shape[0] >= orig_id.shape[0]

# Datasets that I'm NOT going to add to tcat
skip = ['Center', 'ID', 'J', 'Mvir', 'Rvir',
        'Spin', 'Velocity', 'Vmax', 'Vrms']

toadd = [k for k in pkeys if k not in skip and k not in list(tcat.keys())]

if len(toadd) == 0:
    print("No datasets to add!\nExiting.")
    sys.exit(1337)

print("Adding the following datasets:")
for k in toadd:
    print("\t"+k)

dsets = [pcat[k][:] for k in toadd]
# empties = [empty_like(k) for k in dsets]
empties = []
for ii in range(len(dsets)):
    shape = dsets[ii].shape
    if len(shape) == 1:
        empties.append(empty_like(tvmax))
    else:
        empties.append(empty((tvmax.shape[0], shape[1])))

for ii in range(len(empties)):
        empties[ii].fill(-1)

nomatch = 0
multimatch = 0

for ii in range(len(tvmax)):
    myid = orig_id[ii]
    loc = pid == myid

    mypvmax = pvmax[loc]
    mytvmax = tvmax[ii]

    if mypvmax.shape[0] == 0:
        nomatch = nomatch + 1
        continue
    elif mypvmax.shape[0] > 2:
        multimatch = multimatch + 1
        continue

    mypvmax = mypvmax[0]

    print(mypvmax)
    print(mytvmax)
    print(abs(mypvmax - mytvmax)/mytvmax)
    try:
        assert abs(mypvmax - mytvmax)/mytvmax < 2e-2
    except AssertionError:
        print("Failed assertion; exiting.")
        print(abs(mypvmax - mytvmax))
        print(abs(mypvmax - mytvmax)/mytvmax)
        print(ii)
        print(myid)
        f.close()
        sys.exit(1337)

    for jj in range(len(empties)):
        val = dsets[jj][loc][0]
        empties[jj][ii] = val

if nomatch > 0:
    print("There were {0} halos for which I couldn't find a match "
          "-- saving -1 as the values".format(nomatch))
if multimatch > 0:
    print("There were {0} halos for which I found multiple matches -- saving "
          "-1 as the values for those".format(nomatch))

if nomatch > 0 or multimatch > 0:
    from shutil import copyfile
    from time import time
    outname = fname + str(int(time()))
    print("Saving a backup file as {0}, since I found either no matches or "
          "multiple matches in some cases.".format(outname))
    copyfile(fname, outname)

print("Saving those datasets to {0}".format(tcat.name))

for ii in range(len(empties)):
    tcat.create_dataset(toadd[ii], data=empties[ii])

f.close()
print("Done!")
