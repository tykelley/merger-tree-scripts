import os
from subprocess import Popen, PIPE
import sys

import h5py
import numpy as np


if __name__ == "__main__":
    # I want to open a file that has a Rockstar merger tree in it,
    # extract the z = 0 data, save that in a new halo catalog,
    # and also save Vmax_max as a derived quantity in that new catalog

    if len(sys.argv) < 4:
        print("usage:  python {0} <HDF5 file> <snap for catalog>"
              " <tree file> [cat_name = 'HaloCatalog_RockstarMergerTree']"
              " [trees group = 'RockstarMergerTrees']"
              " ['scale to save = 1']".format(sys.argv[0]))
        sys.exit(1337)

    fname = sys.argv[1]
    outsnap = int(sys.argv[2])
    treefile = sys.argv[3]
    if len(sys.argv) >= 5:
        outcatname = sys.argv[4]
        if not outcatname.startswith('HaloCatalog'):
            print("Prepending 'HaloCatalog_' to {0} to "
                  "conform with IRATE specs.".format(outcatname))
            outcatname = 'HaloCatalog_'+outcatname
    else:
        outcatname = 'HaloCatalog_RockstarMergerTree'

    if len(sys.argv) >= 6:
        tname = sys.argv[5]
    else:
        tname = 'RockstarMergerTrees'

    if len(sys.argv) >= 7:
        atosave = float(sys.argv[6])
    else:
        atosave = 1

    print("Saving data at a = {0}".format(atosave))

    # if I'm going to require this, then I should get the trees from it.
    # Probably faster than getting it from the HDF5 file.
    f = open(treefile, 'r')
    for line in f:
        if not line.startswith('#'):
            numtrees = int(line)
            break
    f.close()

    print("Using grep to find the tree numbers in {0}".format(treefile))
    grep_cmd = "grep '#tree' %s | awk '{print $2}'" % (treefile)
    shell_out = Popen(grep_cmd, stdout=PIPE, shell=True).communicate()[0]
    tnums = np.array(shell_out.strip().split("\n"), dtype=int)
    assert len(tnums) == numtrees
    print("Expecting {0} trees in the file.".format(numtrees))

    if not os.path.isfile(fname):
        print("{0} doesn't exist; exiting.".format(fname))
        sys.exit(1337)

    print("Opening {0}".format(fname))

    f = h5py.File(fname, 'a')
    mt = f[tname]

    # Check that I have a place to save the data
    if "Snapshot{0:05}".format(outsnap) in list(f.keys()):
        if outcatname in f['Snapshot{0:05}'.format(outsnap)]:
            print("Existing output catalog (Snapshot{0:05}/{1}) already "
                  "exists.".format(outsnap, outcatname))
            sys.exit(1337)

    tnames = ['Tree_{0}'.format(tnums[i]) for i in range(numtrees)]
    t0 = tnames[0]

    colheads = list(mt[t0].keys())
    arrays = [np.empty(numtrees, dtype=mt[t0][colheads[i]].dtype)
              for i in range(len(colheads))]
    for i in range(len(arrays)):
        arrays[i][:] = -2  # Cause nothing else is using -2

    colheads = np.array(colheads)
    treenums = np.empty(numtrees, dtype=int)
    treenums[:] = -2
    nskip = 0
    print("Beginning loop over trees.")

    for i in range(numtrees):
        tree = mt[tnames[i]]
        treenums[i] = tree.name.split('/')[-1].split('_')[-1]  # do this first

        scale = tree['scale'][:]
        if len(np.where(scale == atosave)[0]) != 1:
            print("Skipping a halo that has no z = {0} information."
                  .format(1./atosave - 1))
            nskip = nskip + 1
            continue
        # z0 is really a misnomer, because it's for the redshift specified
        z0index = np.where(scale == atosave)[0][0]
        for j in range(len(list(tree.keys()))):  # loop over the columns
            # Get the value of the column, in the row where z = 0
            # Save that value in the array
            z0val = tree[colheads[j]][z0index]
            arrays[j][i] = z0val
        if i % 1000 == 0:
            print("Done with tree {0} of {1}".format(i, numtrees))

    print("Skipped {0} halo(s) because they didn't have a redshift {1} entry."
          .format(nskip, (1./atosave) - 1))

    # Now I want to combine together x,y,z into center,
    # vx,vy,vz into velocity, and Jx,Jy,Jz into J
    print("Combining x,y,z into Center; vx,vy,vz into Velocity; "
          "and Jx,Jy,Jz into J")
    x = arrays[np.where(colheads == 'x')[0][0]]
    y = arrays[np.where(colheads == 'y')[0][0]]
    z = arrays[np.where(colheads == 'z')[0][0]]

    vx = arrays[np.where(colheads == 'vx')[0][0]]
    vy = arrays[np.where(colheads == 'vy')[0][0]]
    vz = arrays[np.where(colheads == 'vz')[0][0]]

    jx = arrays[np.where(colheads == 'Jx')[0][0]]
    jy = arrays[np.where(colheads == 'Jy')[0][0]]
    jz = arrays[np.where(colheads == 'Jz')[0][0]]

    cen = np.array(zip(x, y, z))
    vel = np.array(zip(vx, vy, vz))
    J = np.array(zip(jx, jy, jz))

    # Pull out everything that's not x,y,z or vx,vy,vz or jx,jy,jz
    myarrays = [arrays[i] for i in range(len(arrays))
                if (not colheads[i] == 'x'
                    and not colheads[i] == 'y'
                    and not colheads[i] == 'z'
                    and not colheads[i] == 'vx'
                    and not colheads[i] == 'vy'
                    and not colheads[i] == 'vz'
                    and not colheads[i] == 'Jx'
                    and not colheads[i] == 'Jy'
                    and not colheads[i] == 'Jz')]
    mycolheads = [colheads[i] for i in range(len(arrays))
                  if (not colheads[i] == 'x'
                      and not colheads[i] == 'y'
                      and not colheads[i] == 'z'
                      and not colheads[i] == 'vx'
                      and not colheads[i] == 'vy'
                      and not colheads[i] == 'vz'
                      and not colheads[i] == 'Jx'
                      and not colheads[i] == 'Jy'
                      and not colheads[i] == 'Jz')]

    # Then add the arrays that I combined/made on my own
    myarrays.append(cen)
    mycolheads.append('Center')
    myarrays.append(vel)
    mycolheads.append('Velocity')
    myarrays.append(J)
    mycolheads.append('J')

    print("Capitalizing the first letter of each dataset name.")
    for i in range(len(mycolheads)):
        mycolheads[i] = mycolheads[i].capitalize()

    # create a group to hold the z = 0 data (after reading it)
    if "Snapshot{0:05}".format(outsnap) in list(f.keys()):
        print("Creating {0} under group Snapshot{1:05} to save z = 0 data"
              .format(outcatname, outsnap))
        outcat = f['Snapshot{0:05}'.format(outsnap)].create_group(outcatname)
    else:
        print("Creating group Snapshot{0:05}/{1} to save z = 0 data."
              .format(outsnap, outcatname))
        f.create_group("Snapshot{0:05}".format(outsnap))
        outcat = f['Snapshot{0:05}'.format(outsnap)].create_group(outcatname)

    # Now I need to save the information in outcat
    print("Saving {1} columns of z = 0 information to {0}"
          .format(outcatname, len(mycolheads)))
    for i in range(len(mycolheads)):
        outcat.create_dataset(mycolheads[i], data=myarrays[i])
    outcat.create_dataset('TreeNumber', data=treenums)

    print("Finished writing to {}".format(outcatname))
    f.close()
