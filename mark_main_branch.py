from shutil import copyfile
import sys
from time import time

import h5py
from numpy import (empty, unique, sort, where,
                   sum, max, argmax, in1d,
                   zeros, array)

# given my sort of merger tree file, we grab all the trees from the catalog
# (to save time), then loop through those and for each we look at z = 0, pick
# the branch that has the most cumulative mass going all the way up, and so on
# and so on then we save that in a new dataset in each tree

# take two of this script.  let's try to make it faster by starting at the top
# (i.e. the earliest scale factor in the tree), then adding downwards

# that means that I should probably write a function that just takes care of a
# single tree, then call that function for each tree


def make_main_branch(scale, ids, mvir, desc_id, dopbar=False):
    """
    returns an array of length equal to each of the inputs, where a 1 signifies
    that it's the main branch of the tree, as defined by the total mass above
    it.  all inputs must be the same length.

    :param scale:
        an array of scale factors for the halos in the tree
    :param ids:
        an array of the ids of the halos in the tree
    :param mvir:
        an array of the virial masses of the halos in the tree
    :param desc_id:
        an array of the descendant ids of the halos in the tree (i.e. it finds
        the progenitors of each halo

    :returns:
        mmp:  an array that is 1 if the corresponding halo is the main branch,
              and a 0 if it isn't
        mabove: an array of the mass above (and including) that halo in the
                tree

    """
    uniq_a = unique(scale)
    uniq_a = sort(uniq_a)

    if scale.shape[0] > 1e4:
        print("Have to touch all {0} halos in the tree once."
              .format(scale.shape[0]))

    mabove = empty(len(scale))
    mmp = empty(len(scale))

    mabove[:] = 0   # includes the mass at that level
    mmp[:] = 0

    for a in uniq_a:
        ata = scale == a

        stepids = ids[ata]
        stepmvir = mvir[ata]

        for ii in range(len(stepids)):
            thisid = stepids[ii]
            thismass = stepmvir[ii]

            index = where(thisid == ids)[0]
            try:
                assert index.shape[0] == 1
            except AssertionError:
                print("Got two identical indices.")
                print(thisid)
                print(ids[index])
                print(index)
                print(ids[index])
                print(stepids[ii+1])
                sys.exit(1337)
            mabove[index] = thismass + sum(mabove[desc_id == thisid])
            # this works because i've sorted uniq_a,
            # so i'm by design going down the tree

    uniq_a = uniq_a[::-1]       # now it goes from largest to smallest

    for ii in range(len(uniq_a)):
        ata = scale == uniq_a[ii]
        if ii == 0:
            assert where(ata)[0].shape[0] == 1
            mmp[ata] = 1
            continue

        # now, if i'm not at the head of the tree, then I need to find the
        # progenitors of the current mmp and mark the max of those as the mmp

        prevmmpid = ids[(mmp == 1) & (scale == uniq_a[ii-1])]
        assert prevmmpid.shape[0] == 1

        progen_ids = ids[desc_id == prevmmpid]

        if progen_ids.shape[0] == 0:
            # then the previous most massive progenitor's branch didn't extend
            # to the earliest time that's fine--it just means that we break out
            # of this loop
            break

        progen_mabove = mabove[in1d(ids, progen_ids)]

        assert progen_ids.shape[0] == progen_mabove.shape[0]

        thismmp_id = progen_ids[argmax(progen_mabove)]

        index = where(ids == thismmp_id)[0]

        assert index.shape[0] == 1

        mmp[index] = 1

    return mmp, mabove


# this needs to be a recursive function with a for loop, so it will be SLOW:
# return (mass in this level) + sum[mass in each branch above]
def sum_mass_above(id, allids, desc_ids, mvir):
    """
    The main function here, which I'm going to call recursively. It'll return
    the mass at the current level, plus (recursively) the mass of all the
    levels above it.

    So, I'll have to run this function at each timestep essentially, pick the
    halo whose branch has the most mass, mark that as the main prog, and then
    move on to the same with all the halos at that are possible progens of
    that main prog

    so, in other words, this is gonna be fucking slow as fuck

    :param tree:
        an HDF5 group containing the merger tree that I'm analyzing
    :param desc_id:
        the ID of the halo that I'm looking for progenitor branches
        of (i.e. the halo that I want all the progenitors of)

    :returns:
        in the base case (top of the tree), returns the mass of the halo at the
        top of the tree

        in the recursive case, returns the mass of all the halos at the current
        level in the tree, plus the mass of each branch extending out of the
        current level
    """

    loc = where(id == allids)[0]

    assert loc.shape[0] == 1

    thismass = mvir[loc][0]

    progens = where(desc_ids == id)[0]

    if isinstance(progens, type(array([0, 1, 2]))):
        progens = array([progens])
    progen_ids = allids[progens]

    if progen_ids.shape[0] == 0:
        return thismass

    else:
        return (thismass
                + sum([sum_mass_above(pid, allids, desc_ids, mvir)
                       for pid in progen_ids]))


if __name__ == "__main__":
    usage = ("usage:  python {0} <input HDF5 file> [Vmax cut = 7 km/s] "
             "[trees group = 'RockstarMergerTrees']".format(sys.argv[0]))

    if len(sys.argv) < 2:
        print(usage)
        sys.exit(1337)

    vcut = 7
    tname = 'RockstarMergerTrees'

    fname = sys.argv[1]
    if len(sys.argv) > 2:
        vcut = float(sys.argv[2])
    if len(sys.argv) > 3:
        tname = sys.argv[3]

    print("Opening {0}".format(fname))

    f = h5py.File(fname)

    mt = f[tname]

    snaps = [k for k in list(f.keys()) if k.startswith('Snapshot')]
    snaps.sort()
    snap = f[snaps[-1]]
    cat = snap['HaloCatalog_RockstarMergerTree']

    tnums = cat['TreeNumber'][:]

    print("Skipping trees with Vmax < {0} km/s in the last timestep"
          .format(vcut))

    copied = False
    ndel = 0

    for ii in range(len(tnums)):
        tree = mt['Tree_{0}'.format(tnums[ii])]

        if (('TotalMass_mmp' in tree.keys() or 'TotalMassAbove' in tree.keys())
                and not copied):
            outname = fname+".{0}".format(int(time()))
            print("Creating a backup of the file named {0} because I'm "
                  "going to overwrite some datasets".format(outname))
            copyfile(fname, outname)
            copied = True

        if 'TotalMass_mmp' in list(tree.keys()):
            del tree['TotalMass_mmp']
            ndel = ndel + 1

        if 'TotalMassAbove' in list(tree.keys()):
            del tree['TotalMassAbove']

        scale = tree['scale'][:]
        allids = tree['id'][:]
        desc_ids = tree['desc_id'][:]
        mvir = tree['mvir'][:]
        vmax = tree['vmax'][:]

        loc = where(scale == max(scale))[0]

        assert loc.shape[0] == 1

        if vmax[loc] < vcut:
            # then stick in a placeholder dataset filled with minus 2
            mymmp = zeros(len(scale), dtype='int')
            mymmp[:] = -2
            mabove = zeros(len(scale), dtype='float')
            mabove[:] = -2
            tree.create_dataset('TotalMass_mmp', data=mymmp)
            tree.create_dataset('TotalMassAbove', data=mabove)
            continue

        mymmp, mabove = make_main_branch(scale, allids, mvir,
                                         desc_ids, dopbar=True)

        tree.create_dataset('TotalMass_mmp', data=mymmp)
        tree.create_dataset('TotalMassAbove', data=mabove)

    if ndel == 0:
        print("Done!\n")
    else:
        print("Overwrote an existing 'TotalMass_mmp' dataset in {0} of {1} "
              "trees\nBackup can be found as {2}, but finished!\n"
              .format(ndel, tnums.shape[0], outname))
