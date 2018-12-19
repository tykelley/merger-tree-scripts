import re
import sys

import h5py
import numpy as np


def read_tree(fname):
    print("Opening {0}".format(fname))
    thistree = []
    comments = []
    treenums = []
    trees = []
    ex = 0
    startedtrees = False
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                colheads = line.strip('#\n')
                colheads = (re.sub(r'\(\d+\)', '', colheads)
                            .replace('/', '_to_')
                            .split())
            elif line.startswith('#tree'):
                startedtrees = True
                treenums.append(int(line.split()[-1]))
                if not thistree:
                    continue
                else:
                    trees.append(np.vstack(thistree).T)
                    # print("Finished tree number {0}".format(treenums[-2]))
                    thistree = []
            elif startedtrees:
                thistree.append(np.asfarray(line.strip().split()))
            elif line.startswith('#'):
                comments.append(line.strip('#'))
            else:
                ntrees = int(line)
                print("Expecting {0} trees in the file".format(ntrees))
                ex = ex + 1
                if ex > 2:
                    sys.exit(1)
    print("Reading final tree and saving it...")
    trees.append(np.vstack(thistree).T)
    assert len(trees) == ntrees

    return colheads, trees, treenums, comments


def save_to_hdf5(outname, colheads, trees, treenums, comments, tname):
    print("Saving data to {0} in file {1}".format(tname, outname))
    f = h5py.File(outname, 'a', libver='latest')
    if tname in f.keys():
        print("File already contains a group named {0}, so I can't save to it."
              " Exiting.".format(tname))
        sys.exit(1337)
    t = f.create_group(tname)
    for c in comments:
        t.attrs[c] = -1

    for i, tnum in enumerate(treenums):
        tre = t.create_group('Tree_' + str(tnum))
        for j, col in enumerate(colheads):
            dset = tre.create_dataset(col, trees[i][j].shape, 'f8')
            dset = trees[i][j]  # noqa: F841
    f.close()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage:  python {0} <inname> <outname> "
              "[tree name = 'RockstarMergerTrees']".format(sys.argv[0]))
        sys.exit(1337)
    [colheads, trees, treenums, comments] = read_tree(sys.argv[1])
    outname = sys.argv[2]
    if not outname.endswith('.hdf5'):
        outname = outname + '.hdf5'
    if len(sys.argv) > 3:
        tname = sys.argv[3]
    else:
        tname = 'RockstarMergerTrees'
    save_to_hdf5(outname, colheads, trees, treenums, comments, tname)
