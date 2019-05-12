import re
from subprocess import Popen, PIPE
import sys

import h5py
import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:
    have_pbar = False
else:
    have_pbar = True


def get_col_names(fname):
    """Parses the header to get column names and parameters."""
    with open(fname) as f:
        cols = f.readline().strip("#\n").lower()
        cols = (re.sub(r'\(\d+\)', '', cols)
                  .replace('/', '_to_')
                  .split())
        return cols


def get_cosmo(fname):
    """Gets the cosmological parameters from the header."""
    grep_cosmo = Popen(['grep', '^#Omega_', str(fname)], stdout=PIPE)
    grep_box = Popen(['grep', '^#Full box', str(fname)], stdout=PIPE)
    cosmo_str = (grep_cosmo.communicate()[0]
                           .decode("utf-8")
                           .strip("#\n")
                           .split("; "))
    box_str = (grep_box.communicate()[0]
                       .decode("utf-8")
                       .strip("#\n")
                       .split(" = "))
    cosmo = {i.split(' = ')[0]: float(i.split(' = ')[1]) for i in cosmo_str}
    cosmo['Box_size_Mpc/h'] = float(box_str[1].split()[0])
    return cosmo


def get_tree_nums(fname):
    """Uses bash commands to grab the tree number from each row."""
    grep_cmd = Popen(["grep", "#tree", str(fname)], stdout=PIPE)
    awk_cmd = Popen(["awk", "{print $2}"], stdin=grep_cmd.stdout, stdout=PIPE)
    out, err = awk_cmd.communicate()
    out = out.decode("utf-8")  # 'out' is a byte-string in Python 3
    out = out.splitlines()
    return np.array(out, dtype=int)


def main_trees(df, prev_id=None):
    """Assign halos to the main branch.

    Takes in a full merger tree and determines the most massive progenitor
    (MMP) at each timestep.

    Parameters
    ----------
    df : DataFrame
        Full merger tree, scale descending
    prev_id : int, optional
        `id` of the most recent halo. The default will start at the
        last descendant.

    Returns
    -------
    mmp : ndarray, int
        An array containing `id` of the most massive progenitor at each
        timestep.

    Raises
    ------
    RuntimeError
        If `scale` is not descending.
    """
    if not df.scale.is_monotonic_decreasing:
        raise RuntimeError("`df.scale` is not descending.")

    a_uniq = df.scale.unique()
    mmp = np.zeros_like(a_uniq, dtype=int)
    if prev_id is None:
        prev_id = df[df.scale == a_uniq[0]].id.values[0]
    for i, a in enumerate(a_uniq):
        if i == 0 and a == 1:
            mmp[i] = df[df.id == prev_id].index[0]
            continue
        msk = (df.scale == a) & (df.desc_id == prev_id)
        if msk.sum() == 0:
            # We have reached the end of the branch.
            break
        mmp[i] = df[msk].mvir.idxmax()
        prev_id = df.loc[mmp[i]].id
    return mmp[mmp > 0]


def main_trees_quick(df):
    """Lazily mark main branches as most massive halo at each timestep.

    Groups the halos in the catalog by tree and scale. Then designates
    the most massive halo within each group as the most massive progenitor.

    Generally, this function will get most halos correctly marked as main
    progenitor halos. If a complete and accurate merger tree is required,
    combine this with `verify_main_branches` to ensure accuracy.

    The errors occur if there are two halos, both of maximal mass, at the
    same timestep. This cannot be avoided while remaining quick due to the
    equality used within each 'tree-scale' group.

    Parameters
    ----------
    df : DataFrame
        Complete catalog of merger trees

    Returns
    -------
    mmp : ndarray, bool
        True if a halo is the most massive within its tree and scale group

    """
    groups = df.groupby(['tree', 'scale'])
    try:
        mmp = groups.mvir.transform(np.max) == df.mvir
    except (KeyError, AttributeError):
        mmp = groups.Mvir.transform(np.max) == df.Mvir
    return mmp


def correct_mmp(df, prev_id):
    """Identifies the progenitor halos of the main branch.

    Starting at the prev_id, walk through the merger tree
    matching halos at the current timestep (a) with the
    most massive progenitor at a previous timestep (prev_id).

    Parameters
    ----------
    df : DataFrame
        Full merger tree, scale descending
    prev_id : int
        id of the last correct halo in the main branch

    Returns
    -------
    mmp : ndarray
        An array of the corrected most massive progenitor IDs

    Raises
    ------
    RuntimeError
        If `scale` is not descending.
    """
    if not df.scale.is_monotonic_decreasing:
        raise RuntimeError("`df.scale` is not descending.")

    prev_a = df.scale.values[(df.id == prev_id)][0]
    a_uniq = df.scale[(df.scale <= prev_a)].unique()
    mmp = np.zeros_like(a_uniq, dtype=int)
    for i, a in enumerate(a_uniq):
        # Place index of previous halo in mmp to correct for duplicate
        # halos occurring at the same timestep
        if i == 0:
            mmp[i] = df.loc[df.id == prev_id].index.values[0]
            continue

        msk = (df.scale == a) & (df.desc_id == prev_id)
        if msk.sum() == 0:
            # We have reached the end of this main branch.
            mmp = mmp[mmp > 0]
            break
        mmp[i] = df[msk].mvir.idxmax()
        prev_id = df.loc[mmp[i]].id
    return mmp


def check_main_branches(df):
    """Checks that the main branch is comprised of progenitors.

    Examines the desc_id of all halos in the main branch and compares
    that to the next halo's halo_id. The first halo (descendant at a = 1)
    is ignored because it does not have any descendants. Similarly, the
    last halo is ignored because it does not have any progenitors.

    Parameters
    ----------
    df : DataFrame
        Main branch of a merger tree, scale descending

    Returns
    -------
    last_id : int or None
        id of the last correct halo in the main branch. None if all halos are
        progenitors of the first halo

    Raises
    ------
    RuntimeError
        If `scale` is not descending.
    """
    if not df.scale.is_monotonic_decreasing:
        raise RuntimeError("`df.scale` is not descending.")

    # First halo as no descendants and the last halo has no progenitor.
    desc_ids = df.desc_id.values[1:]
    halo_ids = df.id.values[:-1]
    # The desc_id of the last halo should be the next halo_id
    if np.all(desc_ids == halo_ids):
        return None
    else:
        # Return ID of the last correct halo
        return df.id.values[np.argmin(desc_ids == halo_ids)]


def verify_main_branches(df):
    """Verifies the accuracy of the main branches.

    Checks that all of the halos within the main branch are progenitors of
    each other. If they are not, identify the last progenitor and update the
    main branch with the correct halos.

    Parameters
    ----------
    df : DataFrame
        Complete merger tree, scale descending

    Returns
    -------
    df : DataFrame
        Updated merger tree catalog with the correct main branches marked
    """
    tnums = df.tree.unique()
    if have_pbar:
        tnums = tqdm(tnums)

    for tn in tnums:
        tree = df.loc[df.tree == tn]
        last_id = check_main_branches(tree.loc[tree.TotalMass_mmp])
        if last_id is not None:
            last_a = tree.loc[tree.id == last_id].scale.values[0]
            new_mmp = correct_mmp(tree, last_id)
            msk = (df.tree == tn) & (df.scale <= last_a)
            df.loc[msk, 'TotalMass_mmp'] = np.isin(df.loc[msk].index, new_mmp)
    return df


def make_tree_col(df, tnums):
    """Applies tree numbers to corresponding rows (tree members)."""
    tree_inds = np.where(df.scale.values == 1)[0]
    tree_inds = np.append(tree_inds, df.shape[0])
    return np.repeat(tnums, np.diff(tree_inds))


# If the output looks weird, consider replacing .iloc[1:] with .dropna()
def read_data(fname, cols):
    """Loads consistent-trees file into pandas dataframe."""
    df = (pd.read_csv(fname, header=None, sep=r"\s+", comment="#",
                      names=cols, dtype=np.float64)
            .iloc[1:])  # First line is the total number of trees
    # Could reset_index, but we don't shuffle the DataFrame
    return df


def analyze_trees(fname, only_mb=False, slow_mb=False):
    """Processes a complete merger tree catalog."""
    cols = get_col_names(fname)
    tnums = get_tree_nums(fname)
    df = read_data(fname, cols)
    df['tree'] = make_tree_col(df, tnums)
    if slow_mb:
        df['TotalMass_mmp'] = False
        if have_pbar:
            tnums = tqdm(tnums)
        for tn in tnums:
            mmps = main_trees(df.loc[df.tree == tn])
            mmps = np.isin(df.loc[df.tree == tn].index, mmps)
            df.loc[df.tree == tn, 'TotalMass_mmp'] = mmps
    else:
        df['TotalMass_mmp'] = main_trees_quick(df)

    df = verify_main_branches(df)
    if only_mb:
        return df.loc[df.TotalMass_mmp == 1]
    return df


def save_to_hdf5(fname, df, cosmo={}, tname="RockstarMergerTrees", min_vmax=0):
    """Saves the merger trees mostly following IRATE practices.

    Parameters
    ----------
    fname : str
        Path to the output file
    df : DataFrame
        Merger tree catalog to be saved to `fname`
    cosmo : dict-like, optional
        Dictionary containing the simuation parameters used
    tname : str, optional
        HDF5 group name given to the merger trees inside of the file
    min_vmax : float, optional
        Minimum Vmax (km/s) for a tree to be written to the file
    """
    f = h5py.File(fname, 'a', libver='latest')
    colheads = df.columns.values
    treenums = df.loc[df.vmax >= min_vmax].tree.unique()
    if tname in f.keys():
        print("File already contains a group named {0}, so I can't save to it."
              " Exiting.".format(tname))
        sys.exit(1337)
    t = f.create_group(tname)
    if have_pbar:
        treenums = tqdm(treenums)
    for i, tnum in enumerate(treenums):
        tg = t.create_group('Tree_' + str(tnum))
        for j, col in enumerate(colheads):
            col_data = df.loc[(df.tree == tnum), col].values
            tg.create_dataset(col, data=col_data)
    head = f.create_group('Header')
    for param in cosmo:
        head.create_dataset(param, data=cosmo[param])
    f.close()


def add_z0_catalog(fname, df, snap=152, name="HaloCatalog_RockstarMergerTree"):
    """Adds z=0 halos from the merger trees to a separate catalog.

    Parameters
    ----------
    fname : str
        Path to the output file
    df : DataFrame
        Merger tree catalog to pull halos from
    snap : int
        The final snapshot number of the simulation
    name : str, optional
        HDF5 group name given to the merger trees inside of the file
    """
    f = h5py.File(fname, 'a', libver='latest')
    cat = f.create_group("Snapshot{:05d}/{}".format(snap, name))
    z0_cat = df.loc[(df.scale == 1)].drop('scale', axis=1)
    # Combine coordinate quantities into single columns
    z0_cat['center'] = list(z0_cat[['x', 'y', 'z']].values)
    z0_cat['velocity'] = list(z0_cat[['vx', 'vy', 'vz']].values)
    z0_cat['j'] = list(z0_cat[['jx', 'jy', 'jz']].values)
    z0_cat.drop(
        ['x', 'y', 'z', 'vx', 'vy', 'vz', 'jx', 'jy', 'jz'],
        axis=1,
        inplace=True,
    )
    colheads = z0_cat.columns.values
    for col in colheads:
        try:
            cat.create_dataset(col, data=z0_cat[col].values)
        except TypeError:  # No HDF5 datatype for object (e.g. 'center' above)
            cat.create_dataset(col, data=np.vstack(z0_cat[col].values))
    f.close()


if __name__ == '__main__':
    import argparse
    from timeit import default_timer as time

    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Path to the Consistent-Trees file")
    parser.add_argument("outname", help="Path to the output file")
    parser.add_argument(
        "--main-branch",
        action="store_true",
        dest="mb",
        help="Return only the main branches",
    )
    parser.add_argument(
        "-t",
        "--timing",
        action="store_true",
        help="Report timings of each function call",
    )
    parser.add_argument(
        "--slow",
        action="store_true",
        help="Slowly define the main branches (ensures accuracy)"
    )
    args = parser.parse_args()

    start_time = time()
    mt = analyze_trees(args.file, args.mb, args.slow)
    cosmo = get_cosmo(args.file)
    save_to_hdf5(args.outname, mt, cosmo)
    add_z0_catalog(args.outname, mt)
    total = time() - start_time
    if args.timing:
        print("Total time: {:.0f}min {:.1f}s".format(total/60, total % 60))
