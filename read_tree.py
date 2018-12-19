import re
from subprocess import Popen, PIPE

import numpy as np
import pandas as pd


def get_col_names(fname):
    with open(fname) as f:
        cols = f.readline().strip("#\n")
        cols = (re.sub(r'\(\d+\)', '', cols)
                  .replace('/', '_to_')
                  .split())
        return cols


def get_tree_nums(fname):
    grep_cmd = Popen(["grep", "#tree", str(fname)], stdout=PIPE)
    awk_cmd = Popen(["awk", "{print $2}"], stdin=grep_cmd.stdout, stdout=PIPE)
    out, err = awk_cmd.communicate()
    out = out.decode("utf-8")  # 'out' is a byte-string in Python 3
    out = out.splitlines()
    return np.array(out, dtype=int)


def main_trees(df):
    groups = df.groupby(['tree', 'scale'], sort=False)
    try:
        idx = groups.mvir.transform(max) == df.mvir
    except (AttributeError, KeyError):
        idx = groups.Mvir.transform(max) == df.Mvir
    return df[idx]


def make_tree_col(df, tnums):
    tree_inds = np.where(df.scale.values == 1)[0]
    tree_inds = np.append(tree_inds, df.shape[0])
    return np.repeat(tnums, np.diff(tree_inds))


# If the output looks weird, consider replacing .iloc[1:] with .dropna()
def read_data(fname, cols):
    """Loads consistent-trees file into pandas dataframe"""
    df = (pd.read_csv(fname, header=None, sep=r"\s+", comment="#",
                      names=cols, dtype=np.float64)
            .iloc[1:])  # First line is the number of trees
    # Could reset_index, but won't use the index anyways
    return df


def analyze_trees(fname, return_full_tree=False):
    cols = get_col_names(fname)
    tnums = get_tree_nums(fname)
    df = read_data(fname, cols)
    df['tree'] = make_tree_col(df, tnums)
    mt = main_trees(df)
    if return_full_tree:
        return df, mt
    return mt


if __name__ == '__main__':
    import argparse
    from timeit import default_timer as time

    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="path to the Consistent-Trees file")
    parser.add_argument("-a", "--all", action="store_true",
                        help="return all the trees and the main branches")
    parser.add_argument("-t", "--timing", action="store_true",
                        help="report timings of each function call")
    args = parser.parse_args()

    start_time = time()
    mt = analyze_trees(args.file, args.all)
    run_time = time() - start_time
    print("Total time: {:.0f}min {:.1f}s".format(run_time/60, run_time % 60))
