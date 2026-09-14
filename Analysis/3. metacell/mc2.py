#!/usr/bin/env python
"""
metacell2 ("metacells", a.k.a. mc2) one-pass wrapper
A wrapper with a minimal functioning run of metacell2.
Takes as input:
 - .h5ad anndata object of cells x features,
 - number of target_UMIs_metacell
 - (optional): list of genes to be entirely removed from the dataset
 - (optional): list of genes to be marked as lateral genes
 - output file
 
After prepping input, it runs the mc2 divide_and_conquer wrapper.

Returns as output: 
 - .h5ad of cells 
 - .csv with the cell barcode / metacell2 ID.
"""
# standard libraries
import sys                       # For exit function
import os                        # For filesystem operations
import shutil                    # for filesystem operations
import anndata as ad             # For reading/writing AnnData files
import argparse                  # For parsing arguments

# external libraries
import metacells as mc           # The Metacells package
import numpy as np               # For array/matrix operations
import pandas as pd              # For data frames
import scipy.sparse as sp        # For sparse matrices

def main():
    """
    main entry point
    """
    
    # parsing arguments
    parser=argparse.ArgumentParser(description="Metacell2 one-pass run")
    parser.add_argument("-i", help="input h5ad file")
    parser.add_argument("-N", help="target UMIs per metacell")
    parser.add_argument("-b", help="Blacklisted genes to be removed")
    parser.add_argument("-l", help="Lateral genes to be ignored")
    parser.add_argument("-r", help="Min counts for rare gene modules")
    parser.add_argument("-f", help="Max fraction of cells expressing rare gene modules")
    # parser.add_argument("-p", help="Max parallel piles, optional")
    parser.add_argument("-o", help="output prefix. If e.g. 'path/XYZ', it will generate 2 files: path/XYZ.mc2.h5ad and path/XYZ.mc2.csv .")
    # parser.add_argument("-v", help="Verbose. Print messages if set to 'True'.")
    args = parser.parse_args()

    # reading the data
    print("mc2 | Loading h5ad ...")
    cells = ad.read_h5ad(args.i)
    cells.X = cells.X.astype("float32") # see https://github.com/tanaylab/metacells/issues/55
    print("mc2 | " + str(cells.shape[0]) + " cells; " + str(cells.shape[1]) + " genes.")

    # mc toplevel, is this needed?
    mc.ut.top_level(cells)
    mc.ut.set_name(cells, "mc2")

    # remove blacklisted genes if present
    if args.b is not None :
        print("mc2 | Blacklisting genes ...")
        list_bl = [i.replace("\n","") for i in open(args.b, "r")]
        filt_bl = [i for i in cells.var_names if i not in list_bl]
        cells = cells[:,filt_bl]

    # mark lateral genes if present
    if args.l is not None :
        print("mc2 | Marking lateral genes ...")
        LATERAL_GENE_NAMES = [i.replace("\n","") for i in open(args.l, "r")]
        mc.pl.mark_lateral_genes(cells,lateral_gene_names=LATERAL_GENE_NAMES)
        lateral_gene_mask = mc.ut.get_v_numpy(cells, "lateral_gene")
        lateral_gene_names = set(cells.var_names[lateral_gene_mask])

    # set parallelisation piles
    print("mc2 | Setting max piles ...")
    max_parallel_piles = mc.pl.guess_max_parallel_piles(cells)
    mc.pl.set_max_parallel_piles(max_parallel_piles)

    # set parameters not defined by user
    # target mc umis
    if args.N is None:
        # prior definitions
        median_umi_cell = int(np.median(cells.X.sum(axis=1).tolist()))
        n_cells = int(cells.X.shape[0])
        est_ncells_mc = 50 # ~50 cells per metacell, could be even less
        print("mc2 | estimated median umi cell = " + str(median_umi_cell))
        # calculation
        est_N = (median_umi_cell * est_ncells_mc) / np.log10(n_cells)
        est_N = int(round(est_N, -4)) # round to nearest 10000s
        # upper boundary
        if est_N > 300000:
            est_N = 300000
        # lower boundary
        if est_N < 30000:
            est_N = 30000
        print("mc2 | estimated target_metacell_umis = " + str(est_N))
        args.N = est_N
    TARGET_MC_UMIS = int(args.N)
    # rare_min_gene_maximum
    if args.r is None:
        args.r = 7
    MIN_GENE_MAX = int(args.r)
    # rare max gene cell fraction
    if args.f is None:
        args.f = 0.001
    RARE_MAX_CELL_FRAC = float(args.f)

    # run divide_and_conquer
    print("mc2 | mc2 algorithm, target_metacell_umis = " + str(args.N) + " ...")
    with mc.ut.progress_bar():
        mc.pl.divide_and_conquer_pipeline(
            cells,
            target_metacell_umis = TARGET_MC_UMIS,
            rare_min_gene_maximum = MIN_GENE_MAX,
            rare_max_gene_cell_fraction = RARE_MAX_CELL_FRAC,
            random_seed=123456
        )
    
    # give an actual proper name to the metacell, not only the number
    cells.obs["mc_name"] = [ "mc" + str(i) for i in cells.obs['metacell'].tolist() ]

    # write results, .h5ad
    print("mc2 | Writing output .h5ad ... ")
    outfn_h5ad = args.o + ".mc2.h5ad"
    cells.write_h5ad(outfn_h5ad)
    
    # write results, .tsv
    print("mc2 | Writing output .csv ...")
    outfn_csv = args.o + ".mc2.csv"
    cells.obs.to_csv(outfn_csv)
    
    # output number of metacells identified
    n_metacells = len(set(cells.obs['mc_name'].tolist()))
    print("mc2 | Found " + str(n_metacells) + " metacells.")

    return 0 # by convention, return 0 for success, non-zero otherwise


if __name__ == "__main__":
    """
    Main guard
    """
    # used to exit the program with the exit code returned by the main function
    sys.exit(main())
    # After runnung the script check in the command line the result: echo $?
