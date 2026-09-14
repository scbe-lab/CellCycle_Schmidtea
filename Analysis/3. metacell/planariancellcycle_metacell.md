# Planarian cell cycle project ; metacell analysis

First we defined a list of "lateral genes" for Schmidtea mediterranea, which is basically a list of any gene labelled or linked to ribosomal biology, as per our annotation in the "Rosetta" table.
This file is the same that was used in: 
  Pérez-Posada, A., García-Castro, H., Emili, E. et al. 
  Multimodal single cell analyses reveal gene networks of planarian stem cell differentiation. 
  Nat Commun 16, 10683 (2025). https://doi.org/10.1038/s41467-025-65712-0
  supplementary data 2

```sh
cat Smed_Rosetta.tsv | grep "ibosom" | grep h1SM | cut -f1 | sort | uniq > smed_ribo.txt
```

Second, we ran the python script `mc2.py` which is a barebones implementation of the one-pass workflow of the `metacells` package (Ben-Kiki et al., 2022).

We used the following parameters:

 - `-i` to provide the input h5ad with raw counts of genes x cells
 - `-l` to provide the list of lateral genes
 - `-N` to provide the target_umi_metacell which we fixed at ~3,000 UMIs per metacell (we saw similar results with higher UMIs as well),
 - `-o` to provide the basename for the output files

```sh
mc2.py -i Smed_L78-L47_20250523_unprocessed.h5ad -l smed_ribo.txt -N 3000 -o smed_cdh1_3k
```

The resulting table `smed_cdh1_3k.mc2.csv` was imported to the adata object using scanpy and pandas.

The code of the python script is provided in the `mc2.py` script alongside this markdown.
