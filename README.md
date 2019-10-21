# Access to the data used for Kraken benchmarking

This repository contains summary data tables and scripts we used to processes them.

* [my-data-output_true_H.csv](https://github.com/noraracht/kraken_scripts/blob/master/my-data-output_true_H.csv) contains information about genomic distances estimated before and after filtering for simulated Drosophila genomes with overlapping contaminants. Data also include actual percent overlap between contaminant portions computed throught * *k* *\-mer analysis. [dros-overlap-exp-Tol-db.R](https://github.com/noraracht/kraken_scripts/blob/master/dros-overlap-exp-Tol-db.R) script takes this data table as an input [my-data-output_true_H.csv](https://github.com/noraracht/kraken_scripts/blob/master/my-data-output_true_H.csv) and generates results for contamination overlap experiment represented in Fig. 4b in a manuscript.


* Query summary reports and distance matrices used for simulation experiment with non overlapping contaminants:
    - [Drosophila_contam_non_overlap_exp.zip](https://github.com/noraracht/kraken_raw_data/blob/master/Drosophila_contam_non_overlap_exp.zip)



