# Access to the data used for Kraken benchmarking

This repository contains summary data tables and scripts we used to processes them.


* Filtering on simulated Drosophila genome skims with non overlapping contaminant portions.
    - [Drosophila_contam_both_species_formatted.xls](https://github.com/noraracht/kraken_scripts/blob/master/Drosophila_contam_both_species_formatted.xls) contains information about genomic distances estimated before and after filtering for simulated Drosophila genomes where only single genomes was contaminated. [E4_script_NO_conf.R](https://github.com/noraracht/kraken_scripts/blob/master/E4_script_NO_conf.R) script takes this data table as an input and generates results for experiment represented in Fig. 4a.
     - [Drosophila_contam_both_species_formatted.xls](https://github.com/noraracht/kraken_scripts/blob/master/Drosophila_contam_both_species_formatted.xls) contains information about genomic distances estimated before and after filtering for simulated Drosophila genomes where only single genomes was contaminated. [multiple *k* and \al    ](https://github.com/noraracht/kraken_scripts/blob/master/E3_script_two_conf.R) script takes this data table as an input and generates results for experiment represented in Fig. 4a.





* Filtering on simulated Drosophila genome skims with overlapping contaminant portions.
    - [my-data-output_true_H.csv](https://github.com/noraracht/kraken_scripts/blob/master/my-data-output_true_H.csv) contains information about genomic distances estimated before and after filtering for simulated Drosophila genomes with overlapping contaminants. Data also include actual percent overlap between contaminant portions computed throught *k*\-mer analysis. [dros-overlap-exp-Tol-db.R](https://github.com/noraracht/kraken_scripts/blob/master/dros-overlap-exp-Tol-db.R) script takes this data table as an input and generates results for experiment shown in Fig. 4b.


* [my-data-output_true_H.csv](https://github.com/noraracht/kraken_scripts/blob/master/my-data-output_true_H.csv) contains information about genomic distances estimated before and after filtering for simulated Drosophila genomes with overlapping contaminants. Data also include actual percent overlap between contaminant portions computed throught *k*\-mer analysis. [dros-overlap-exp-Tol-db.R](https://github.com/noraracht/kraken_scripts/blob/master/dros-overlap-exp-Tol-db.R) script takes this data table as an input and generates results for contamination overlap experiment represented in Fig. 4b in the manuscript.



