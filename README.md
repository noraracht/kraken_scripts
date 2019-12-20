# Access to the data used for Kraken benchmarking

This repository contains summary data tables and scripts we used to processes them.


* Filtering on simulated Drosophila genome skims with non overlapping contaminant portions.
    - [Drosophila_contam_both_species_formatted.xls](https://github.com/noraracht/kraken_scripts/blob/master/Drosophila_contam_both_species_formatted.xls) contains information about genomic distances estimated before and after filtering for simulated Drosophila genomes where only one of the species was contaminated. [E4_script_NO_conf.R](https://github.com/noraracht/kraken_scripts/blob/master/E4_script_NO_conf.R) script takes this data table as an input and generates results for experiment represented in Fig. 4a.
     - [Drosophila_contam_both_species_formatted_withAlpha.xls](https://github.com/noraracht/kraken_scripts/blob/master/Drosophila_contam_both_species_formatted_withAlpha.xls) summarizes information about genomic distances estimated before and after filtering for the same experiment as above but reported for *k* = {28, 32, 35} and *α* = {0.00, 0.05}. [E3_script_two_conf.R](https://github.com/noraracht/kraken_scripts/blob/master/E3_script_two_conf.R) script takes this data table as an input and generates results plot in Fig. S8.

* Filtering on simulated Drosophila genome skims with overlapping contaminant portions.
    - [my-data-output_true_H.csv](https://github.com/noraracht/kraken_scripts/blob/master/my-data-output_true_H.csv) includes information about genomic distances estimated before and after filtering for simulated Drosophila skims with overlapping contaminants. Data also include actual percent overlap between contaminant portions computed throught *k*\-mer analysis. [dros-overlap-exp-Tol-db.R](https://github.com/noraracht/kraken_scripts/blob/master/dros-overlap-exp-Tol-db.R) script takes this data table as an input and generates results for experiment shown in Fig. 4b.


* Filtering real Drosophila genome skims.
    - [sum_report_skmer_14dros_real_cleaned_kraken_std_unmasked_0.0_k35.csv](sum_report_skmer_14dros_real_cleaned_kraken_std_unmasked_0.0_k35.csv) has information about percentages of reads classified by Kraken.
    - [drosophilaskims.csv](https://github.com/noraracht/kraken_scripts/blob/master/drosophilaskims.csv) contains relative distance error before and after filtering for every pair of Drosophila skims.
    - [E3_script_dros_real_noABS2.R](https://github.com/noraracht/kraken_scripts/blob/master/E3_script_dros_real_noABS2.R) script takes [sum_report_skmer_14dros_real_cleaned_kraken_std_unmasked_0.0_k35.csv](sum_report_skmer_14dros_real_cleaned_kraken_std_unmasked_0.0_k35.csv) and [drosophilaskims.csv](https://github.com/noraracht/kraken_scripts/blob/master/drosophilaskims.csv) tables above as an input and generates results for Figs. 5a, 5c, S9, S10.
    - [Drosophila_C_to_different_domains.csv](https://github.com/noraracht/kraken_scripts/blob/master/Drosophila_C_to_different_domains.csv) has percentages classified by Kraken for every species at different domain levels.
    - [Stacked_bar_plot_dros_dif_domains2.R](https://github.com/noraracht/kraken_scripts/blob/master/Stacked_bar_plot_dros_dif_domains2.R) script takes [Drosophila_C_to_different_domains.csv](https://github.com/noraracht/kraken_scripts/blob/master/Drosophila_C_to_different_domains.csv) as an input and generates plot in Fig. 5b.
    
    
* Sensitivity analysis of Kraken
    - [summary_report_kraken_viral_dist_unmasked_0.0_corrected_k35.csv](https://github.com/noraracht/kraken_scripts/blob/master/summary_report_kraken_viral_dist_unmasked_0.0_corrected_k35.csv) contains percent read classified at domain level or lower for every query sequence. It serves as an input for  [Cmnds_inScript_E1.R](https://github.com/noraracht/kraken_scripts/blob/master/Cmnds_inScript_E1.R) script to generate results depicted in Fig. S5 and S6. 
    
     - [export_df_all_data_extended_v2.csv](https://github.com/noraracht/kraken_scripts/blob/master/export_df_all_data_extended_v2.csv) contains FN, FP, TP, TN, FPR, recall computed per bin. It is an input for  [Cmnds_inScript_E2_extended_git_v2.R](https://github.com/noraracht/kraken_scripts/blob/master/Cmnds_inScript_E2_extended_git_v2.R) script to generate ROC curves depicted in Fig. 3. Fig. S7 was constructed based on [export_df_all_data_extended.csv](https://github.com/noraracht/kraken_scripts/blob/master/export_df_all_data_extended.csv) using [Cmnds_inScript_E2_extended.R](https://github.com/noraracht/kraken_scripts/blob/master/Cmnds_inScript_E2_extended.R) script.
   

   - [E1_exp_parser_viralDB.ipynb](https://github.com/noraracht/kraken_scripts/blob/master/E1_exp_parser_viralDB.ipynb) script to extract data about percent reads classified at domain level or below.
   
    - [kraken_report_parser_for binned output](https://github.com/noraracht/kraken_scripts/blob/master/kraken_report_parser_v4_binned_conf_matrix_-checkpoint.ipynb) reads output reports and computes statistics at the bin level. It requires [10Kmatrix](https://github.com/noraracht/kraken_scripts/blob/master/10K_dist_stats.csv), [excluded_genomes](https://github.com/noraracht/kraken_scripts/blob/master/Keep_exclude_corrected.csv) and [metadata](https://github.com/noraracht/kraken_scripts/blob/master/metadata.tsv).
   
* Theoretical exposition
    - [backenv.R](https://github.com/noraracht/kraken_scripts/blob/master/backenv.R)

* Phylogenetic analysis
    - [Drosophila-FME.R](https://github.com/noraracht/kraken_scripts/blob/master/Drosophila-FME.R) script generates Figure 5d and S11. It requires [FME_per_species.csv](https://github.com/noraracht/kraken_scripts/blob/master/FME_per_species.csv) and [amount.csv](https://github.com/noraracht/kraken_scripts/blob/master/amount.csv).

    - [jc.py](https://github.com/noraracht/kraken_scripts/blob/master/jc.py)
      performes Jukes-Cantor tranformation of distances output by Skmer.

    - [wrf.R](https://github.com/noraracht/kraken_scripts/blob/master/wrf.R)
      computes weighted and unweighted RF distance between two trees. 
      
* Mammalian genome analysis
    - Since mammalian genomes might have higher similarity with human reads we recommend two step filtering. [split_cseq_human_from_bac_v2.py](https://github.com/noraracht/kraken_scripts/blob/master/split_cseq_human_from_bac_v2.py) script performs normal filtering, then isolates reads classified as human and performs second filtering at higher confidence level. Output FASTQ file contains a combinations of reads that didn't match to the database after both rounds of filtering. Script take as an input one uncompressed *.fastq/.fq/* file per each sample. We provide [CerSimSim1.0_0.1x.fastq](https://drive.google.com/file/d/1MfXECtGZVj_dWSF_wiw26E84NLa_WwOE/view?usp=sharing) simulated  white rhinoceros reads as a test example. For a test run, put script into the same directory as test file and run
      ```
      python split_cseq_human_from_bac_v2.py -i CerSimSim1.0_0.1x.fastq -threads 24 -db $DBNAME -confidence1 0.0 -confidence2 0.5
      ```

      [Kraken-II](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual) should be installed and database constructed. User can specify a threshold for first (`-confidence1` option) and second (`-confidence2` option) confidence scores in the [0.0 1.0] interval. Default confidence values are 0.0 and 0.5 for first and second filtering steps, correspondingly.
      Replace "$DBNAME" above with your preferred database name/location. If multiple processing cores are available, number of threads can be modified with `-threads` option. Results can be assessed based on output reports.
 
  


