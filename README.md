# biocon_SAR

These are the code and data files associated with "Elevated CO2 and enriched nitrogen decrease species richness most at small spatial scales in a grassland experiment."

There are two code files: SAR_biocon_data.R processes the original data file hi_res_SR16_biocon.csv (with associated meta-data hi_res_SR16_biocon_meta.txt) into two data files (sr_df.csv and perm_boot.RData) which are used as input to SAR_biocon_models.R. The processed data are included here for convenience (the bootstrap generated estimates can be time consuming) and precise reproduction of th uncertainty estimates provided in the manuscript.
