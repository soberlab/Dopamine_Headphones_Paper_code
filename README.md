# Dopamine_Headphones_Paper_code
This repository has the code used to generate figures and perform data analysis for the Dopamine Headphones paper by Saravanan et al, 2019. Note that the raw data is not provided and so these codes cannot be used to replicate the figures as is.
Files present in this repository:

1. make_final_ohda_figures_corrected_bootstrap.m
This file contains the code to produce all the figures for this paper. It depends on the following two files for analysis which have also been made available. To reiterate though, the datasets have not been made available. This file has been commented for clarity.

1.1 load_bird_params_2018.m
This file is a switch case file that allows me to collect all the relevant data from individual birds for analysis by referencing them with a single number. Note that this file includes birds whose data was used for separate experiments not reported in this paper.

1.2. get_full_bird_phones_data_vs.m
In conjunction with load_bird_params_2018.m, this file is used to load all the data from birds that were subjected to headphones experiments.

2.  stats_for_ohda_figures_direct_probs.m
This file contains the code to perform all statistical analysis for the paper. This file has been commented and seeded for reproducibility of stats. It depends on the above two files in addition to the one mentioned below.

2.1. get_direct_prob.m
A small function written to compute the direct posterior probability of one sample of bootstrapped values being greater than or equal to a second sample.

3. make_figures_to_quantify_6ohda_lesion_in_area_x_2018.m
A file I used to quantify the extent of lesions in 6-OHDA injected birds at a population level. Note that this code has not been commented well and may be difficult to follow. The raw data for this code has also not been provided.

4. get_data_from_sober_brainard_2009_supp_fig_6.m
This file was used to transform data points extracted from an eps version of the Supplementary figure 6 from Sober and Brainard 2009 that was reported in Figure 3a of this paper. 

Data can be provided to interested parties upon reasonable request.
