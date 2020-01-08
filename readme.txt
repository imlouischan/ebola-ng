readme

- 0a_main_mcmc.R
The main script to create mcmc samples of all alternative models. 
- 0b_main_plot.R
The main script to plot and present mcmc results. 

- 1a_data.R
The script to read data, including time events and infectors. 
- 1b_bound.R
The script to create the functions for boundaries of missing time points. 
- 2a_incubation.R
The script to create the function for incubation period (and plot). 
- 2b_distribution.R
The script to create the function for time lag distributions using MLE (and plot). 
- 3a_pmf_r.R
The script to create the function for offspring distributions. 
- 3b_pmf_s.R
The script to create the function for serial interval distributions. 
- 4a_L_ri.R
The script to create the likelihood of offspring number of each case. 
- 4b_L_si.R
The script to create the likelihood of serial interval of each case. 
- 4c_logL.R
The script to create the likelihood of all cases. 
- 4d_pij.R
The script to create the transmission probability used in each mcmc step to select an infector. 
- 5_mcmc.R
The script to create mcmc samples using metropolis hastings. 
- 6a_parallel.R
The script to run several mcmc chains in parallel using the package called "parallel". 
- 6b_combine_chains.R
The script to combine mcmc chains. 

- 7a_plot_mixing.R
The script to plot mcmc chain samples. 
- 7b_plot_dependence.R (Figure 2)
The script to plot dependence between parameters. 
- 8a_plot_pmf_r.R (Figure 3B)
The script to plot offspring distribution. 
- 8b_plot_pmf_s.R (Figure 3A)
The script to plot serial interval distribution. 
- 9a_cal_pmf_r.R
The script to calculate mean and SD of offspring distribution. 
- 9b_cal_pmf_s.R
The script to calculate mean and SD of serial interval distribution. 
- 10a_comp_model.R (Table 1)
The script to calculate criterion values to compare alternative distributions. 
- 10b_comp_para.R (Figure 4)
The script to plot comparison of parameters in two scenarios. 
- 11_efficacy.R (Figure 3D)
The script to plot efficacy (individual protect effect) of all cases. 
- 12_infector.R (Figure 3C)
The script to plot distribution of infector of Case 20. 
- 13_timelines.R (Figure 1)
The script to plot observed and reconstructed time events of all cases. 

- ebola_ng.xlsx
The data including time events and infectors of 20 cases. 

- samples_alt_N1e+06.Rdata
The mcmc samples of all alternative models, of each one million iterations. 

- figure
The folder storing all figures generated using the codes. 
