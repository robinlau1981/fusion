Here is a copy of Matlab codes involved in the following paper

@inproceedings{liu2021robust,

title={Robust Dynamic Multi-Modal Data Fusion: A Model Uncertainty Perspective},

author={Liu, Bin},

booktitle={to be filled},

pages={to be filled},

year={2021},

organization={to be filled} }

Appreciate if you cite the above paper.

The corresponce betwee algorithm names in the code and that in the paper is as follows:

Code                     Paper

pf        <----------->  PF

dmmpf     <----------->  DMA

pf_df     <----------->  SMA

pf_alpha  <----------->  TS

main_alg_compare.m              : the main function for reproducing experimental results in the paper

simu_data.mat and simu_data2.mat: the two data sets used in the experiment 

simulate_data_gen.m             : code for generating simu_data.mat

simulate_data_gen2.m            : code for generating simu_data2.mat

residualR.m                     : residual resampling
