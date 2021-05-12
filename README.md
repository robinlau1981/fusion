Here is a copy of data Matlab codes involved in the following paper

@inproceedings{liu2021robust,

title={Robust Dynamic Multi-Modal Data Fusion: A Model Uncertainty Perspective},

author={Liu, Bin},

booktitle={to be filled},

pages={to be filled},

year={2021},

organization={to be filled} }

Appreciate if you cite the above paper.

The correspondance betwee algorithm names in the code and that in the paper is as follows:

"pf" in the code       <----------->  "PF" in the paper

"dmmpf" in the code     <----------->  "DMA" in the paper

"pf_df" in the code     <----------->  "SMA" in the paper

"pf_alpha" in the code  <----------->  "TS" in the paper

A brief description of the code files is as follows

main_alg_compare.m              : the main function for reproducing experimental results in the paper

simu_data.mat and simu_data2.mat: the two data sets used in the experiment 

simulate_data_gen.m             : code for generating simu_data.mat

simulate_data_gen2.m            : code for generating simu_data2.mat

residualR.m                     : residual resampling
