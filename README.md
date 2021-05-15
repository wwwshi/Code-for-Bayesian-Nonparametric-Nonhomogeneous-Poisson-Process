# Code-for-Bayesian-Nonparametric-Nonhomogeneous-Poisson-Process

This repository contains code for [Bayesian Nonparametric Nonhomogeneous Poisson Process with Applications to USGS Earthquake Data](https://www.sciencedirect.com/science/article/pii/S2211675321000051?casa_token=G3JaDdwC1IQAAAAA:xEPCJck2j5vPH5SQmniH40nH6gORp6LnAvnEIHZixf4YxUwZEIyuVOjxuuwq8-A2Z0jOghsgdm0) by Junxian Geng, Wei Shi and Guanyu Hu, published in *Spatial Statistics*.

`data_train.Rdata` contains training data set for earthquake real data.

`data_test.Rdata` contains test data set for earthquake real data.

`PP_Rcpp_update_z.cpp` contains Rcpp functions for calculating VN and updating cluster assignment z.

`function.R` contains helper functions for initialization and updating.

`main_realdata.R` contains functions for our Collapsed Gibbs sampler algorithm and Dahl's method for summarizing MCMC samples, the code to fit training data set and MAE comparison based on test data set.
