# 2011-spark-model
Calcium sparks are events at the cellular level characterized by the release of calcium ions from the sarcoplasmic reticulum via ryanodine receptors.
In cardiac myocytes, upon membrane depolarization, voltage-gated L-type calcium channels open and allow for the entrance of calcium ions into the cell, raising the concentration of intracellular calcium. 
Calcium ions bind to ryanodine receptors on the membrane of the junctional sarcoplasmic reticulum and cause the release of calcium ions from the sarcoplasmic reticulum (a positive feedback mechanism called calcium-induced calcium release).

Calcium sparks can also be spontaneous (occuring without any electrical stimulus). Spontaneous sparks are caused by the stochastic opening of ryanodine receptors.
Upon discoveries about the geometry and behavior of ryanodine receptors, Sobie et al. created the sticky cluster model (see publication [here](https://www.ncbi.nlm.nih.gov/pubmed/12080100)).
In 2011, parameters in the sticky cluster model were adjusted to reflect new experimental data by [Ramay et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3156908/).

The 2011 sticky cluster model was written in MATLAB. This project aims to improve the speed and efficiency of the spark simulation using CUDA. Plotting of data is still done in MATLAB using CSV files from C++ (.csv files in this repo are C++ output, .eps files are MATLAB plots, sticky_plot.m is MATLAB code used to plot). Provided plots and data from the C++ program reflect a simulation with 3 trials and 3,000,000 iterations. 

This project was done in the Cardiac Systems Pharmacology Lab of Dr. Eric Sobie, Department of Pharmacological Sciences at the Icahn School of Medicine, Mount Sinai.

To run the C++ version of the model:

``` 
g++ CaSparkModel2011.cpp
./a.out
```
