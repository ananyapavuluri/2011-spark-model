# 2011-spark-model
Calcium sparks are events at the cellular level characterized by the release of calcium ions from the sarcoplasmic reticulum via ryanodine receptors.
In cardiac myocytes, upon membrane depolarization, voltage-gated L-type calcium channels open and allow for the entrance of calcium ions into the cell, raising the concentration of intracellular calcium. 
Calcium ions bind to ryanodine receptors on the membrane of the junctional sarcoplasmic reticulum and cause the release of calcium ions from the sarcoplasmic reticulum (a positive feedback mechanism called calcium-induced calcium release).

Calcium sparks can also be spontaneous (occuring without any electrical stimulus). Spontaneous sparks are caused by the stochastic opening of ryanodine receptors.
Upon discoveries about the geometry and behavior of ryanodine receptors, Sobie et al. created the sticky cluster model (see publication [here](https://www.ncbi.nlm.nih.gov/pubmed/12080100)).
In 2011, parameters in the sticky cluster model were adjusted to reflect new experimental data by [Ramay et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3156908/).

The 2011 sticky cluster model was written in MATLAB. This project aims to improve the speed and efficiency of the spark simulation using [CUDA], (https://www.nvidia.com/object/io_69526.html) a parallel computing platform that uses a graphics processing unit to perform general purpose computations. Running several simulations in parallel drastically decreases the time spent waiting to acquire data from the mathematical model. 

The CUDA and C++ programs in this repository output data in the form of CSV files, which are then read into MATLAB for plotting and analysis. The MATLAB script used to generate figures for this model is also provided. 

This project was done in the Cardiac Systems Pharmacology Lab of Dr. Eric Sobie, Department of Pharmacological Sciences at the Icahn School of Medicine of Mount Sinai (New York, NY), and with the guidance of M.D./Ph.D. candidate Deanalisa Jones.

To run the C++ version of the model on the command line:

``` 
g++ CaSparkModel2011.cpp
./a.out
```

The CUDA version of the model was created using Microsoft Visual Studio 2015.
