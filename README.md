# Code developed to investigate new approaches to cerebrovascular reactivity mapping using MRI

The code in this repository was used to analyse two approaches to cerebrovascular reactivity mapping. The data acquired in this study can be found at https://doi.org/10.5287/bodleian:Xk48adQAO. Please reference this code using the forthcoming Zenodo DOI, or for the dataset the DOI above, if you use either in your work.

Data were analysed using MATLAB version 8.4.0.150421 (R2014b) and FSL (FMRIB Software Library) version 5.0.7 on Mac OS X version 10.11.6. Additional MATLAB toolboxes are included to read JSON files (https://github.com/mathbiol/couch4mat) and to plot points with X and Y errorbars (https://uk.mathworks.com/matlabcentral/fileexchange/40221-plot-data-with-error-bars-on-both-x-and-y-axes).

This code assumes the data are organised along the lines of the BIDS format (http://bids.neuroimaging.io). Calculated images should be placed in a folder within the data folder with the name "derivatives". Within this folder the code in this repository should be placed in a folder named "code".
* cvr_study - data folder
* cvr_study/derivatives - processed data folder
* cvr_study/derivatives/code - code for processing data folder