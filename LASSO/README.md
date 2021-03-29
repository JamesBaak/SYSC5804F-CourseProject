# LASSO
This file contains the needed code for the LASSO channel reconstruction.

## Running the LASSO_channel_estimation script
Ensure that the Dependencies folder has been added to Matlab's path. This can be done by right-clicking on the Dependencies folder from within Matlab and then "Add to Path" > "Selected Folders and Subfolders". This will import the dependencies for Mark Schmidt's code.

## LASSO Variations
Due to the fact that the values being used for MIMO reconstruction are all complex, the normal MATLAB Lasso function cannot be used. In order to use LASSO regression a new function needs to be imported.

1. Matlab Code by Mark Schmidt:
    Within Mark Schmidt's Matlab code available online there is a provided method for LASSO of complex values.
    Source https://www.cs.ubc.ca/~schmidtm/Software/code.html


2. Matlab LASSO: Another method discussed at https://stats.stackexchange.com/questions/469653/implementing-complex-lasso-in-matlab is a method similar to Mark's (1) where the real and complex components are transformed into a real matrix which includes both the real and complex coefficients. 
