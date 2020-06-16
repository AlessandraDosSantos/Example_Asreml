# Example_Asreml
Comparing Banded structure using corb and own function

It was created two functions to model the residual correlation as banded structures. The idea was to compare these functions with the default corb function in asreml package. 
I did that to see with the results were similar with less singularity problems.

It was simulated two scenarios, with 1000 data sets each, for a design grid of 10 columns and 6 rows.
The models had genetic random effects and spatial residual effects model as AR(1) x corb(1) in column and row direction.  
