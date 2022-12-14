# proteomePhyto

## The model

Our coarse-grained model is built on a constrained optimization problem that aims to maximize the steady-state growth rate of a generic phytoplankton cell considering temperature, irradiance levels, and external dissolved inorganic nitrogen and carbon concentrations. 

The model is written in Julia and contains three files: 
1) par.jl defines the parameter values
2) eqn.jl defines the model equations
3) save.jl runs the model over environmental conditions and saves model output

Julia packages required to run the code above: JuMP, Ipopt, Plots 

## Model analyses 

Model analyses were carriedout in R. R code to vizualize and reproduce the figures in the paper are given in:
1) val.R for model validation
2) sensi.R for all other figures and sensitivity analyses



