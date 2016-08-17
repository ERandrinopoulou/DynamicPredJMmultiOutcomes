# DynamicPredJMmultiOutcomes
## Dynamic predictions using joint models of two longitudinal outcomes and competing risk data.

This repository includes the code for (1) modeling multiple longitudinal (1 continuous and 1 ordinal) and a competing risk model and (2) dynamic predictions assuming the value parameterization. 

Specifically:
* "**data.id1.Rdata**" / "**data1.Rdata**": include simulated data for the longitudinal and survival outcomes.
* "**Functions**": includes functions needed for the code.
* "**ModelJAGS**": includes the joint model for jags.
* "**PrepareData**": includes the preparation of the data.
* "**Fit**": includes the main code. Specifically, it loads the packages and the data, runs the functions, creates the model, runs the model in JAGS and saves the results.
* (Not available yet) "**DynPred**": has the code to run the dynamic preditions assuming the value parameterization. 

How does it work:
* Download all files and place them in one folder.
* Set as working directory in R this folder.
* Run the code in "Fit" for fitting the joint model.
* (Not available yet) Run the code in "DynPred" for obtaining dynamic predictions.

