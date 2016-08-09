# DynamicPredJMmultiOutcomes
## Dynamic predictions using joint models of two longitudinal outcomes and competing risk data.

The code for modeling multiple longitudinal (2 continuous and 1 ordinal) and a competing risk model together with the dynamic predition assuming the value parameterization is presented. Specifically:
* "data.id1.Rdata" / "data1.Rdata": include simulated data for the longitudinal and survival outcomes.
* "Functions": includes functions needed for the code.
* "ModelJAGS": includes the joint model for jags.
* "PrepareData": includes the preparation of the data.
* "Fit": includes the main code. Specifically, it loads the data, runs the functions, creates the model, runs the model in JAGS and saves the results.
* (Not available yet) "DynPred": has the code to run the dynamic preditions assiming the value parameterization. 

How does it work:
* Download all files and place them in one folder.
* Set as working directory in R this folder.
* Run the code in "Fit" for fitting the joint model.
* (Not available yet) Run the code in "DynPred" for obtaining dynamic predictions.

