# DynamicPredJMmultiOutcomes
Dynamic predictions using joint models of two longitudinal outcomes and competing risk data.

The code for modeling multiple longgitudinal (2 continuous and 1 ordinal) and a competing risk model together with the dynamic predition assuming the value parameterization is presented.

"data.id1.Rdata" / "data1.Rdata": include simulated data for the longitudinal and survival outcomes.
"Functions": includes all the functions needed for the code.
"ModelJAGS": includes the joint model for jags.
"PrepareData": includes the code in order to prepare the data.
"Fit": includes the main code. Specifically, it loads the data, runs the fucntions, save the model, run the model in JAGS and saves the results.
"DynPred": includes the code to run the preditions. In order to run this code, results from "Fit" should be first obtained.
