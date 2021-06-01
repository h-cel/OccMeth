# OccMeth
Comparison of occurrence methods

Cite as

## General overview

This is the code used for the calculations and evaluation in Van de Velde et al. (2021). 
In this article, 2 precipitation occurrence bias adjustment methods (also known as bias correction methods) are compared with the application of Quantile Delta Mapping (QDM) without any precipitation occurrence bias adjustment step.

## Structure

In this code, you can find the following files. An explanation and the dependency is included in each file.

* a_loadClimateData
* b_configurationBiasCorrection
  * prepareBiasdata
  * BiasCorrection
      * occAdj_TDA
        * T
      * occAdj_SSRmonthly
      * QDM
      * postprocessingSSR
* c_BCEvaluation
  * TruncateObs
  * BC_Evaluation
    * fullfig
* Visualisation

## Detailed overview

Four main (or configuration) files are included in this code:
* a_loadClimateData: loads the climate data
* b_configurationBiasCorrection: launches the bias correction/bias adjustment
* c_BCEvaluation: calculates all indices and the RB_O and RB_MB values as used in the article
* Visualisation: uses the RB_O and RB_MB values calculated to make the graphs used in the article

Last update on 01/06/'21 by J. Van de Velde
