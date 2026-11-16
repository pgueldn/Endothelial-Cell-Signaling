Endothelial signaling Monte Carlo models
========================================

MATLAB code and numerical data for the endothelial-cell manuscript.

Quick start
===========

Open this folder as the MATLAB Current Folder, then run:
run_EC_export_publication_figures

This recreates Figures 2-5 from saved numerical data, without running
Monte Carlo simulations. Edit the options at the beginning of the script to change formats or request
numerical labels instead of check/cross marks. 

Run test_package for data-consistency checks and one baseline ODE integration.

Tested with MATLAB R2025a. 

Run new searches
================

Run these scripts:

run_EC_MonteCarlo_constant_parameter_validation_fit
run_EC_MonteCarlo_parameter_sweep_validation_fit
run_EC_MonteCarlo_constantspineperturbedmodel

The first two scripts each default to 100,000 samples and save their top 50
parameter sets. Change numMC near the top for a short trial. They write to
results/constant/ and results/reaction_specific/, respectively. Full baseline
trajectory storage is disabled to reduce memory usage; this does not change
sampling, ODE solutions, or validation scoring.

The hybrid script reads those new results, reranks the source models using the
current validation definition, and searches binary donor-substitution masks.
It retains tied-best constant spines and selects near-top perturbed donors.
Each hybrid uses one constant spine and one coherent donor, not independently
chosen donors for each edge. Defaults are 250 individuals, 200 generations,
five finalist pairs and five seeds.  

Data and parameter definitions
==============================

Location: Contents
data/constant/: All 100,000 scores and global triplets, best parameters, top scores, validation responses and sampling bounds
data/reaction_specific/: All 100,000 scores, best parameters, top-50 parameter sets, edge mapping, sampling bounds and validation responses
data/hybrid/: Archived ensemble rankings, parameter ranges, configuration and representative model
data/figures/: Figure 2 numerical reconstruction, Figure 3 statistics and all 121 Figure 5 WSS samples
data/model/: Source network workbooks for AngII_v5 and the historical v6 reconstruction
data/validation/: Literature bibliography

CreateValidationMatrix.m is the executable validation definition. Entries are
1 (increase), -1 (decrease), 0 (no change), and +/-0.5 (unknown/mixed, not scored).
The current definition contains 79 scored comparisons.

ODE_master.m integrates the supplied ODE expressions with ode45 over 0-35
model-time units, returning 100 requested time samples. Inputs start at 0.2,
with Stress set to 0.5 and exogenous AngII inactive. Each validation perturbation
adds 0.5 to one input in the final baseline state. Scoring uses the final
perturbed activity minus the final baseline activity, not a fold change.
New searches classify abs(difference) < 0.02 as unchanged. A difference exactly
at the threshold retains its sign. The solver endpoint is not an independent
proof of convergence to steady state.

