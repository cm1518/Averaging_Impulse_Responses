REPLICATION FILES FOR "AVERAGING IMPULSE RESPONSES USING PREDICTION POOLS"
Paul Ho, Thomas Lubik, and Christian Matthes
10/17/2023

--------------------------------------------------------------------------

<< MAIN FILES FOR COMPUTATIONS AND PLOTS >>

main_XYZ.m files: estimation, predictive densities, and optimal weights
>> start with simulation/application specific preamble
>> followed by main_calc_simulation or main_calc_application to carry out computations
>> see main_app_TFP for prototypical empirical application setup

plot_XYZ.m files: plot all figures in paper

suffixes:
>> sim_HJ       : Monte Carlo with AR(1) model (Section 3.1 and Appendix B)
>> sim_peristent: Monte Carlo with persistent shock (Section 3.2)
>> sim_SW       : Monte Carlo with Smets Wouters (2007) model (Appendix C)
>> app_Ramey    : monetary shock empirical application (Section 4.1 and Appendix D)
>> app_TFP      : TFP shock empirical application (Section 4.2)


<< utilities FOLDER >>
- calc_simulation/calc_application: carry out impulse response averaging calculations
- FAVAR/LP/sgl eqn/VAR subfolder  : functions for estimation, IRFs, and predictive densities for each model
- weights subfolder               : functions for averaging


<< out FOLDER >>
output from Monte Carlos and empirical applications for plotting