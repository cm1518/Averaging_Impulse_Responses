function data = sim_ARMA(sim_par)

% ------------------------------------------------------------------------
% SIMULATION
% ------------------------------------------------------------------------

% LOAD PARAMETERS
rho_y = sim_par.rho_y;
rho_u = sim_par.rho_u;
T_sim = sim_par.T_sim;
T_pre = sim_par.T_pre;

% TOTAL PERIODS TO SIMULATE
T_tot = T_pre + T_sim;

% GENERATE SHOCKS
e1 = randn(T_tot,1);
e2 = randn(T_tot,1)*sqrt(1-rho_u^2);

% SIMULATION
y  = zeros(T_tot,1);
u2 = zeros(T_tot,1);
for t = 2:T_tot
    u2(t) = rho_u*u2(t-1) + e2(t);
    y(t)  = rho_y*y(t-1)  + e1(t) + u2(t);
end

% OBSERVABLES FOR SUBSAMPLE
y =  y(end-T_sim+1:end);
v = e1(end-T_sim+1:end);


% ------------------------------------------------------------------------
% SET UP FOR ESTIMATION
% ------------------------------------------------------------------------

% DATA
data.tab = array2table([y,v],'VariableNames',{'y','v'});

% FULL LIST OF VARIABLE ABBREVIATIONS
data.var_list_all = {'y','v'};

% % NUMBER OF VARIABLES
% data.n_var_total = 2;

% NUMBER OF PERIODS
data.T = T_sim;

end