function data = sim_HJ_OV(sim_par)

% ------------------------------------------------------------------------
% SIMULATION
% ------------------------------------------------------------------------
%    y(t) = rho*y(t-1) + e(t) + v(t)
%    e(t) ~ N(0,s_e^2); v(t) ~ N(0,s_v^2)
% ------------------------------------------------------------------------

% LOAD PARAMETERS
s_e = 1; %sim_par.s_e;
s_v = 1; %sim_par.s_e;
rho_x = 0.75;
rho = sim_par.rho;
T   = sim_par.T;
s_eu = sim_par.s_eu;

% SIMULATION\
eu = mvnrnd(zeros(1,2),[1,s_eu;s_eu,1],T);
e = eu(:,1);
u = eu(:,2);
v = normrnd(0,s_v,T,1);
x = NaN(T,1);
y = NaN(T, 1);
x(1) = normrnd(0, sqrt(1/(1-rho_x^2))); % initialize at stationary distribution
y(1) = normrnd(0, sqrt((s_e^2 + s_v^2)/(1-rho^2))); % initialize at stationary distribution
for t = 2:T
    x(t) = rho_x*x(t-1) + u(t);
    y(t) = rho*y(t-1) + x(t) + e(t) + v(t);
end


% ------------------------------------------------------------------------
% SET UP FOR ESTIMATION
% ------------------------------------------------------------------------

% DATA
data.tab = array2table([y,e],'VariableNames',{'y','e'});

% FULL LIST OF VARIABLE ABBREVIATIONS
data.var_list_all = {'y','e'};

% % NUMBER OF VARIABLES
% data.n_var_total = 2;

% NUMBER OF PERIODS
data.T = T;


end