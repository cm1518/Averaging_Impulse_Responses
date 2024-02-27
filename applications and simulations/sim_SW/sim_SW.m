function data = sim_SW(sim_par)

% ------------------------------------------------------------------------
% SOLVE MODEL
% X(t) = T1*X(t-1) + T0*e(t), e(t)~N(0,Se)
% ------------------------------------------------------------------------
[~, T1, T0, Se, pp, ii] = SW_sol;
[n_state, n_shock] = size(T0);


% ------------------------------------------------------------------------
% COMPUTE IMPLIED OBSERVABLES
% Y(t) = A + H*X(t)
% ------------------------------------------------------------------------
n_series = 7;
H = zeros(n_series,n_state);
A = zeros(n_series,1);

% 1. GDP GROWTH
A(1) = pp.gam_bar;  H(1,ii.y) = 1;  H(1,ii.y_1) = -1;

% 2. INFLATION
A(2) = pp.pii_bar;  H(2,ii.pii) = 1;

% 3. FED FUNDS RATE
A(3) = pp.r_bar;    H(3,ii.r) = 1;

% 4. CONSUMPTION GROWTH
A(4) = pp.gam_bar;  H(4,ii.c) = 1;  H(4,ii.c_1) = -1;

% 5. INVESTMENT GROWTH
A(5) = pp.gam_bar;  H(5,ii.i) = 1;  H(5,ii.i_1) = -1;

% 6. WAGE GROWTH
A(6) = pp.gam_bar;  H(6,ii.w) = 1;  H(6,ii.w_1) = -1;

% 7. HOURS
A(7) = pp.l_bar;    H(7,ii.l) = 1;


% ------------------------------------------------------------------------
% SIMULATE MODEL FORWARD
% ------------------------------------------------------------------------

% LOAD PARAMETERS
T_sim = sim_par.T_sim;
T_pre = sim_par.T_pre;

% TOTAL PERIODS TO SIMULATE
T_tot = T_pre + T_sim;

% SIMULATION
e = mvnrnd(zeros(n_shock,1),Se,T_tot)'; % innovations
X = zeros(n_state, T_tot);              % model variables
for t = 2:T_tot
    X(:,t) = T1*X(:,t-1) + T0*e(:,t);
end
X = X(:,T_pre+1:end);   % extract observed periods


% ------------------------------------------------------------------------
% OBSERVED DATA
% ------------------------------------------------------------------------
Y = A + H*X;    % data
ind = 5;        % index for monetary shock
nse = 0.0;      % amount of noise in instrument relative to std dev of shock (=0 for perfectly observed)
m = e(ind,T_pre+1:end) + normrnd(0,nse*pp.sig_r,1,T_sim);   % instrument


% ------------------------------------------------------------------------
% SET UP FOR ESTIMATION
% ------------------------------------------------------------------------

% FULL LIST OF VARIABLE ABBREVIATIONS
data.var_list_all = {'GDP', 'Inflation', 'FFR', 'Consumption', 'Investment', 'Wage', 'Hours', 'Instrument'};

% DATA
data.tab = array2table([Y',m'],'VariableNames',data.var_list_all);

% NUMBER OF PERIODS
data.T = T_sim;

% 
% % ------------------------------------------------------------------------
% % COMPUTE TRUE IRF
% % ------------------------------------------------------------------------
% HH = 20; % number of horizons
% 
% % IRF for model variables
% irf_X = NaN(n_state, HH+1); 
% irf_X(:,1) = T0(:,ind)*sqrt(Se(ind,ind)); % XXX PAUL UPDATE 6/6/22: added std dev
% for h = 1:HH
%     irf_X(:,h+1) = T1*irf_X(:,h);
% end
% 
% % IRF for data
% irf_Y = H*irf_X;


end