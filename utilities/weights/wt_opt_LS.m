function  w = wt_opt_LS(Y, varargin)

% ========================================================================
% SOLVE FOR OPTIMAL WEIGHTS TO MAXIMIZE LOG SCORE LINEAR POOL
% ------------------------------------------------------------------------
%
% -- INPUT --
% either (pd_1, pd_2, ...) each a vector of predictive densities
%     or PD a cell with PD{i} = pd_i for i = 1, 2, ...
%
% -- OUPUT --
% w: weights for each model
%
% ------------------------------------------------------------------------

% COMBINE PREDICTIVE DENSITIES
if nargin > 1
    Yhat = cat(2,varargin{:}{:});   	% models entered seperately
else
    Yhat = cat(2,varargin{1}{:});	% models entered as single cell
end
Y = Y(all(~isnan(Yhat),2));
Yhat = Yhat(all(~isnan(Yhat),2),:);
n_model = size(Yhat,2);

% PREPARE FOR fmincon
fun = @(w) sum(( Y - Yhat*w' ).^2);      % objective function: log score
w_init = ones(1,n_model)/n_model;   % initial guess for optimization: equal weights on all models
lb = zeros(n_model,1);              % lower bound of zero
ub = ones(n_model,1);               % upper bound of one
Aeq = ones(1,n_model); beq = 1;     % weights sum to one
A = []; b = [];                     % no additional inequality constraints for weights

% SOLVE FOR OPTIMAL WEIGHTS
options = optimoptions('fmincon','Display','off');              % suppress printing
if min(log(Yhat'),[],'all') > -Inf
    w = fmincon(fun, w_init, A, b, Aeq, beq, lb, ub, [], options);  % solve
else
    w = NaN;
end
        
end