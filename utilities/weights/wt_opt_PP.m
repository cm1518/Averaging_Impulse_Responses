function  w = wt_opt_PP(varargin)

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
    pd = cat(2,varargin{:});   	% models entered seperately
else
    pd = cat(2,varargin{1}{:});	% models entered as single cell
end
pd = pd(all(~isnan(pd),2),:);
n_model = size(pd,2);

% PREPARE FOR fmincon
fun = @(w) -sum(log( w*pd' ));      % objective function: log score
w_init = ones(1,n_model)/n_model;   % initial guess for optimization: equal weights on all models
lb = zeros(n_model,1);              % lower bound of zero
ub = ones(n_model,1);               % upper bound of one
Aeq = ones(1,n_model); beq = 1;     % weights sum to one
A = []; b = [];                     % no additional inequality constraints for weights

% SOLVE FOR OPTIMAL WEIGHTS
options = optimoptions('fmincon','Display','off');              % suppress printing
w = fmincon(fun, w_init, A, b, Aeq, beq, lb, ub, [], options);  % solve
        
end