function out = OLSfreq_est(Y, X, varargin)

% ========================================================================
% FREQUENTIST OLS ESTIMATION
% ========================================================================
%
% MODEL:
%   Y = X*b + e
%
% INPUT:
%   Y (T-by-1): dependent variable
%   X (T-by-N): independent variables
%   se: type of standard errors ('white', 'nw'); homoskedastic if omitted
%
% OUTPUT:
%   sqrt(T)*(beta-b) ~ N(0, V*T)
%                  e ~ N(0, s2)
%
% ========================================================================


% ------------------------------------------------------------------------
% PARSE INPUTS
% ------------------------------------------------------------------------

p = inputParser;
addOptional(p, 'constant', 'true', ... 
    @(x) any(validatestring(x,{'true', 'false'})));
addOptional(p, 'sample', 1:size(Y,1), ...
    @(x) isnumeric(x) && min(floor(x)==x));
addOptional(p, 'se', 'homoskedastic',...
    @(x) any(validatestring(x,{'homoskedastic','white','nw'})));
parse(p, varargin{:})
out = p.Results;


% ------------------------------------------------------------------------
% CLEAN DATA
% ------------------------------------------------------------------------

% select sample
Y = Y(out.sample,:);
X = X(out.sample,:);

% remove rows with NaN
no_NaN = all(~isnan([Y,X]),2);
Y = Y(no_NaN,:);
X = X(no_NaN,:);
out.sample_no_NaN = out.sample(no_NaN); %AG added

% add constant if needed
if strcmp(out.constant, 'true')
    X = [X, ones(size(X,1),1)];
end

% dimensions
[T, N] = size(X);


% ------------------------------------------------------------------------
% POINT ESTIMATES
% ------------------------------------------------------------------------

b = X\Y;
e = Y - X*b;


% ------------------------------------------------------------------------
% STANDARD ERRORS
% ------------------------------------------------------------------------

XX_i = eye(N)/(X'*X);

if strcmp(out.se, 'homoskedastic')
    % homoskedastic std errors
    V = sum(e.^2)/(T-N)*XX_i;

elseif strcmp(out.se, 'white')
    % White heteroskedasticity-consistent std errors
    Q = X'*diag(e.^2)*X;
    V = XX_i*Q*XX_i * T/(T-N);

elseif strcmp(out.se, 'nw')
    % Newey-West std errors
    M = floor(0.75*T^(1/3)); % Andrews (1991)
    % M = floor(1.3*sqrt(T)); % Lazarus, Lewis, Stock, Watson (2018)
    Q1 = zeros(N,N,T,M);
    for j = 1:M
        for i = j+1:T
            Q1(:,:,i,j) = (1-j/(M+1))*e(i)*e(i-j)...
                *(X(i,:)'*X(i-j,:) + X(i-j,:)'*X(i,:));
        end
    end
    Q = X'*diag(e.^2)*X + sum(Q1,[3,4]);
    V = XX_i*Q*XX_i * T/(T-N);
end


% ------------------------------------------------------------------------
% SAVE OUTPUT
% ------------------------------------------------------------------------

out.B = b;
out.V = (V+V')/2;	% make sure V is symmetric (in case of numerical error)


end