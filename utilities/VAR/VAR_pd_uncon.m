function pd_out = VAR_pd_uncon(data, draw, test, ind, H, identification)

% ========================================================================
% COMPUTE PREDICTIVE DENSITY FROM VAR
% ========================================================================
%
% -- MODEL --
%   y(t) = con + B*[y(t-1),...,y(t-n_lag)] + C*u(t)
%
% -- INPUT --
%   draw: structure with coefficient draws
%   test: observations to compute predictive density for
%   ind: variable to compute predictive density for
%   H: horizon over which to compute predictive densities
%   identification: identification scheme
%       >> 'shock-augmented': Cholesky, instrument ordered first (Plagborg-Moller and Wolf, 2021)
%       >> 'proxy': project instrument on residuals (Mertens and Ravn, 2013)
%   
% -- OUTPUT --
%   pd(T_test-by-H): predictive densities
%
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% PRELIMINARIES
% ------------------------------------------------------------------------

% LOAD DATA
if strcmp(identification,'internal IV')
    Y = [data.m, data.Y]; 
elseif strcmp(identification,'Cholesky') || strcmp(identification,'proxy') 
    Y = data.Y;
end

% DIMENSIONS
[T, N] = size(Y);           % number of periods, variables
n_draw = size(draw.B,3);	% number of draws
n_lag = size(draw.B,2)/N;   % number of lags
T_test = length(test);      % number of observations in test set

% LAGGED Y
X = lagmatrix(Y, 1:n_lag);

% MATRICES FOR STORAGE
pd      = NaN(H+1, T_test, length(ind), n_draw);
EY_draw = NaN(H+1, T_test, length(ind), n_draw);
VY_draw = NaN(H+1, T_test, length(ind), n_draw);



% ------------------------------------------------------------------------
% PREDICTIVE DENSITIES
% ------------------------------------------------------------------------

for ii = 1:n_draw
    
    % EXTRACT PARAMETER DRAWS
    cons = draw.cons(:,ii);
    B = draw.B(:,:,ii);
    V = draw.V(:,:,ii);
    
    % CONSTRUCT NECESSARY OBJECTS
    B_comp = [B; eye(N*(n_lag-1)),  zeros(N*(n_lag-1),N)];	% companion form coefficient
    
    % PRE-COMPUTE VARIANCES
    VY = zeros(N*n_lag);        % matrix for conditional variance
    for hh = 1:H+1
        VY = B_comp*VY*B_comp';
        VY(1:N,1:N) = VY(1:N,1:N) + V;
        for jj = 1:length(ind)
            VY_draw(hh,:,jj,ii) = VY(ind(jj),ind(jj));      % variance of interest
        end
    end
    
    % COMPUTE EXPECTATIONS
    for tt = 1:T_test        
        for hh = 1:H+1
            if test(tt)+hh-1 <= T
            if hh == 1 
                EY_prev = X(test(tt),:)';    % previous observations
                EY = NaN(N*n_lag,1);        % vector for conditional expectation
                EY(1:N) = cons + B*EY_prev;
            else
                % HORIZON >0: ITERATE FORWARD
                % mean
                EY_prev = EY;
                EY(1:N) = cons + B*EY_prev;
            end
            EY(N+1:end) = EY_prev(1:N*(n_lag-1));	% lagged entries from previous expectation
            for jj = 1:length(ind)
                EY_draw(hh,tt,jj,ii) = EY(ind(jj));     % expectation of interest
            end
            end
        end
    end

end

% PREDICTIVE DENSITY
for hh = 1:H+1
    tmp = test+hh-1 <= T;
    test_sub = test(tmp);
    for jj = 1:length(ind)
        pd(hh,tmp,jj,:) = normpdf(Y(test_sub+hh-1,ind(jj)),...
            squeeze(EY_draw(hh,tmp,jj,:)), squeeze(sqrt(VY_draw(hh,tmp,jj,:))));
    end
end
pd = mean(pd,4);
% for hh = 1:H+1
%     tmp = test+hh-1 <= T;
%     test_sub = test(tmp);
%     pd(hh,tmp,:) = normpdf(Y(test_sub+hh-1,ind),...
%         squeeze(EY_draw(hh,tmp,:)), squeeze(sqrt(VY_draw(hh,tmp,:))));
% end
% pd = mean(pd,3)';

% EY = mean(EY_draw,3); %AG added. Think these lines can be deleted (AG
% 4/11/2022)
% VY = mean(VY_draw,3); %AG added so that the checking function can work
pd_out.pd = permute(pd,[2,3,1]);
% pd_out.pd = pd; pd_out.EY_draw=EY_draw; pd_out.VY_draw=VY_draw; pd_out.EY=EY; pd_out.VY=VY;

end