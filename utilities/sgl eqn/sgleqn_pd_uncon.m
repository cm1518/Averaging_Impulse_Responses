function pd_out = sgleqn_pd_uncon(est, draw, test, H)

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
%   ind: variables to compute predictive density for
%   H: horizon over which to compute predictive densities
%   identification: identification scheme
%       >> 'Cholesky': Cholesky, last shock
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
% Y = data.Y;
% if strcmp(identification,'shock-augmented')
%     Y = [data.m, data.Y]; % ind = ind + 1; IMPORTANT: We are assuming the ind fed in is already accounting for being shock augmented 
% % *** NEW - ALLOW FOR CHOLESKY **** %
% elseif strcmp(identification,'Cholesky') || strcmp(identification,'proxy') 
%     Y = data.Y;
% end
Y = est.data.Y;
z = est.data.m;

% DIMENSIONS
T = length(Y);              % number of periods, variables
n_draw  = size(draw.B,2);	% number of draws
n_lag_y = est.n_lag;
n_lag_z = size(draw.B,1) - n_lag_y - 1;
T_test = length(test);      % number of observations in test set

% LAGGED Y AND X
Y_lag = lagmatrix(Y, 1:n_lag_y+1);
Z_lag = lagmatrix(z, 0:n_lag_z-1); % assume we include contemporaneous z
% Y = lagmatrix();
% X = lagmatrix(Y, 1:n_lag);

% MATRICES FOR STORAGE
pd      = NaN(H+1, T_test, 1, n_draw);
EY_draw = NaN(H+1, T_test, 1, n_draw);
VY_draw = NaN(H+1, T_test, 1, n_draw);


% ------------------------------------------------------------------------
% PREDICTIVE DENSITIES
% expectations and variances are H x T_test x n_ind x n_draw
% ------------------------------------------------------------------------

vz = var(z);
mz = mean(z);
for ii = 1:n_draw
    
    % EXTRACT PARAMETER DRAWS
    cons = draw.B(          end      ,ii);
    Bdy  = draw.B(n_lag_z+1:end-1    ,ii)';
    Bz   = draw.B(        1:n_lag_z,ii)';
    vy   = draw.V(ii);
    
    % CONSTRUCT NECESSARY OBJECTS
    By = [Bdy(1)+1, Bdy(2:end)-Bdy(1:end-1), -Bdy(end)];
    By_comp = [By; eye(n_lag_y),  zeros(n_lag_y,1)];	% companion form coefficient
    Bz_comp = [Bz; zeros(n_lag_y, n_lag_z)];

    % PRE-COMPUTE VARIANCES
    VY = zeros(n_lag_y+1);        % matrix for conditional variance
    VZ = zeros(n_lag_z);
    for hh = 1:H+1
        VZ(hh,hh) = vz;
        VY = By_comp*VY*By_comp' + Bz_comp*VZ*Bz_comp';
        VY(1,1) = VY(1,1) + vy;
        VY_draw(hh,:,1,ii) = VY(1,1);      % variance of interest
    end
    
    % COMPUTE EXPECTATIONS
    for tt = 1:T_test
        for hh = 1:H+1
            if test(tt)+hh-1 <= T
            if hh == 1 
                % HORIZON 0: CONDITION ON LAGGED DATA AND SHOCK
                EY_prev = Y_lag(test(tt),:)';       % previous observations
                EZ      = Z_lag(test(tt),:)';
                EY      = NaN(1+n_lag_y,1);         % vector for conditional expectation
            else
                % HORIZON>0: ITERATE FORWARD
                EY_prev   = EY;
                EZ(2:end) = EZ(1:end-1);    % lagged shocks
                EZ(1)     = mz;             % use sample mean
            end
            EY(1)     = cons + By*EY_prev + Bz*EZ;
            EY(2:end) = EY_prev(1:end-1);   % lagged entries from previous expectation
            EY_draw(hh,tt,1,ii) = EY(1);    % expectation of interest
            end
        end
    end
    
end

% PREDICTIVE DENSITY
for hh = 1:H+1
    tmp = test+hh-1 <= T;
    test_sub = test(tmp);
    pd(hh,tmp,1,:) = normpdf(Y(test_sub+hh-1),...
        squeeze(EY_draw(hh,tmp,1,:)), squeeze(sqrt(VY_draw(hh,tmp,1,:))));
end
pd = mean(pd,4);

% EY = mean(EY_draw,3); %AG added. Think these lines can be deleted (AG
% 4/11/2022)
% VY = mean(VY_draw,3); %AG added so that the checking function can work

% SAVE OUTPUT
% dimensions: T x n_ind x H x n_draw
pd_out.pd      = permute(     pd,[2,3,1]  );
% pd_out.EY_draw = permute(EY_draw,[2,3,1,4]);
% pd_out.VY_draw = permute(VY_draw,[2,3,1,4]);

end