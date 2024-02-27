function pd_out = VAR_pd(data, draw, test, ind, H, identification)

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
%       >> 'internal IV': Cholesky, instrument ordered first (Plagborg-Moller and Wolf, 2021)
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
if strcmp(identification,'internal IV')
    Y = [data.m, data.Y]; % ind = ind + 1; IMPORTANT: We are assuming the ind fed in is already accounting for being shock augmented 
% *** NEW - ALLOW FOR CHOLESKY **** %
elseif strcmp(identification,'Cholesky') || strcmp(identification,'proxy') 
    Y = data.Y;
end

% DIMENSIONS
[T, N] = size(Y);           % number of periods, variables
n_draw = size(draw.B,3);	% number of draws
n_lag  = size(draw.B,2)/N;  % number of lags
T_test = length(test);      % number of observations in test set

% LAGGED Y
X = lagmatrix(Y, 1:n_lag);

% MATRICES FOR STORAGE
pd      = NaN(H+1, T_test, length(ind), n_draw);
EY_draw = NaN(H+1, T_test, length(ind), n_draw);
VY_draw = NaN(H+1, T_test, length(ind), n_draw);


% ------------------------------------------------------------------------
% PREDICTIVE DENSITIES
% expectations and variances are H x T_test x n_ind x n_draw
% ------------------------------------------------------------------------

for ii = 1:n_draw
    
    % EXTRACT PARAMETER DRAWS
    cons = draw.cons(:,ii);
    B = draw.B(:,:,ii);
    V = draw.V(:,:,ii);
    
    % CONSTRUCT NECESSARY OBJECTS
    B_comp = [B; eye(N*(n_lag-1)),  zeros(N*(n_lag-1),N)];	% companion form coefficient
    % *** NEW - ALLOW FOR CHOLESKY **** %
    if strcmp(identification,'Cholesky')
        C = chol(V,'lower');
        C_N = C(:,end);
        C_i = eye(N)/C;
    elseif strcmp(identification,'proxy')
        out_tmp = draw.est_out;
        out_tmp.B = B; out_tmp.V = V;
        C = VAR_proxy(data, out_tmp);
        C1 = C(:,1);
    end

    % PRE-COMPUTE VARIANCES
    VY = zeros(N*n_lag);        % matrix for conditional variance
    for hh = 1:H+1
        if hh == 1
            if strcmp(identification,'internal IV')
                VY(2:N,2:N) = V(2:N,2:N) - V(2:N,1)/V(1,1)*V(1,2:N);
            elseif strcmp(identification,'Cholesky')
                VY(1:N,1:N) = V - C_N*C_N';
            elseif strcmp(identification,'proxy')
                VY(1:N,1:N) = V - C1*C1';                                
            end
        else
            VY = B_comp*VY*B_comp';
            VY(1:N,1:N) = VY(1:N,1:N) + V;
        end
        for jj = 1:length(ind)
            VY_draw(hh,:,jj,ii) = VY(ind(jj),ind(jj));      % variance of interest
        end
    end
    
    % COMPUTE EXPECTATIONS
    for tt = 1:T_test        
        for hh = 1:H+1
            if test(tt)+hh-1 <= T
            if hh == 1 
                % HORIZON 0: CONDITION ON LAGGED DATA AND SHOCK
                EY_prev = X(test(tt),:)';    % previous observations
                EY = NaN(N*n_lag,1);        % vector for conditional expectation
                if strcmp(identification,'internal IV')
                    EY(1) = Y(test(tt),1);
                    EY(2:N) = cons(2:N) + B(2:N,:)*EY_prev + ...
                        V(2:N,1)/V(1,1)*(Y(test(tt),1) - cons(1) - B(1,:)*EY_prev);
                elseif strcmp(identification,'Cholesky')
                    uu = Y(test(tt)) - cons - B*EY_prev;     % reduced-form residual
                    ee = C_i*uu;                             % structural shock
                    EY(1:N) = cons + B*EY_prev + C_N*ee(N);
                elseif strcmp(identification,'proxy')
                    EY(1:N) = cons + B*EY_prev + C1*data.m(test(tt));                         
                end
            else
                % HORIZON >0: ITERATE FORWARD
                EY_prev = EY;
                EY(1:N) = cons + B*EY_prev;
            end
            EY(N+1:end) = EY_prev(1:N*(n_lag-1));   % lagged entries from previous expectation
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

% EY = mean(EY_draw,3); %AG added. Think these lines can be deleted (AG
% 4/11/2022)
% VY = mean(VY_draw,3); %AG added so that the checking function can work

% SAVE OUTPUT
% dimensions: T x n_ind x H x n_draw
pd_out.pd      = permute(     pd,[2,3,1]  );
% pd_out.EY_draw = permute(EY_draw,[2,3,1,4]);
% pd_out.VY_draw = permute(VY_draw,[2,3,1,4]);

end