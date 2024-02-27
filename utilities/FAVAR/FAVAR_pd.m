function F_draw = FAVAR_pd(data, draw_F, test, H, identification)


% ========================================================================
%
%                  DRAW FACTORS FROM PREDICTIVE DENSITY
%
% ========================================================================

% ------------------------------------------------------------------------
% PRELIMINARIES
% ------------------------------------------------------------------------

% LOAD DATA
% Y = data.Y;
if strcmp(identification,'internal IV')
    F = [data.m, data.F]; % ind = ind + 1; IMPORTANT: We are assuming the ind fed in is already accounting for being shock augmented 
elseif strcmp(identification,'Cholesky') || strcmp(identification,'proxy') 
    F = data.F;
end

% DIMENSIONS
[T, N] = size(F);           % number of periods, variables
n_draw = size(draw_F.B,3);	% number of draws
n_lag  = size(draw_F.B,2)/N;  % number of lags
T_test = length(test);      % number of observations in test set

% LAGGED Y
X = lagmatrix(F, 1:n_lag);

% MATRICES FOR STORAGE
% pd      = NaN(H+1, T_test, length(ind), n_draw);
EY_draw = NaN(H+1, T_test, N, n_draw);
VY_draw = NaN(H+1, N, N, n_draw);
F_draw  = NaN(H+1, T_test, N, n_draw);


% ------------------------------------------------------------------------
% PREDICTIVE DENSITIES
% expectations and variances are H x T_test x n_ind x n_draw
% ------------------------------------------------------------------------

for ii = 1:n_draw
    
    % EXTRACT PARAMETER DRAWS
    cons = draw_F.cons(:,ii);
    B = draw_F.B(:,:,ii);
    V = draw_F.V(:,:,ii);
    
    % CONSTRUCT NECESSARY OBJECTS
    B_comp = [B; eye(N*(n_lag-1)),  zeros(N*(n_lag-1),N)];	% companion form coefficient
    % *** NEW - ALLOW FOR CHOLESKY **** %
    if strcmp(identification,'Cholesky')
        C = chol(V,'lower');
        C_N = C(:,end);
        C_i = eye(N)/C;
    elseif strcmp(identification,'proxy')
        out_tmp = draw_F.est_out;
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
        VY_draw(hh,:,:,ii) = VY(1:N,1:N);
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
                    EY(1) = F(test(tt),1);
                    EY(2:N) = cons(2:N) + B(2:N,:)*EY_prev + ...
                        V(2:N,1)/V(1,1)*(F(test(tt),1) - cons(1) - B(1,:)*EY_prev);
                elseif strcmp(identification,'Cholesky')
                    uu = F(test(tt)) - cons - B*EY_prev;     % reduced-form residual
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
            EY_draw(hh,tt,:,ii) = EY(1:N);
            end
        end
    end

    % DRAW FACTORS FROM PREDICTIVE DENSITY
    for tt = 1:T_test
        for hh = 1:H+1
            F_draw(hh,tt,:,:) = mvnrnd(squeeze(EY_draw(hh,tt,:,:))',squeeze(VY_draw(hh,:,:,:)));
        end
    end
    
end


end