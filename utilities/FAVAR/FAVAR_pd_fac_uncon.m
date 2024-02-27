function F_draw = FAVAR_pd_fac_uncon(data, draw_F, test, H, identification)

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
    F = [data.m, data.fac]; % ind = ind + 1; IMPORTANT: We are assuming the ind fed in is already accounting for being shock augmented 
elseif strcmp(identification,'Cholesky') || strcmp(identification,'proxy') 
    F = data.fac;
end

% DIMENSIONS
[T, N] = size(F);           % number of periods, variables
n_draw = size(draw_F.B,3);	% number of draws
n_lag  = size(draw_F.B,2)/N;  % number of lags
T_test = length(test);      % number of observations in test set

% LAGGED F
X = lagmatrix(F, 1:n_lag);

% MATRICES FOR STORAGE
EY_draw = NaN(H+1, T_test, N, n_draw);
VY_draw = NaN(H+1, N, N, n_draw);
% F_draw  = NaN(H+1, T_test, N, n_draw);


% ------------------------------------------------------------------------
% PREDICTIVE DENSITIES
% ------------------------------------------------------------------------

for ii = 1:n_draw
    
    % EXTRACT PARAMETER DRAWS
    cons = draw_F.cons(:,ii);
    B = draw_F.B(:,:,ii);
    V = draw_F.V(:,:,ii);
    
    % CONSTRUCT NECESSARY OBJECTS
    B_comp = [B; eye(N*(n_lag-1)),  zeros(N*(n_lag-1),N)];	% companion form coefficient
    
    % PRE-COMPUTE VARIANCES
    VY = zeros(N*n_lag);        % matrix for conditional variance
    for hh = 1:H+1
        VY = B_comp*VY*B_comp';
        VY(1:N,1:N) = VY(1:N,1:N) + V;
        VY_draw(hh,:,:,ii) = VY(1:N,1:N);
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
            EY_draw(hh,tt,:,ii) = EY(1:N);
            end
        end
    end

end


% ------------------------------------------------------------------------
% DRAW FACTORS FROM PREDICTIVE DENSITY
% ------------------------------------------------------------------------
% ignore memory requirements
MM = reshape(...
    permute(EY_draw,[1,4,2,3]), ...
    (H+1)*n_draw*T_test, N);
VV = repmat(reshape(...
    permute(VY_draw,[2,3,1,4]), N, N, (H+1)*n_draw),...
    1, 1, T_test);
F_draw = mvnrnd(MM, VV)';
F_draw = permute(reshape(F_draw,N,H+1,n_draw,T_test),[2,4,1,3]);
% for tt = 1:T_test
%     for hh = 1:H+1
%         F_draw(hh,tt,:,:) = mvnrnd(squeeze(EY_draw(hh,tt,:,:))',squeeze(VY_draw(hh,:,:,:)))';
%     end
% end


end