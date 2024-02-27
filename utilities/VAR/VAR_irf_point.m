function irf = VAR_irf_point(B, C, hor, ind)

% ========================================================================
% COMPUTE VAR IMPULSE RESPONSE FOR GIVEN COEFFICIENTS
% ========================================================================
%
% -- MODEL --
%   y(t) = con + B*[y(t-1),...,y(t-n_lag)] + C*u(t)
%   var(u(t)) = I
%
% -- INPUT --
%   B(N*nlag-by-N): VAR coefficient
%   C(N-by-N): effect of structural shock on impact
%   hor : horizon for impulse response
%   ind : indices for shocks of interest
% 
% -- OUTPUT --
%   irf(i,j,t): impulse response of variable i to shock j at horizon t-1.
%
% ------------------------------------------------------------------------


% CONSTRUCT COMPANION FORM
N = size(C, 1);
n_lag = size(B,2)/N;
B_comp = [B; eye(N*(n_lag-1)),  zeros(N*(n_lag-1),N)];
C_comp = zeros(N*n_lag); C_comp(1:N,1:N) = C;

% % SHOCKS TO COMPUTE IMPULSE RESPONSE FOR
% if nargin == 2
%     n_shocks = size(C,2);	 % number of shocks to compute impulse response for
%     ind = 1:n_shocks;        % compute irf for all shocks if none selected
% else
%     n_shocks = length(ind);
% end

% IMPULSE RESPONSES
n_shocks = length(ind);
irf = NaN(N, n_shocks, hor+1);
for j = 1:n_shocks
    i_shock = ind(j);
    for h = 1:hor+1
        if h == 1
            y = C_comp(:,i_shock);
        else
            y = B_comp*y;
        end
        irf(:,j,h) = y(1:N);
    end
end


end