function irf = VAR_irf(draw, hor, ind, identification, data)

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
%   irf(i,j,t,k): impulse response of variable i to shock j at horizon t-1, draw k.
%
% ------------------------------------------------------------------------

[N, n_draw] = size(draw.V,[2,3]);
n_shocks = length(ind);
irf = NaN(N, n_shocks, hor+1, n_draw);
for i = 1:n_draw
    % IDENTIFICATION
    if strcmp(identification,'internal IV')
        C = chol(draw.V(:,:,i),'lower');
    elseif strcmp(identification,'proxy')
        out_tmp = draw.est_out;
        out_tmp.B = draw.B(:,:,i); out_tmp.V = draw.V(:,:,i);
        C = VAR_proxy(data, out_tmp);
    elseif strcmp(identification,'Cholesky')
        C = chol(draw.V(:,:,i),'lower');
    end
    
    % COMPUTE IRF
    irf(:,:,:,i) = VAR_irf_point(draw.B(:,:,i), C, hor, ind);
end


end