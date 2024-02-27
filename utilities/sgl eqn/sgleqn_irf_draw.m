function irf_draws = sgleqn_irf_draw(par_draws, n_lag_y, n_lag_z)

% B_draw = mvnrnd(out.B, out.V, n_draw);
n_draw = size(par_draws.B,2);
irf_draws = NaN(n_lag_z+1, n_draw);
for i = 1:n_draw
    irf_draws(:,i) = sgleqn_irf_pt(par_draws.B(:,i), n_lag_y, n_lag_z);
end

end

function psi = sgleqn_irf_pt(B, n_lag_y, n_lag_z)

psi = NaN(n_lag_z+1,1); % irf for growth rate
psi(1) = B(1);
for h = 1:n_lag_z
    K = min(h, n_lag_y);
    psi(h+1) = B(h+1) + ...
        sum(B(h:-1:h-K+1).*B(end-n_lag_y:end-n_lag_y+K-1));
end
psi = cumsum(psi); % irf for levels

end