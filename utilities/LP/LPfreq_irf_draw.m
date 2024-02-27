function draws = LPfreq_irf_draw(out, n_draw)

H = length(out);
draws = NaN(H, n_draw);

for h = 1:H
    draws(h,:) = normrnd(out{h}.B(1), sqrt(out{h}.V(1,1)), n_draw, 1);
%     shocks = normrnd(0, sqrt(out{h}.V(1,1)), n_draw, 1 );
%     draws(h,:) = shocks + out{h}.B(1);
end

end
