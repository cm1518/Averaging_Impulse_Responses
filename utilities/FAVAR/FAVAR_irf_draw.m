function irf_Y_draw = FAVAR_irf_draw(irf_F_draw, out_Y)

[n_fac, n_sh, H, n_draw] = size(irf_F_draw);
B_draw = mvnrnd(out_Y.B', out_Y.V, n_draw);
irf_Y_draw = NaN(n_sh, H, n_draw);
for ii = 1:n_draw
    irf_Y_draw(:,:,ii) = sum(irf_F_draw(:,:,:,ii).*B_draw(ii,1:n_fac)');
end

end