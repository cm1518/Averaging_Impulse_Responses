function irf_qtl = LPfreq_irf(out, qtl)

H = length(out);
n_qtl = length(qtl);
z_qtl = norminv(qtl);

irf_qtl = NaN(H, n_qtl);
for h = 1:H
    irf_qtl(h,:) = z_qtl*sqrt(out{h}.V(1,1)) + out{h}.B(1);
end

end