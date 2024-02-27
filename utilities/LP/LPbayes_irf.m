function irf = LPbayes_irf(draw, ind)

[N, n_draw] = size(draw{1}.V, [2,3]);
n_lag = size(draw{1}.B,2)/N;
H = length(draw);
n_shocks = length(ind);

irf = NaN(N, n_shocks, H, n_draw);
for i = 1:n_draw
    for j = 1:n_shocks
        i_shock = ind(j);
        for h = 1:H
            if h == 1
                C_comp = zeros(N*n_lag);
                C_comp(1:N,1:N) = chol(draw{h}.V(:,:,i),'lower');
                y = C_comp(:,i_shock);
            else
                B_comp = [draw{h}.B(:,:,i); ...
                    eye(N*(n_lag-1)),  zeros(N*(n_lag-1),N)];  
                y = B_comp*C_comp(:,i_shock);
            end
            irf(:,j,h,i) = y(1:N);
        end
   end    
end

end