function draw = VARbayes_draw(post, n_draw)


% ------------------------------------------------------------------------
% LOAD ESTIMATES
% ------------------------------------------------------------------------

[NP, N] = size(post.B_cons);
b = reshape(post.B_cons,NP*N,1);
OMEGA = post.OMEGA;
PSI = post.PSI;
df = post.df;


% ------------------------------------------------------------------------
% MAKE DRAWS
% ------------------------------------------------------------------------

% VARIANCE
SIGMA_draw = wishrnd_bart(inv(PSI), df, n_draw);
for i = 1:n_draw
    SIGMA_draw(:,:,i) = eye(N)/SIGMA_draw(:,:,i);
end

% COEFFICIENTS
SIGMA_OMEGA_draw = kron(reshape(SIGMA_draw, N, N*n_draw), OMEGA);
SIGMA_OMEGA_draw = reshape(SIGMA_OMEGA_draw, N*NP, N*NP, n_draw);
b_draw = mvnrnd(b, SIGMA_OMEGA_draw)';

% ------------------------------------------------------------------------
% SAVE DRAWS
% ------------------------------------------------------------------------

B_draw = reshape(b_draw, NP, N, n_draw);
draw.cons = squeeze(B_draw(end,:,:));
draw.B    = permute(B_draw(1:end-1,:,:), [2,1,3]);
draw.V    = SIGMA_draw;
draw.est_out = post;


end