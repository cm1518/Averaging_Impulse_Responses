function draw = LPbayes_draw(post, n_draw)

H = size(post,1) - 1;
draw = cell(H+1,1);

for h = 0:H
    
    % --------------------------------------------------------------------
    % LOAD ESTIMATES
    % --------------------------------------------------------------------
    
    [NP, N] = size(post{h+1}.B_cons);
    b = reshape(post{h+1}.B_cons,NP*N,1);
    OMEGA = post{h+1}.OMEGA;
    PSI = post{h+1}.PSI;
    df = post{h+1}.df;


    % --------------------------------------------------------------------
    % MAKE DRAWS
    % --------------------------------------------------------------------

    % VARIANCE
    SIGMA_draw = wishrnd_bart(inv(PSI), df, n_draw);
    for i = 1:n_draw
        SIGMA_draw(:,:,i) = eye(N)/SIGMA_draw(:,:,i);
    end

    % COEFFICIENTS
    SIGMA_OMEGA_draw = kron(reshape(SIGMA_draw, N, N*n_draw), OMEGA);
    SIGMA_OMEGA_draw = reshape(SIGMA_OMEGA_draw, N*NP, N*NP, n_draw);
    b_draw = mvnrnd(b, SIGMA_OMEGA_draw)';

    % --------------------------------------------------------------------
    % SAVE DRAWS
    % --------------------------------------------------------------------

    B_draw = reshape(b_draw, NP, N, n_draw);
    draw{h+1}.cons = squeeze(B_draw(end,:,:));
    draw{h+1}.B    = permute(B_draw(1:end-1,:,:), [2,1,3]);
    draw{h+1}.V    = SIGMA_draw;
    draw{h+1}.est_out = post{h+1};

end

end