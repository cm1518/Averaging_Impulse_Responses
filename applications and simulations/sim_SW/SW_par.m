function pp = SW_par

% ========================================================================
%
% EXTRACT PARAMETERS
% take mode from Smets and Wouters
%
% ========================================================================

% ------------------------------------------------------------------------
% ENDOGENOUS
% ------------------------------------------------------------------------
pp.varphi = 5.48;	% SS elasticity of capital adjustment cost function
pp.sig_c  = 1.39;  % inverse of elasticity of intertemporal substitution
pp.h      = 0.71;  % habit formation
pp.xi_w   = 0.73;  % wage rigidity (0 for flexible wages)
pp.sig_l  = 1.92;  % elasticity of labor supply wrt real wage
pp.xi_p   = 0.65;  % price rigidity (0 for flexible prices)
pp.iota_w = 0.59;  % wage indexation
pp.iota_p = 0.22;  % price indexation
pp.psi    = 0.54;  % positive func of elasticity of capital util adjust cost (in [0,1])
pp.Phi    = 1.61;  % 1 + share of fixed costs in production
pp.r_pii  = 2.03;  % Taylor rule coeff for inflation
pp.rho    = 0.81;  % Taylor rule persistence
pp.r_y    = 0.08;  % Taylor rule coeff for output gap
pp.r_Dy   = 0.22;  % Taylor rule coeff for change in output gap

% ------------------------------------------------------------------------
% STEADY STATE
% ------------------------------------------------------------------------
pp.pii_bar = 0.81; % quarterly inflation
pp.bet_inv = 0.16; % 100*(beta^{-1} -1)
pp.l_bar   = -0.1; % hours
pp.gam_bar = 0.43; % quarterly GDP growth
pp.alp     = 0.19;	% capital share

% ------------------------------------------------------------------------
% EXOGENOUS
% ------------------------------------------------------------------------
% < PERSISTENCE >	< VARIANCE >
pp.rho_a  = 0.95;   pp.sig_a = 0.45;    % TFP
pp.rho_b  = 0.18;   pp.sig_b = 0.24;    % net worth
pp.rho_g  = 0.97;   pp.sig_g = 0.52;    % exogenous spending
pp.rho_i  = 0.71;   pp.sig_i = 0.45;    % investment-specific
pp.rho_r  = 0.12;   pp.sig_r = 0.24;    % interest rate
pp.rho_p  = 0.90;   pp.sig_p = 0.14;    % price markup
pp.rho_w  = 0.97;   pp.sig_w = 0.24;    % wage markup
pp.mu_p   = 0.74;                       % minus price markup lagged innovation
pp.mu_w   = 0.88;                       % minus wage markup lagged innovation
pp.rho_ga = 0.52;                       % exogenous spending exposure to TFP

% ------------------------------------------------------------------------
% CALIBRATED PARAMETERS
% ------------------------------------------------------------------------
pp.del      = 0.025;   % depreciation rate
pp.gy       = 0.18;    % spending-GDP ratio
pp.lambda_w = 1.5;     % steady-state markup in labor market
pp.eps_p    = 10;      % Kimball aggregator for goods market
pp.eps_w    = 10;      % Kimball aggregator for labor market

% ------------------------------------------------------------------------
% IMPLIED PARAMETERS
% ------------------------------------------------------------------------
pp.bet       = 1/(pp.bet_inv/100 + 1);
pp.gam       = pp.gam_bar/100 + 1;
pp.pii_st    = pp.pii_bar/100 + 1;
pp.gam_sig_c = pp.gam^pp.sig_c;
pp.r_bar     = 100*(pp.gam_sig_c/pp.bet*pp.pii_st - 1);
pp.rk_ss     = pp.gam_sig_c/pp.bet - 1 + pp.del;
pp.w_ss      = ( pp.alp^pp.alp * (1-pp.alp)^(1-pp.alp) /...
                (pp.Phi*pp.rk_ss^pp.alp) )^(1/(1-pp.alp));
pp.ik        = pp.gam - 1 + pp.del;
pp.lk        = (1-pp.alp) / pp.alp * pp.rk_ss / pp.w_ss;
pp.ky        = pp.Phi*pp.lk^(pp.alp-1);
pp.iy        = (pp.gam-1+pp.del)*pp.ky;
pp.cy        = 1 - pp.gy - pp.iy;
pp.zy        = pp.rk_ss*pp.ky;
pp.wlc       = (1-pp.alp)/pp.alp * pp.rk_ss * pp.ky/pp.cy/pp.lambda_w;

end