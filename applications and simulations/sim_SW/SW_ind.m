function ii = SW_ind

% ========================================================================
%
% ORGANIZE VARIABLES
%
% ========================================================================
% name indices for each variable

% ------------------------------------------------------------------------
% ENDOGENOUS VARIABLES (y(t))
% ------------------------------------------------------------------------

% < RIGID PRICE VARIABLES >
ii.y   = 1;      % output
ii.c   = 2;      % consumption
ii.i   = 3;      % investment
ii.q   = 4;      % real value of capital stock
ii.l   = 5;      % labor
ii.ks  = 6;      % capital services used in production
ii.k   = 7;      % capital stock
ii.z   = 8;      % capital utilization
ii.mup = 9;      % price markup
ii.pii = 10;     % inflation
ii.rk  = 11;     % rental rate of capital
ii.muw = 12;     % wage markup
ii.w   = 13;     % wage
ii.r   = 14;     % interest rate

% < RIGID PRICE EXPECTATIONS >
ii.E.c   = 15;
ii.E.q   = 16;
ii.E.i   = 17;
ii.E.pii = 18;
ii.E.l   = 19;
ii.E.rk  = 20;
ii.E.w   = 21;

% < FLEXIBLE PRICE VARIABLES >
ii.y_f  = 22;
ii.c_f  = 23;
ii.i_f  = 24;
ii.q_f  = 25;
ii.l_f  = 26;
ii.ks_f = 27;
ii.z_f  = 28;
ii.k_f  = 29;
ii.w_f  = 30;
ii.rk_f = 31;
ii.r_f  = 32;

% < FLEXIBLE PRICE EXPECTATIONS >
ii.E.c_f  = 33;
ii.E.q_f  = 34;
ii.E.i_f  = 35;
ii.E.l_f  = 36;
ii.E.rk_f = 37;

% < EXOGENOUS PROCESSES >
ii.eps.a = 38;    % TFP
ii.eps.b = 39;    % net worth
ii.eps.g = 40;    % exogenous spending
ii.eps.i = 41;    % investment-specific
ii.eps.r = 42;    % interest rate
ii.eps.p = 43;    % price markup
ii.eps.w = 44;    % wage markup

% < CURRENT INNOVATIONS >
ii.eta.p  = 45;   % price markup
ii.eta.w  = 46;   % wage markup


% ------------------------------------------------------------------------
% EXOGENOUS VARIABLES (z(t))
% ------------------------------------------------------------------------
ii.eta_sh.a = 1;    % TFP
ii.eta_sh.b = 2;    % net worth
ii.eta_sh.g = 3;    % exogenous spending
ii.eta_sh.i = 4;    % investment-specific
ii.eta_sh.r = 5;    % interest rate
ii.eta_sh.p = 6;    % price markup
ii.eta_sh.w = 7;    % wage markup


% ------------------------------------------------------------------------
% EXPECTATION ERRORS (eta(t))
% ------------------------------------------------------------------------
ii.E_sh.c    = 1;
ii.E_sh.q    = 2;
ii.E_sh.i    = 3;
ii.E_sh.pii  = 4;
ii.E_sh.l    = 5;
ii.E_sh.rk   = 6;
ii.E_sh.w    = 7;
ii.E_sh.c_f  = 8;
ii.E_sh.q_f  = 9;
ii.E_sh.i_f  = 10;
ii.E_sh.l_f  = 11;
ii.E_sh.rk_f = 12;

end