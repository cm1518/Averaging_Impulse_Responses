function [retcode, T1, T0, Se, pp, ii] = SW_sol

% ========================================================================
%
%           COMPUTE LIKLELIHOOD FOR  SMETS & WOUTERS (2007) MODEL
% ------------------------------------------------------------------------
%
% A. SOLVE MODEL USING GENSYS
%       1. cast Smets & Wouters (2007) model in the form
%           g0*y(t) = g1*y(t-1) + c + psi*z(t) + pi*eta(t)
%       2. use gensys to generate system
%           y(t) = G1*y(t-1) + C + impact*z(t)
%
% B. CAST MODEL INTO STATE SPACE REPRESENTATION
% 
% C. COMPUTE LIKELIHOOD
% 
% ========================================================================


% EXTRACT PARAMETERS
pp = SW_par;

% INDEX VARIABLES
ii = SW_ind;


% ========================================================================
%
%                       EQUILIBRIUM CONDITIONS
%
% ========================================================================
% prepare matrices to input into gensys
% GAM0*y(t) = GAM1*y(t-1) + C + PSI*z(t) + PII*eta(t)

% ------------------------------------------------------------------------
% INITIALIZE MATRICES
% ------------------------------------------------------------------------
GAM0 = zeros(46);           GAM1 = zeros(46);           C = zeros(46,1);
PSI  = zeros(46,7);         PII  = zeros(46,12);


% ------------------------------------------------------------------------
% RIGID PRICE EQUATIONS
% ------------------------------------------------------------------------

% 1. RESOURCE CONSTRAINT
GAM0(1,ii.y) = 1;       GAM0(1,ii.c) = -pp.cy;
GAM0(1,ii.i) = -pp.iy;  GAM0(1,ii.z) = -pp.zy;  GAM0(1,ii.eps.g) = -1;

% 2. CONSUMPTION EULER EQUATION
hg = pp.h/pp.gam;
GAM0(2,ii.c)     = 1;               GAM1(2,ii.c) = hg/(1+hg);
GAM0(2,ii.E.c)   = -1/(1+hg);       GAM0(2,ii.l) = -pp.wlc*(pp.sig_c-1) / ( pp.sig_c*(1+hg) );
GAM0(2,ii.E.l)   = -GAM0(2,ii.l);   GAM0(2,[ii.r,ii.eps.b]) = (1-hg) / ( pp.sig_c*(1+hg) );
GAM0(2,ii.E.pii) = -GAM0(2,ii.r);

% 3. INVESTMENT EULER EQUATION
bg_sig_c = pp.bet*pp.gam^(1-pp.sig_c);
GAM0(3,ii.i)     = 1;                         GAM1(3,ii.i) = 1/(1+bg_sig_c);
GAM0(3,ii.E.i)   = -bg_sig_c/(1+bg_sig_c);    GAM0(3,ii.q) = -1/( pp.varphi*pp.gam^2*(1+bg_sig_c) );
GAM0(3,ii.eps.i) = -1;

% 4. VALUE OF CAPITAL
GAM0(4,[ii.q,ii.r,ii.eps.b]) = 1;
GAM0(4,ii.E.q) = -(1-pp.del)*pp.bet*pp.gam^(-pp.sig_c);
GAM0(4,ii.E.pii) = -1;
GAM0(4,ii.E.rk) = -(1+GAM0(4,ii.E.q));

% 5. PRODUCTION FUNCTION
GAM0(5,ii.y) = 1;                   GAM0(5,ii.ks)    = -pp.Phi*pp.alp;
GAM0(5,ii.l) = -pp.Phi*(1-pp.alp);  GAM0(5,ii.eps.a) = -pp.Phi;

% 6. CAPITAL UTILIZATION DEFINITION
GAM0(6,ii.ks) = 1;      GAM1(6,ii.k) = 1;       GAM0(6,ii.z)  = -1;

% 7. OPTIMAL CAPITAL UTILIZATION
GAM0(7,ii.z) = 1;       GAM0(7,ii.rk) = -(1-pp.psi)/pp.psi;

% 8. EVOLUTION OF CAPITAL
GAM0(8,ii.k) = 1;                   GAM1(8,ii.k)     = (1-pp.del)/pp.gam;
GAM0(8,ii.i) = -(1-GAM1(8,ii.k));   GAM0(8,ii.eps.i) = GAM0(8,ii.i)*pp.varphi*pp.gam^2*(1+bg_sig_c);

% 9. PRICE MARKUP
GAM0(9,[ii.mup,ii.w]) = 1;          GAM0(9,ii.ks) = -pp.alp;
GAM0(9,ii.l)          = pp.alp;     GAM0(9,ii.eps.a) = -1;

% 10. PHILLIPS CURVE
GAM0(10,ii.pii) = 1;                         GAM0(10,ii.E.pii) = -bg_sig_c/(1 + pp.iota_p*bg_sig_c);
GAM1(10,ii.pii) = pp.iota_p/(1 + bg_sig_c);  GAM0(10,ii.eps.p) = -1;
GAM0(10,ii.mup) = (1 - bg_sig_c*pp.xi_p)*(1 - pp.xi_p) / ...
    ( (1 + pp.iota_p*bg_sig_c)*(1 + (pp.Phi-1)*pp.eps_p)*pp.xi_p );

% 11. RENTAL RATE OF CAPITAL
GAM0(11,[ii.rk,ii.ks]) = 1;     GAM0(11,[ii.l,ii.w])  = -1;

% 12. MARGINAL SUBSTITUTION
GAM0(12,ii.muw) = 1;            GAM0(12,ii.w) = -1;
GAM0(12,ii.l)   = pp.sig_l;     GAM0(12,ii.c) = 1/(1-hg);
GAM1(12,ii.c)   = hg/(1-hg);

% 13. EVOLUTION OF WAGES
GAM0(13,ii.w)   = 1;                GAM0(13,[ii.E.w,ii.E.pii]) = -bg_sig_c/(1+bg_sig_c);
GAM1(13,ii.w)   = 1/(1+bg_sig_c);   GAM1(13,ii.pii) = -pp.iota_w*GAM1(13,ii.w);
GAM0(13,ii.pii) = (1+bg_sig_c*pp.iota_w)/(1+bg_sig_c);
GAM0(13,ii.muw) = (1-bg_sig_c*pp.xi_w)*(1-pp.xi_w) / ( (1+bg_sig_c)*(1+(pp.lambda_w-1)*pp.eps_w)*pp.xi_w );
GAM0(13,ii.eps.w) = -1;

% 14. MONETARY POLICY RULE
GAM0(14,ii.r)   = 1;
GAM1(14,ii.r) = pp.rho;
GAM0(14,ii.pii) = -(1-pp.rho)*pp.r_pii; 
GAM0(14,ii.y)   = -(1-pp.rho)*pp.r_y - pp.r_Dy;
GAM0(14,ii.y_f) = -GAM0(14,ii.y);       GAM1(14,ii.y) = -pp.r_Dy;
GAM1(14,ii.y_f) = -GAM1(14,ii.y);       GAM0(14,ii.eps.r) = -1;


% ------------------------------------------------------------------------
% RIGID PRICE EXPECTATIONS
% ------------------------------------------------------------------------

% 15. CONSUMPTION (c)
GAM0(15,ii.c) = 1;      GAM1(15,ii.E.c) = 1;    PII(15,ii.E_sh.c) = 1;

% 16. REAL VALUE OF CAPITAL STOCK (q)
GAM0(16,ii.q) = 1;      GAM1(16,ii.E.q) = 1;	PII(16,ii.E_sh.q) = 1;

% 17. INVESTMENT (i)
GAM0(17,ii.i) = 1;      GAM1(17,ii.E.i) = 1;    PII(17,ii.E_sh.i) = 1;

% 18. INFLATION (pii)
GAM0(18,ii.pii) = 1;	GAM1(18,ii.E.pii) = 1;  PII(18,ii.E_sh.pii) = 1;

% 19. LABOR (l)
GAM0(19,ii.l) = 1;      GAM1(19,ii.E.l) = 1;    PII(19,ii.E_sh.l) = 1;

% 20. RENTAL RATE OF CAPITAL (rk)
GAM0(20,ii.rk) = 1;     GAM1(20,ii.E.rk) = 1;   PII(20,ii.E_sh.rk) = 1;

% 21. WAGES (w)
GAM0(21,ii.w) = 1;      GAM1(21,ii.E.w) = 1;    PII(21,ii.E_sh.w) = 1;


% ------------------------------------------------------------------------
% FLEXIBLE PRICE EQUATIONS
% ------------------------------------------------------------------------

% 22. RESOURCE CONSTRAINT
GAM0(22,ii.y_f) = 1;            GAM0(22,ii.c_f) = -pp.cy;
GAM0(22,ii.i_f) = -pp.iy;       GAM0(22,ii.z_f) = -pp.zy;
GAM0(22,ii.eps.g) = -1;

% 23. CONSUMPTION EULER EQUATION
GAM0(23,ii.c_f) = 1;            GAM1(23,ii.c_f) = hg/(1+hg);
GAM0(23,ii.E.c_f) = -1/(1+hg);	GAM0(23,ii.l_f) = -pp.wlc*(pp.sig_c-1) / ( pp.sig_c*(1+hg) );
GAM0(23,ii.E.l_f) = -GAM0(23,ii.l_f);
GAM0(23,[ii.r_f,ii.eps.b]) = (1-hg) / ( pp.sig_c*(1+hg) );
% GAM0(23,ii.eps.b) = GAM0(23,ii.r_f);

% 24. INVESTMENT EULER EQUATION
GAM0(24,ii.i_f)   = 1;                      GAM1(24,ii.i_f) = 1/(1+bg_sig_c);
GAM0(24,ii.E.i_f) = -bg_sig_c/(1+bg_sig_c); GAM0(24,ii.q_f) = -1/( pp.varphi*pp.gam^2*(1+bg_sig_c) );
GAM0(24,ii.eps.i) = -1;

% 25. VALUE OF CAPITAL
GAM0(25,[ii.q_f,ii.r_f,ii.eps.b]) = 1;
GAM0(25,ii.E.q_f) = -(1-pp.del)*pp.bet*pp.gam^(-pp.sig_c);
GAM0(25,ii.E.rk_f) = -(1+GAM0(25,ii.E.q_f)); 

% 26. PRODUCTION FUNCTION
GAM0(26,ii.y_f) = 1;                    GAM0(26,ii.ks_f)  = -pp.Phi*pp.alp;
GAM0(26,ii.l_f) = -pp.Phi*(1-pp.alp);   GAM0(26,ii.eps.a) = -pp.Phi;

% 27. CAPITAL UTILIZATION DEFINITION
GAM0(27,ii.ks_f) = 1;	GAM1(27,ii.k_f) = 1;	GAM0(27,ii.z_f)  = -1;

% 28. OPTIMAL CAPITAL UTILIZATION
GAM0(28,ii.z_f) = 1;                    GAM0(28,ii.rk_f) = -(1-pp.psi)/pp.psi;

% 29. EVOLUTION OF CAPITAL
GAM0(29,ii.k_f) = 1;                    GAM1(29,ii.k_f)   = (1-pp.del)/pp.gam;
GAM0(29,ii.i_f) = -(1-GAM1(29,ii.k_f));	GAM0(29,ii.eps.i) = GAM0(29,ii.i_f)*pp.varphi*pp.gam^2*(1+bg_sig_c);

% 30. PRICE MARKUP
GAM0(30,ii.ks_f) = -pp.alp;             GAM0(30,ii.l_f)   = pp.alp;
GAM0(30,ii.w_f)  = 1;                   GAM0(30,ii.eps.a) = -1;

% 31. RENTAL RATE OF CAPITAL
GAM0(31,[ii.rk_f,ii.ks_f]) = 1;         GAM0(31,[ii.w_f,ii.l_f]) = -1;

% 32. MARGINAL SUBSTITUTION
GAM0(32,ii.w_f) = -1;                   GAM0(32,ii.l_f) = pp.sig_l;
GAM0(32,ii.c_f) = 1/(1-hg);             GAM1(32,ii.c_f) = hg/(1-hg);


% ------------------------------------------------------------------------
% FLEXIBLE PRICE EXPECTATIONS
% ------------------------------------------------------------------------

% 33. CONSUMPTION (c)
GAM0(33,ii.c_f) = 1;    GAM1(33,ii.E.c_f) = 1;      PII(33,ii.E_sh.c_f) = 1;

% 34. REAL VALUE OF CAPITAL STOCK (q)
GAM0(34,ii.q_f) = 1;	GAM1(34,ii.E.q_f) = 1;      PII(34,ii.E_sh.q_f) = 1;

% 35. INVESTMENT (i)
GAM0(35,ii.i_f) = 1;    GAM1(35,ii.E.i_f) = 1;      PII(35,ii.E_sh.i_f) = 1;

% 36. LABOR (l)
GAM0(36,ii.l_f) = 1;    GAM1(36,ii.E.l_f) = 1;      PII(36,ii.E_sh.l_f) = 1;

% 37. RENTAL RATE OF CAPITAL (rk)
GAM0(37,ii.rk_f) = 1;   GAM1(37,ii.E.rk_f) = 1;     PII(37,ii.E_sh.rk_f) = 1;


% ------------------------------------------------------------------------
% EXOGENOUS LAWS-OF-MOTION
% ------------------------------------------------------------------------

% 38. TFP (epsa)
GAM0(38,ii.eps.a) = 1;	GAM1(38,ii.eps.a) = pp.rho_a;	PSI(38,ii.eta_sh.a) = 1;

% 39. NET WORTH (epsb)
GAM0(39,ii.eps.b) = 1;	GAM1(39,ii.eps.b) = pp.rho_b;   PSI(39,ii.eta_sh.b) = 1;

% 40. EXOGENOUS SPENDING (epsg)
GAM0(40,ii.eps.g) = 1;  GAM1(40,ii.eps.g) = pp.rho_g;   PSI(40,ii.eta_sh.g) = 1;
PSI(40,ii.eta_sh.a) = pp.rho_ga;

% 41. INVESTMENT-SPECIFIC (epsi)
GAM0(41,ii.eps.i) = 1;  GAM1(41,ii.eps.i) = pp.rho_i;   PSI(41,ii.eta_sh.i) = 1;

% 42. INTEREST RATE (epsr)
GAM0(42,ii.eps.r) = 1;  GAM1(42,ii.eps.r) = pp.rho_r;	PSI(42,ii.eta_sh.r) = 1;

% 43. PRICE MARKUP (epsp)
GAM0(43,ii.eps.p) = 1;  GAM1(43,ii.eps.p) = pp.rho_p;   PSI(43,ii.eta_sh.p) = 1;
GAM1(43,ii.eta.p) = -pp.mu_p;

% 44. WAGE MARKUP (epsw)
GAM0(44,ii.eps.w) = 1;  GAM1(44,ii.eps.w) = pp.rho_w;   PSI(44,ii.eta_sh.w) = 1;
GAM1(44,ii.eta.w) = -pp.mu_w;


% ------------------------------------------------------------------------
% DEFINE CURRENT INNOVATIONS
% ------------------------------------------------------------------------

% 45. CURRENT PRICE SHOCK
GAM0(45,ii.eta.p) = 1;	PSI(45,ii.eta_sh.p) = 1;

% 46. CURRENT WAGE SHOCK
GAM0(46,ii.eta.w) = 1;  PSI(46,ii.eta_sh.w) = 1;



% ========================================================================
%
%                   SOLVE SYSTEM USING GENSYS
%
%   [G1,C,impact,fmat,fwt,ywt,gev,eu,loose]
%       y(t) = G1*y(t-1)+C+impact*z(t)+ywt*inv(I-fmat*inv(L))*fwt*z(t+1)
%
% ========================================================================

[T1,~,T0,~,~,~,~,~,rc] = gensys(GAM0,GAM1,C,PSI,PII,1+1E-8);
retcode = (max(max(abs(rc)))>1e-10);

% ========================================================================
%
%                           AUGMENT STATES
%
% ========================================================================
% augment lagged variables to be used as states in Kalman filter

[nstates,nshocks] = size(T0);

% ------------------------------------------------------------------------
% INDEX LAGGED VARIABLES
% ------------------------------------------------------------------------
ii.y_1 = 47;
ii.c_1 = 48;
ii.i_1 = 49;
ii.w_1 = 50;


% ------------------------------------------------------------------------
% AUGMENTING T0 AND T1
% ------------------------------------------------------------------------
% recall: y(t) = T1*y(t-1) + T0*z(t)

% EXPAND MATRICES
T1 = [T1,zeros(nstates,4);zeros(4,nstates+4)];
T0 = [T0;zeros(4,nshocks)];

% DEFINITIONS FOR LAGGED VARIABLES
T1(ii.y_1, ii.y) = 1;
T1(ii.c_1, ii.c) = 1;
T1(ii.i_1, ii.i) = 1;
T1(ii.w_1, ii.w) = 1;

% UPDATE NUMBER OF STATES AND SHOCKS
[~,nshocks] = size(T0);


if retcode==0

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%         B. EXPRESS SOLUTION OF MODEL IN STATE SPACE REPRESENTATION
% 
% 1. MEASUREMENT EQUATION:      
%       y(t) = A + H*s(t) + u(t) 
%       u(t) ~ N(0,Su)
% 
% 2. TRANSITION EQUATION:
%       s(t) = F*s(t-1) + B*e(t)
%       e(t) ~ iid N(0,Se)
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ------------------------------------------------------------------------
%                        MEASUREMENT EQUATION
% ------------------------------------------------------------------------

%     nseries = 7;
%     H = zeros(nseries,nstates);
%     A = zeros(nseries,1);
% 
%     % 1. GDP GROWTH
%     A(1) = gam_bar;
%     H(1,y_t) = 1;
%     H(1,y_t1) = -1;
% 
%     % 2. INFLATION
%     A(2) = pii_bar;
%     H(2,pii_t) = 1;
% 
%     % 3. FED FUNDS RATE
%     A(3) = r_bar;
%     H(3,r_t) = 1;
% 
%     % 4. CONSUMPTION GROWTH
%     A(4) = gam_bar;
%     H(4,c_t) = 1;
%     H(4,c_t1) = -1;
% 
%     % 5. INVESTMENT GROWTH
%     A(5) = gam_bar;
%     H(5,i_t) = 1;
%     H(5,i_t1) = -1;
% 
%     % 6. WAGE GROWTH
%     A(6) = gam_bar;
%     H(6,w_t) = 1;
%     H(6,w_t1) = -1;
% 
%     % 7. HOURS
%     A(7) = l_bar;
%     H(7,l_t) = 1;

%     % ASSUME NO MEASUREMENT ERROR
%     Su = zeros(nseries,nseries);


% ------------------------------------------------------------------------
%                        TRANSITION EQUATION
% ------------------------------------------------------------------------

    % VARIANCE OF SHOCKS
    Se = zeros(nshocks,nshocks);
    Se(1,1) = pp.sig_a^2;
    Se(2,2) = pp.sig_b^2;
    Se(3,3) = pp.sig_g^2;
    Se(4,4) = pp.sig_i^2;
    Se(5,5) = pp.sig_r^2;
    Se(6,6) = pp.sig_p^2;
    Se(7,7) = pp.sig_w^2;
    
else
    retcode = 0;
    T1 = NaN;
    T0 = NaN;
    Se = NaN;
end


end