function C = VAR_proxy(data, est_out)

% ========================================================================
% PROXY IDENTIFICATION FROM MERTENS AND RAVN (AER, 2013)
% ========================================================================
%
% -- INPUT --
%   est_out: output from estimation
%   m: instrument
%
% -- OUTPUT --
%   C: matrix with impact of identified shock in first column
%
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% LOAD INPUTS
% ------------------------------------------------------------------------

% VARIANCE
V = est_out.V;
k = 1; N = size(V,1);

% DATA
Y = data.Y; n_lag = size(est_out.B,2)/size(est_out.B,1);
X = lagmatrix(Y,1:n_lag);
Y = Y(est_out.sample_no_NaN,:);
X = X(est_out.sample_no_NaN,:);
e = Y - est_out.cons' - X*est_out.B';

% INSTRUMENT
m = data.m(est_out.sample_no_NaN);


% ------------------------------------------------------------------------
% PROXY IDENTIFICATION
% ------------------------------------------------------------------------

% REGRESSION TO OBTAIN CORRELATION BETWEEN PROXY AND SHOCK
PHI_C = [ones(length(m),1) m]\e;
PHI_C = PHI_C(2:end,:);
PHI_C_11 = PHI_C(1:k,1:k); PHI_C_21 = PHI_C(1:k,k+1:N);
C21i_C11 = (PHI_C_11\PHI_C_21)';

% BLOCKS OF SHOCK COVARIANCE MATRIX
V11 = V(1:k,1:k); V21 = V(k+1:N,1:k); V22 = V(k+1:N,k+1:N);

% SOLVE FOR IMPULSE RESPONSE ON IMPACT
ZZp = C21i_C11*V11*C21i_C11' - (V21*C21i_C11'+C21i_C11*V21') + V22;
C12_C12p = (V21-C21i_C11*V11)'*(ZZp\(V21- C21i_C11*V11));
C11_C11p = V11-C12_C12p;
C11 = sqrt(C11_C11p);
C1 = [C11; C21i_C11*C11];

% ORDER IDENTIFIED SHOCK FIRST
C = [C1, zeros(N, N-1)];

end