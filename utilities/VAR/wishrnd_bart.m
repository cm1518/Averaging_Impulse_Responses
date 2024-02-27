function A = wishrnd_bart(S, d, n)

% ========================================================================
% DRAW FROM WISHART DISTRIBUTION USING BARTLETT DECOMPOSITION
% >> similar to rwwish function from Sims VAR package
% >> faster than Matlab functions for large number of draws
% ========================================================================
%
% -- INPUT --
%   S: shape
%   d: degrees of freedom
%   n: number of draws
%
% -- OUTPUT --
%   A(d-by-d-by-n_draw): array with n draws from Wishart(S,df)
%
% ------------------------------------------------------------------------

% NUMBER OF DIMENSIONS
p = size(S,1);

% SET UP MATRIX FOR WISHART DRAWS
A = zeros(p*p,n);

% i'TH DIAGONAL ELEMENT: CHI-SQUARE WITH df+1-i DEGREES OF FREEDOM
A(p*(0:p-1) + (1:p),:) = sqrt(chi2rnd(repmat(d+1-(1:p), n, 1))');

% LOWER TRIANGLE ELEMENTS: STANDARD NORMAL
i_low = NaN(p*(p-1)/2, 1);
ii = 0;
for i = 1:p-1
    i_low(ii + (1:p-i)) = p*(i-1) + (i+1:p);
    ii = ii + p - i;
end
A(i_low,:) = randn(p*(p-1)/2, n);

% MULTIPLY BY CHOLESKY TO GET CORRECT SCALE
C = chol(S,'lower');
A = C*reshape(A,p,p*n);

% RESHAPE AND MULTIPLY BY TRANSPOSE
A = reshape(A,p,p,n);
for i = 1:n
    A(:,:,i) = A(:,:,i)*A(:,:,i)';
end

end