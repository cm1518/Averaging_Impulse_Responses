function pd_out = LPbayes_pd_uncond(data, draw, test, ind, H, identification)


% ------------------------------------------------------------------------
% PRELIMINARIES
% ------------------------------------------------------------------------

% LOAD DATA
Y = data.Y;
if strcmp(identification,'shock-augmented')
    Y = [data.m, data.Y];  
    ind = ind + 1; 
elseif strcmp(identification,'proxy')
    %Y = data.Y;
    Y = [data.m, data.Y];  
    ind = ind + 1; 
    % We don't know how to do a proxy Bayesian VAR, so we treat it like
    % shock-augmented
end

% DIMENSIONS
[T, N] = size(Y);           % number of periods, variables
n_draw = size(draw{1}.B,3);	% number of draws
n_lag = size(draw{1}.B,2)/N;   % number of lags
T_test = length(test);      % number of observations in test set

% LAGGED Y
X = lagmatrix(Y, 1:n_lag);

% MATRICES FOR STORAGE
pd = NaN(H+1, T_test, n_draw);
EY_draw = NaN(H+1, T_test, n_draw);
VY_draw = NaN(H+1, T_test, n_draw);


% ------------------------------------------------------------------------
% PREDICTIVE DENSITIES
% ------------------------------------------------------------------------

for i = 1:n_draw
    
%     % EXTRACT PARAMETER DRAWS
%     cons1 = draw{1}.cons(:,i);
%     B1 = draw{1}.B(:,:,i);
%     V1 = draw{1}.V(:,:,i);
    
    % COMPUTE PREDICTIVE DENSITIES
    for t = 1:T_test        
        for h = 1:H+1
            if test(t)+h-1 <= T
            if h == 1 
                % HORIZON 0: Define first horizon EY and VY
                EY1 = X(test(t),:)';        % previous observations
                EY = NaN(N*n_lag,1);        % vector for conditional expectation
                VY = zeros(N);              % matrix for conditional variance (not companion form)
                
                % mean
                EY(1:N) = draw{h}.cons(:,i) + draw{h}.B(:,:,i)*EY1;
                EY(N+1:end) = EY1(1:N*(n_lag-1));
                
                % variance
                VY = draw{h}.V(:,:,i);
                
                EY1 = EY;
                VY1 = VY;
            else
                % HORIZON >0: ITERATE FORWARD
                % mean
                EY(1:N) = draw{h-1}.cons(:,i) + draw{h-1}.B(:,:,i)*EY1;
                % variance
                VY = VY1 + draw{h-1}.V(:,:,i);
            end
            EY_draw(h,t,i) = EY(ind);               % expectation of interest
            VY_draw(h,t,i) = VY(ind,ind);           % variance of interest
            end
        end
    end

end

% PREDICTIVE DENSITY
for h = 1:H+1
    tmp = test+h-1 <= T;
    test_sub = test(tmp);
    pd(h,tmp,:) = normpdf(Y(test_sub+h-1,ind),...
        squeeze(EY_draw(h,tmp,:)), squeeze(sqrt(VY_draw(h,tmp,:))));
end
pd = mean(pd,3)';
pd_out.pd = pd; pd_out.EY_draw=EY_draw; pd_out.VY_draw=VY_draw;

end