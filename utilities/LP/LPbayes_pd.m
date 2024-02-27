function pd_out = LPbayes_pd(data, draw, test, ind, H, identification)

% ------------------------------------------------------------------------
% PRELIMINARIES
% ------------------------------------------------------------------------

% LOAD DATA
if strcmp(identification,'internal IV')
    Y = [data.m, data.Y];  %ind = ind + 1;
else
    Y = data.Y;
end

% DIMENSIONS
[T, N] = size(Y);           % number of periods, variables
n_draw = size(draw{1}.B,3);	% number of draws
n_lag = size(draw{1}.B,2)/N;   % number of lags
T_test = length(test);      % number of observations in test set

% LAGGED Y
X = lagmatrix(Y, 1:n_lag);

% MATRICES FOR STORAGE
pd      = NaN(H+1, T_test, length(ind), n_draw);
EY_draw = NaN(H+1, T_test, length(ind), n_draw);
VY_draw = NaN(H+1, T_test, length(ind), n_draw);
% pd = NaN(H+1, T_test, n_draw);
% EY_draw = NaN(H+1, T_test, n_draw);
% VY_draw = NaN(H+1, T_test, n_draw);


% ------------------------------------------------------------------------
% PREDICTIVE DENSITIES
% ------------------------------------------------------------------------

for ii = 1:n_draw
    
    % EXTRACT PARAMETER DRAWS
    cons1 = draw{1}.cons(:,ii);
    B1 = draw{1}.B(:,:,ii);
    V1 = draw{1}.V(:,:,ii);
    
    % COMPUTE PREDICTIVE DENSITIES
    for tt = 1:T_test        
        for hh = 1:H+1
            if test(tt)+hh-1 <= T
                if hh == 1 
                    % HORIZON 0: CONDITION ON LAGGED DATA AND SHOCK
                    EY_prev = X(test(tt),:)';   % previous observations
                    EY = NaN(N*n_lag,1);        % vector for conditional expectation
                    VY = zeros(N);              % matrix for conditional variance (not companion form)
                    EY(1) = Y(test(tt),1);
                    
                    EY(2:N) = cons1(2:N) + B1(2:N,:)*EY_prev + ...
                        V1(2:N,1)/V1(1,1)*(Y(test(tt),1) - cons1(1) - B1(1,:)*EY_prev);
                    VY(2:N,2:N) = V1(2:N,2:N) - V1(2:N,1)/V1(1,1)*V1(1,2:N);
                    EY(N+1:end) = EY_prev(1:N*(n_lag-1));
                    EY_prev = EY;
                    VY1 = VY;
                else
                    % HORIZON >0: ITERATE FORWARD
                    % mean
                    EY(1:N) = draw{hh}.cons(:,ii) + draw{hh}.B(:,:,ii)*EY_prev;
                    % variance
                    VY = VY1 + draw{hh}.V(:,:,ii);
                end
                for jj = 1:length(ind)
                    EY_draw(hh,tt,jj,ii) = EY(ind(jj));               % expectation of interest
                    VY_draw(hh,tt,jj,ii) = VY(ind(jj),ind(jj));       % variance of interest
                end
            end
        end
    end

end

% PREDICTIVE DENSITY
for hh = 1:H+1
    tmp = test+hh-1 <= T;
    test_sub = test(tmp);
    for jj = 1:length(ind)
        pd(hh,tmp,jj,:) = normpdf(Y(test_sub+hh-1,ind(jj)),...
            squeeze(EY_draw(hh,tmp,jj,:)), squeeze(sqrt(VY_draw(hh,tmp,jj,:))));
    end
end
pd = mean(pd,4);
pd_out.pd      = permute(     pd,[2,3,1]  );
% for hh = 1:H+1
%     tmp = test+hh-1 <= T;
%     test_sub = test(tmp);
%     pd(hh,tmp,:) = normpdf(Y(test_sub+hh-1,ind),...
%         squeeze(EY_draw(hh,tmp,:)), squeeze(sqrt(VY_draw(hh,tmp,:))));
% end
% pd = mean(pd,3)';
% pd_out.pd = pd; %pd_out.EY_draw=EY_draw; pd_out.VY_draw=VY_draw;

end