function test_data = create_test_data(spec, est, Y_dropped, z_dropped, ind, h)

n_lag = est.n_lag;

% TESTING SAMPLE (AG ADDED)
Y_test = lagmatrix(Y_dropped(:,ind), -h+1);
X_test = [ z_dropped, lagmatrix([z_dropped, Y_dropped], 1:n_lag)];
  
%Added by CGL
if strcmp(est.id,'non-recursive')
    X_test = [ z_dropped, lagmatrix([z_dropped, Y_dropped], 1:n_lag)];
elseif  strcmp(est.id,'recursive')
    endog_var_i = 1:length(est.var);
    endog_var_i(ind) = []; %Don't include current LHS var
    X_test = [ z_dropped, Y_dropped(:, endog_var_i), lagmatrix([z_dropped, Y_dropped], 1:n_lag)];
end

% INCLUDE CONSTANT IF NEEDED
if strcmp(spec.constant,'true')
    X_test = [X_test, ones(size(X_test,1),1)];
end

test_data.Y = Y_test;
test_data.X = X_test;

end