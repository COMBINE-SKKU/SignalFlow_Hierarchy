function gc = run_gc(timeseries, model_order_select, lambda)
% RUN_GC Performs Granger Causality analysis on multivariate time series data
%
% This function estimates Granger causality between variables in a multivariate
% time series using a state space model approach.
%
% Inputs:
%   timeseries        - Matrix of time series data (time x variables)
%   model_order_select- Boolean flag to determine if model order should be
%                       automatically selected (true) or set to 1 (false)
%   lambda            - Regularization parameter for ridge regression
% Outputs:
%   gc               - Matrix of pairwise Granger causality values
%
% The function uses either AIC or BIC criteria for model order selection if
% model_order_select is true, otherwise uses order 1. Estimation is performed
% using Ordinary Least Squares (OLS) regression.

% Default parameters
reg_mode = 'OLS';    % VAR model estimation regression mode
ic_mode = 'OLS';     % Information criteria regression mode
max_order = 10;      % Maximum model order for estimation
order_criterion = 'AIC'; % Criterion for model order selection

% Model order selection
if model_order_select
    [~, ~, aic_order, ~] = tsdata_to_infocrit(timeseries, ...
        max_order, ic_mode);
    
    % Select order based on AIC criterion
    model_order = aic_order;
    fprintf('\nUsing AIC best model order = %d\n', model_order);
else
    model_order = 1;
    fprintf('\nUsing fixed model order = 1\n');
end

% Estimate VAR model
[var_coeffs, residuals_cov] = tsdata_to_var(timeseries, model_order, ...
    reg_mode, lambda);

% Calculate pairwise Granger causality

gc = var_to_pwcgc(var_coeffs, residuals_cov, timeseries);

% Replace NaN values with zeros
gc(isnan(gc)) = 0;


end

