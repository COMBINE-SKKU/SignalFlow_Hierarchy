function [bestLambda] = VAR_optimization(BOLD, Target)
% VAR_OPTIMIZATION performs ridge regression optimization for VAR model
%   [bestLambda] = VAR_OPTIMIZATION(BOLD, Target) performs ridge regression
%   optimization for a VAR (Vector Autoregressive) model to find the best
%   regularization parameter lambda.
%
%   Inputs:
%   BOLD        : Input BOLD data matrix of size (n x m x p), where n is the
%                 number of variables, m is the number of time points, and p
%                 is the number of trials.
%   Target      : Target matrix for comparison, typically a vector (n*m x 1).
%
%   Output:
%   bestLambda  : Best regularization parameter lambda that minimizes the
%                 difference between predicted and actual target.
%

    lambda_values = 0:500:1e+5;  % The lambda search space
    nLambda = length(lambda_values);  % Number of lambda values to test
    performanceValues = Inf * ones(1, nLambda);  % Initialize array to store performance for each lambda

    flatTarget = Target(:);  % Flatten target matrix into a vector

    % Parallel loop to compute performance for each lambda
    parfor i = 1:nLambda
        lambda = lambda_values(i);  % Current lambda value
        VAR = var_ols_ridge(BOLD, lambda);  % Compute VAR model using ridge regression
        predictedTarget = VAR(:);  % Flatten VAR matrix into a vector
        CC = corr(flatTarget, predictedTarget);  % Compute Pearson's correlation coefficient
        performanceValues(i) = 1 - CC;  % Compute performance measure
    end

    % Find the best performance and corresponding lambda
    [~, bestIdx] = min(performanceValues);  % Find index of minimum performance
    bestLambda = lambda_values(bestIdx);  % Get corresponding lambda value
end
