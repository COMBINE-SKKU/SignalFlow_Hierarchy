function bestLambda = VAR_optimization2(data, atlas_name, state)
% VAR_OPTIMIZATION2 calculates optimal lambda value for VAR model
%   bestLambda = VAR_OPTIMIZATION2(data, atlas_name, state) calculates the
%   optimal regularization parameter lambda for a Vector Autoregressive
%   (VAR) model using multiple timeseries data.
%
%   Inputs:
%   data        : Cell array of timeseries data from multiple subjects
%                 Each cell should contain a matrix of size (n x m),
%                 where n is the number of variables and m is the number
%                 of time points.
%   atlas_name  : String representing the brain atlas used ('Schaefer',
%                 'MMP', etc.).
%   state       : String indicating the state of fMRI data ('resting',
%                 'movie', 'pain', etc.).
%
%   Output:
%   bestLambda  : Optimal regularization parameter lambda that maximizes
%                 functional connectivity (FC) recovery across subjects.

    % Search range of lambda values
    lambda_values = 0:200:1e+3;

    % Initialize an array to store correlation coefficient (CC) values
    CC_values = zeros(size(lambda_values));

    % Compute functional connectivity (FC) matrix
    FC = corr(data);
    FC = FC - diag(diag(FC));  % Remove diagonal elements

    % Iterate over lambda values
    for i = 1:length(lambda_values)
        lambda = lambda_values(i);

        % Calculate VAR with current lambda
        VAR = var_ols_ridge(data', lambda);

        % Normalize VAR matrix and scale for optimal G calculation
        VAR = VAR / max(VAR(:)) * 0.1;

        % Calculate optimal G matrix
        VAR_G = optimal_G(VAR, FC, atlas_name, state);

        % Calculate correlation coefficient (CC) for current VAR
        CC_values(i) = calculateCC(VAR, FC, VAR_G, atlas_name, state);
    end

    % Find lambda value corresponding to maximum CC
    [~, maxIndex] = max(CC_values);
    bestLambda = lambda_values(maxIndex);
end
