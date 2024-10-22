function A = var_ols_ridge(X, lambda)
% VAR_OLS_RIDGE computes coefficients using OLS with ridge regularization
%   A = VAR_OLS_RIDGE(X, lambda) computes coefficient matrix A for a VAR
%   model using ordinary least squares (OLS) with ridge regularization.
%
%   Inputs:
%   X       : Input data matrix of size (n x m x p), where n is the number
%             of variables, m is the number of time points, and p is the number
%             of trials.
%   lambda  : Regularization parameter for ridge regression.
%
%   Output:
%   A       : Coefficient matrix of size (n x n x p), where A(:,:,k) is the
%             coefficient matrix for lag k.

[n, m, p] = size(X);  % Dimensions of input matrix X
p1 = 1 + p;
pn = p * n;

M = m - p;

% Stack lags
X0 = reshape(X(:, p1:m, :), n, M * p);  % Concatenate unlagged observations
XL = zeros(n, p, M);
for k = 1:p
    XL(:, k, :) = reshape(X(:, p1 - k:m - k, :), n, M);  % Concatenate k-lagged observations
end
XL = reshape(XL, pn, M);  % Stack lags

XXL = XL * XL';
A = (X0 * XL') / (XXL + lambda * eye(n));  % Compute coefficients using ridge regression

A = reshape(A, n, n, p);  % Reshape A to n x n x p
A = A - diag(diag(A));  % Remove diagonal elements
end
