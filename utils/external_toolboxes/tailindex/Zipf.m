function [QQ SS1 SS2 BF] = Zipf(X,k);
% Kratz and Resnick (1996) estimator of tail index
% X is a vector of observations in descending order, the k largest ones are
% used to estimate gamma, one over the tail-index
% Also returns the first two Danielsson, Jansen, De Vries (1996) estimators
%
% First version: Richard Tol, 2 November 2019
% This version: Richard Tol, 2 November 2019

lnX = log(X);
n = length(X);
for j=1:k,
    lnj(j) = log((k+1)/j);
    lnn(j) = log(n/j);
    lnn2(j) = lnn(j)*lnn(j);
end
SS1 = lnn*lnX(1:k) / sum(lnn2);

W = [ones(k,1) lnj'];
q = inv(W'*W)*W'*lnX(1:k);
QQ = q(2);
W = [ones(k,1) lnX(1:k)];
s = inv(W'*W)*W'*lnj';
SS2 = s(2);

BF = sqrt(QQ*SS2);