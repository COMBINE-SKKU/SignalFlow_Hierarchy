function [H Hsd] = TailHuisman(Hill,min,max,N,Hm);
% Huisman, Koedijk, Kool, Palm (J Bus Econ Stat, 2001) estimator of tail index
%
% Hill is a vector of estimates (e.g., Hill's) of the tail-index
% min is the smallest sample size used to estimate the tail-index
% max is the largest sample size used to estimate the tail-index
% Hm is the estimate of the tail-index for a tail sample size of max+1
% N is the sample size
%
% First version: Richard Tol, 8 March 2020
% This version: Richard Tol, 10 March 2020

W = zeros(max-min+1,max-min+1);
for i=1:max-min+1,
    trend(i) = i+min-1;
    W(i,i) = sqrt(trend(i));
end

intercept = ones(1,length(trend));
X = [intercept' trend'];

A = zeros(max-min+1,max-min+2);
for i=1:max-min+1,
    A(i,max-min+2-i) = -1;
    for j=(max-min+3-i):(max-min+2)
        A(i,j) = 1/i;
    end
end

alpha = [Hill; Hm];
for i=1:max-min+2,
    mu(i) = log((1-(i+min-1)/N)^(-1/alpha(i)));
end

for i=1:max-min+2,
    for j=1:max-min+2,
        if i <= j,
            sigma(i,j)= (i+min-1)*(1-(j+min-1)/N)*mu(i)^(alpha(i)+1)*mu(j)^(alpha(j)+1)/alpha(i)/alpha(j);
        else
            sigma(i,j)=0;
        end
    end
end
omega = A*sigma*A';

beta  = inv(X'*W'*W*X)*X'*W'*W*Hill;
H = beta(1);
resid = Hill - X*beta;
SSR = resid'*W'*W*resid;
sigsq = SSR/(length(resid)-2);

cov = sigsq*inv(X'*W'*W*X)*X'*W'*W*omega*W'*W*X*inv(X'*W'*W*X);
Hsd = sqrt(cov(1,1));