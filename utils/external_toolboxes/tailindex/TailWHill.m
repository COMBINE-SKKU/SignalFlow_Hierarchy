function [Hill Hillsd DJV1 DJV2 AM AMsd T1 T1sd T2 T3 D Dsd] = TailWHill(X,W,k);
% H = Hill (Annals Stats, 1975) weighted estimator of tail index
% DJV = Danielsson, Jansen, De Vries (Comm Stats, 1996) weighted estimators
% AM = Aban, Meerschaert (J Stat Plan Inf, 2004) weighted estimator
% T = Tripathi (Stat Meth, 2014) weighted estimator
% D = Dekkers, Einmahl, de Haan (Ann Stat, 1989) weighted estimator
%
% X is a vector of observations in descending order, the k largest ones are
% used to estimate the tail-index
%
% W is the vector of weights. If W=1, this returns the original
% estimators
%
% First version: Richard Tol, 8 November 2011
% This version: Richard Tol, 29 February 2020

lnX = log(X);
kw = sum(W(1:k));
W = W*k/kw;
M1 = sum(W(1:k).*(lnX(1:k) - lnX(k+1)))/k;
M2 = sum(W(1:k).*(lnX(1:k) - lnX(k+1)).^2)/k;
M3 = sum(W(1:k).*(lnX(1:k) - lnX(k+1)).^3)/k;

Hill = 1/M1;
Hillsd = k*Hill/(k-1)/sqrt(k-1);
DJV1 = 2*M1/M2;
DJV2 = 3*M2/M3; 
AM = Hill*(k-1)/k;
AMsd = Hill/sqrt(k-1);
T1 = (k-3)/M1/k;
T1sd = Hillsd*(k-3)/k;
T2 = (k-2)/(M1*k-max(0,lnX(k)));
T3 = (k-3)/(M1*k-max(0,lnX(k)));
D = M1 + 1 -0.5/(1-M1*M1/M2);
D = 1/D;
Dsd = D/sqrt(k);