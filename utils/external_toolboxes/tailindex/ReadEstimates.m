%% ReadEstimates
%Reads estimates of social cost of carbon and their weights
%
% First version: Richard Tol, 8 November 2011
% This version: Richard Tol, 25 February 2020

display('Read data');

vFileToOpen = 'C:\Users\rtol\Google Drive\MargCostMeta\socialcostcarbon.xlsx';
SCC = xlsread(vFileToOpen, 'Data', 'L4:L2789');
AuthorWeight = xlsread(vFileToOpen, 'Data', 'K4:K2789');
QualWeight= xlsread(vFileToOpen, 'Data', 'P4:P2789');
Censor= xlsread(vFileToOpen, 'Data', 'I4:I2789');

[SCCs vindx] = sort(SCC,'descend');
AuthorWeight = AuthorWeight(vindx);
QualWeight = QualWeight(vindx);
Censor = Censor(vindx);
TotWeight = AuthorWeight.*QualWeight;
Censored = Censor.*TotWeight;

clear v*