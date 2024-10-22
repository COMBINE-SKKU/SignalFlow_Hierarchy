%% Generalize to Test Set
% This script estimates correlation values and ratios for different
% algorithms using a test set, and saves the results to .mat files.

clear; close all; clc;

% Load data
rDCMs = importdata(fullfile('fig2(mac)data', 'rDCMs(macaque).mat'));
VARs = importdata(fullfile('fig2(mac)data', 'VARs(macaque).mat'));
FASKs = importdata(fullfile('fig2(mac)data', 'FASKs(macaque).mat'));
FLN = importdata(fullfile('fig2(mac)data', 'FLN40.mat'));
SLN = importdata(fullfile('fig2(mac)data', 'SLN40.mat'));
EC_betas_FLN = importdata(fullfile('fig2(mac)results', 'beta_values_macaque.mat'));

% Allocate the size of the result matrices by adding threshold dimension
EC_correlation_values_FLN = zeros(10, 4);
EC_correlation_values_SLN = zeros(10, 4);

% Initialize arrays to store p-values for correlations
EC_p_values_FLN = zeros(10, 4);
EC_p_values_SLN = zeros(10, 4);

% Initialize arrays to store the ratios
ratio_array = zeros(10, 4);

% Define indices for feedback and feedforward
FBindex = SLN < 0.5;
FFindex = SLN >= 0.5;

% Calculate mean betas
Betas = mean(EC_betas_FLN);

% Initialize array to store iEC values
iEC_all = [];
iter = 1;

% Iterate over the test set
for sub = 9:18
    % Normalize matrices and combine using Betas
    rDCM = rDCMs(:,:,sub);
    VAR = VARs(:,:,sub);
    FASK = FASKs(:,:,sub);
    
    rDCM = rDCM/max(rDCM(:)) * 0.2;
    VAR = VAR/max(VAR(:)) * 0.2;
    FASK = FASK/max(FASK(:)) * 0.2;
    
    iEC = (Betas(1) * rDCM + Betas(2) * VAR + Betas(3) * FASK) / 3;
%     iEC_all = cat(3, iEC_all, iEC);
    
    % Calculate correlation and p-values for FLN
    [FLN_rDCM_corr, FLN_rDCM_p] = corr(FLN(:), rDCM(:));
    [FLN_VAR_corr, FLN_VAR_p] = corr(FLN(:), VAR(:));
    [FLN_FASK_corr, FLN_FASK_p] = corr(FLN(:), FASK(:));
    [FLN_iEC_corr, FLN_iEC_p] = corr(FLN(:), iEC(:));
    
    % Calculate correlation and p-values for SLN
    [SLN_rDCM_corr, SLN_rDCM_p] = corr(SLN(:), rDCM(:));
    [SLN_VAR_corr, SLN_VAR_p] = corr(SLN(:), VAR(:));
    [SLN_FASK_corr, SLN_FASK_p] = corr(SLN(:), FASK(:));
    [SLN_iEC_corr, SLN_iEC_p] = corr(SLN(:), iEC(:));
    
    % Store correlation values
    EC_correlation_values_FLN(iter, :) = [FLN_rDCM_corr, FLN_VAR_corr, FLN_FASK_corr, FLN_iEC_corr];
    EC_correlation_values_SLN(iter, :) = [SLN_rDCM_corr, SLN_VAR_corr, SLN_FASK_corr, SLN_iEC_corr];
    
    % Store p-values
    EC_p_values_FLN(iter, :) = [FLN_rDCM_p, FLN_VAR_p, FLN_FASK_p, FLN_iEC_p];
    EC_p_values_SLN(iter, :) = [SLN_rDCM_p, SLN_VAR_p, SLN_FASK_p, SLN_iEC_p];
    
    % Calculate and store the ratios
    iEC_FF = iEC(FFindex);
    iEC_FB = iEC(FBindex);
    
    ratio_neg_FF = sum(iEC_FF < 0) / length(iEC_FF);
    ratio_pos_FF = sum(iEC_FF > 0) / length(iEC_FF);
    ratio_neg_FB = sum(iEC_FB < 0) / length(iEC_FB);
    ratio_pos_FB = sum(iEC_FB > 0) / length(iEC_FB);
    
    ratio_array(iter, :) = [ratio_neg_FF, ratio_pos_FF, ratio_neg_FB, ratio_pos_FB];
    
    iter = iter + 1;
end

% Save the estimated correlation values
save('fig2(mac)results/EC_correlation_FLN.mat', 'EC_correlation_values_FLN');
save('fig2(mac)results/EC_FF_FB.mat', 'ratio_array');
%% save edge strength correlation data for box plot

% Extract correlation values for the selected atlas
correlation_values = EC_correlation_values_FLN;
num_iter = size(correlation_values,1);

% Reshape data vectors for each measure
data_correlation = reshape(correlation_values, [], 1);

% Create groups for each measure
group_correlation = repmat({'FLN'}, size(data_correlation));

% Create legends for each measure
legend = repmat({'rDCM', 'VAR', 'FASK', 'iEC'}, num_iter, 1);
legend_all = reshape(legend, [], 1);

% Save as table
result = table(data_all, legend_all, group_all, 'VariableNames', {'Data', 'Legend', 'Group'});

save('fig2(mac)results/Figure2e_edge_strength_correlation','result')

%% save feedforward/back data for box plot

% Prepare data for table
temp1 = reshape(ratio_array, [], 1);
temp1origin = repmat({'neg', 'pos', 'neg', 'pos'}, 10, 1);
temp1origin = reshape(temp1origin, [], 1);
color = repmat({'FF', 'FF', 'FB', 'FB'}, 10, 1);
color = reshape(color, [], 1);

% Store results in table and save
result = [array2table(temp1), cell2table([color, temp1origin])];
save('fig2(mac)results/Figure2f_EC_strength_for_FF_FB.mat', 'result');

