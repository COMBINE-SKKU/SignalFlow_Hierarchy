%% Estimate Beta Values Using Training Set
% This script estimates beta values for rDCM, VAR, and FASK algorithms using
% the training set and saves the results to a .mat file.

clear; close all; clc;

% Load data
rDCMs = importdata(fullfile('fig2(mac)data', 'rDCMs(macaque).mat'));
VARs = importdata(fullfile('fig2(mac)data', 'VARs(macaque).mat'));
FASKs = importdata(fullfile('fig2(mac)data', 'FASKs(macaque).mat'));
FLN = importdata(fullfile('fig2(mac)data', 'FLN40.mat'));
SLN = importdata(fullfile('fig2(mac)data', 'SLN40.mat'));

% Apply log-linear transformation to preserve the connectivity information
% while compressing the connection values to a reasonable range.
% 1.2 and 0.3 are values from Mejias et al., 2016
FLN = 1.2 * FLN.^0.3;

% Initialize arrays to store Betas and beta histories
EC_betas_FLN = zeros(8, 3);
betaHistories_FLN = [];

% Training sample
num_iter = 8;
for sub = 1:num_iter
    % Extract matrices for the current iteration
    rDCM = rDCMs(:, :, sub);
    VAR = VARs(:, :, sub);
    FASK = FASKs(:, :, sub);

    % Normalize matrices and store in cell array
    input_matrices = cell(3, 1);
    input_matrices{1} = rDCM / max(rDCM(:)) * 0.2;
    input_matrices{2} = VAR / max(VAR(:)) * 0.2;
    input_matrices{3} = FASK / max(FASK(:)) * 0.2;

    % Estimate Betas using the findBetas function
    [Betas, ~] = findBetas(input_matrices, 'Target', FLN);
    EC_betas_FLN(sub, :) = Betas;
end

% Save the estimated beta values to a .mat file
save(fullfile('fig2(mac)results', 'beta_values_macaque.mat'), 'EC_betas_FLN');
