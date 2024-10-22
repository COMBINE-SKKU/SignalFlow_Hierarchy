%% Estimate Beta Values Using Training Set
% This script estimates beta values for rDCM, VAR, and FASK algorithms using
% a subset of the training set and saves the results to a .mat file.
% Note: This example uses 20 sample subjects from the actual dataset.

clear; close all; clc;

% Load data
FCs = importdata(fullfile('fig2(FCrec)data', 'MMP_resting_FCs.mat'));
rDCMs = importdata(fullfile('fig2(FCrec)data', 'MMP_resting_rDCMs.mat'));
VARs = importdata(fullfile('fig2(FCrec)data', 'MMP_resting_VARs.mat'));
FASKs = importdata(fullfile('fig2(FCrec)data', 'MMP_resting_FASKs.mat'));

% Set variables
nSub = 10; % Number of subjects used for training
nAlgos = 3; % Number of algorithms
EC_betas = zeros(nSub, nAlgos);
atlas_name = 'MMP';
state = 'resting';

% Train beta values
for sub = 1:nSub
    % Extract and normalize connectivity matrices
    rDCM = rDCMs(:,:,sub) / max(rDCMs(:,:,sub), [], 'all') * 0.1;
    VAR = VARs(:,:,sub) / max(VARs(:,:,sub), [], 'all') * 0.1;
    FASK = FASKs(:,:,sub) / max(FASKs(:,:,sub), [], 'all') * 0.1;
    FC = FCs(:,:,sub);

    % Prepare input matrices
    input_matrices = {rDCM, VAR, FASK};

    % Estimate Betas
    [Betas, ~] = findBetas(input_matrices, 'Target', FC, ...
        'Atlas', atlas_name, 'State', state);
    
    % Store beta values
    EC_betas(sub, :) = Betas;
end

%% Generalizability of Beta Values
% This section tests the generalizability of the estimated beta values on a test set.

% Load actual beta values from the original training set
EC_betas = importdata(fullfile('fig2(FCrec)data', 'MMP_iEC_betas_training.mat'));

% Calculate median Betas across the training set
Betas = median(EC_betas, 1);

% Initialize variables
nTestSub = 10; % Number of subjects in the test set
EC_correlation_values_test = zeros(nTestSub, 4); % To store correlation values

% Test beta values on new subjects
for sub = 11:20
    % Extract and normalize connectivity matrices
    rDCM = rDCMs(:, :, sub) / max(rDCMs(:, :, sub), [], 'all') * 0.1;
    VAR = VARs(:, :, sub) / max(VARs(:, :, sub), [], 'all') * 0.1;
    FASK = FASKs(:, :, sub) / max(FASKs(:, :, sub), [], 'all') * 0.1;
    FC = FCs(:, :, sub);

    % Compute integrated effective connectivity (iEC)
    iEC = (Betas(1) * rDCM + Betas(2) * VAR + Betas(3) * FASK) / 3;

    % Calculate correlations
    rDCM_corr = calculateCC(rDCM, FC, 1.5, atlas_name, state);
    VAR_corr = calculateCC(VAR, FC, 3.5, atlas_name, state);
    FASK_corr = calculateCC(FASK, FC, 10, atlas_name, state);
    iEC_corr = calculateCC(iEC, FC, 1, atlas_name, state);

    % Store the correlation values
    EC_correlation_values_test(sub-10, :) = [rDCM_corr, VAR_corr, FASK_corr, iEC_corr];
end

%% Simple Visualization of Figure 2h
% This script generates boxplots for correlation values from Figure 2h.
% One should use BoxPlotPro from MATLAB FileExchange to reproduce figures
% from the manuscript with the same quality.

data = importdata(fullfile('fig2(FCrec)data', 'MMP_iEC_across_subjects.mat'));

% Generate legends for each measure
num_iter = size(data, 1);
legend = repmat({'rDCM', 'VAR', 'FASK', 'iEC'}, num_iter, 1);
legend = reshape(legend, [], 1);

% Create the boxplot and set box colors
figure;
h = boxplot(data, legend, 'Colors', 'k', 'Symbol', ''); % 'k' sets the box color to black

% Define colors for each box
col_vals = [...
    'fd'; 'e3'; 'a2';
    'f9'; 'dc'; 'e1';
    '9b'; 'c5'; 'db';
    'f6'; '69'; '45'];
col_vals = reshape(hex2dec(col_vals), [3 numel(col_vals)/6])' ./ 255;
P = size(col_vals, 1);
colors = interp1(1:size(col_vals, 1), col_vals, linspace(1, P, 4), 'linear');

% Set the colors for each box
hGroups = findobj(gca, 'Tag', 'Box');
for j = 1:length(hGroups)
    patch(get(hGroups(j), 'XData'), get(hGroups(j), 'YData'), colors(j, :), 'FaceAlpha', 0.5); % Fill the box with color
end

% Hold the current plot
hold on;

% Add scatter plot points in grey color
unique_groups = {'rDCM', 'VAR', 'FASK', 'iEC'};
num_groups = length(unique_groups);

% Adjust jitter for better visibility
jitterAmount = 0.1;

for i = 1:num_groups
    % Find the x-position for each group
    x_pos = repmat(i, sum(strcmp(legend, unique_groups{i})), 1) + (rand(sum(strcmp(legend, unique_groups{i})), 1) - 0.5) * jitterAmount;
    y_pos = data(strcmp(legend, unique_groups{i}));
    scatter(x_pos, y_pos, 5, 'filled', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', [0.5, 0.5, 0.5]); % Grey color
end

ylim([0.3 1])
