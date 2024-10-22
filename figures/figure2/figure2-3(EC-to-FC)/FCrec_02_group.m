%% Estimate Beta Values Using Training Set (group level)
% This script estimates beta values for rDCM, VAR, and FASK algorithms using
% a subset of the training set and saves the results to a .mat file.

clear; close all; clc;

atlas_name = 'MMP';
state = 'resting';

% Load data
FC_training = importdata(fullfile('fig2(FCrec)data', 'MMP_FC_group1.mat'));
rDCM_training = importdata(fullfile('fig2(FCrec)data', 'MMP_rDCM_group1.mat'));
VAR_training = importdata(fullfile('fig2(FCrec)data', 'MMP_VAR_group1.mat'));
FASK_training = importdata(fullfile('fig2(FCrec)data', 'MMP_FASK_group1.mat'));

input_matrices = cell(3,1);
input_matrices{1} = rDCM_training/max(rDCM_training(:))*0.1;
input_matrices{2} = VAR_training/max(VAR_training(:))*0.1;
input_matrices{3} = FASK_training/max(FASK_training(:))*0.1;

% Find optimal beta values
Betas = zeros(10,3);
for i = 1:10
    [Betas(i,:), ~] = findBetas(input_matrices, 'Target', FC, ...
        'Atlas', atlas_name, 'State', state);
end

%% Generalize to testing set

Betas = importdata(fullfile('fig2(FCrec)data', 'MMP_resting_betas.mat'));
Beta = median(Betas);

% Load data
FC = importdata(fullfile('fig2(FCrec)data', 'MMP_FC_group2.mat'));
rDCM = importdata(fullfile('fig2(FCrec)data', 'MMP_rDCM_group2.mat')); rDCM = rDCM/max(rDCM(:))*0.1;
VAR = importdata(fullfile('fig2(FCrec)data', 'MMP_VAR_group2.mat')); VAR = VAR/max(VAR(:))*0.1;
FASK = importdata(fullfile('fig2(FCrec)data', 'MMP_FASK_group2.mat')); FASK = FASK/max(FASK(:))*0.1;

iEC = (Beta(1)*rDCM + Beta(2)*VAR +Beta(3)*FASK)/3;

% Constants
maxiter = 100;

% Initialization
methods = {'rDCM', 'VAR', 'FASK', 'iEC'};
results = zeros(maxiter, length(methods));

% Parallel computation
parfor i = 1:maxiter
    iter_results = zeros(1, length(methods));
    
    iter_results(1) = FC_approx_Hopf(rDCM, 2, FC, atlas_name, state);
    iter_results(2) = FC_approx_Hopf(VAR, 3, FC, atlas_name, state);
    iter_results(3) = FC_approx_Hopf(FASK, 10, FC, atlas_name, state);
    iter_results(4) = FC_approx_Hopf(iEC, 1, FC, atlas_name, state);
    
    results(i, :) = iter_results;
end

%% Simple Visualization of Figure 2h
% This script generates boxplots for correlation values from Figure 2h.
% One should use BoxPlotPro from MATLAB FileExchange to reproduce figures
% from the manuscript with the same quality.

data = importdata(fullfile('fig2(FCrec)data', 'MMP_iEC_group_testing.mat'));

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

ylim([0.4 0.9])

