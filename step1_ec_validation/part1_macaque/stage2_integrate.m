%% Stage 2: Integrate EC results
% This script integrates effective connectivity results from multiple algorithms
% to create the integrated effective connectivity (iEC) using macaque data.
%
% Main workflow:
% 1. Loads EC results from both original and deconvolved BOLD signals
% 2. Beta estimation: Uses training samples (9 subjects) to estimate integration weights
%    for all 9 algorithms and VFL subset (VAR, FASK, LiNGAM)
% 3. Optimization: Applies findBetas_FLN to find optimal linear combination weights
% 4. Visualization: Creates box plots showing beta coefficient distributions with significance testing
% 5. Saves integration coefficients for validation analysis
%
% Author: Younghyun Oh
% Date: 2025-03-04
%
%--------------------------------------------------------------------------%
% Clear workspace
clear; close all; clc

% Load FLN
base_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy';
main_dir = fullfile(base_dir, 'step1_ec_validation/part1_macaque');

% Load EC results
result_dir = fullfile(main_dir, 'results');
ec_results = load(fullfile(result_dir, 'ec_results.mat')).ec_results;
ec_results_deconv = load(fullfile(result_dir, 'ec_results_deconv.mat')).ec_results;

%--------------------------------------------------------------------------%
%% Integrate EC results
%--------------------------------------------------------------------------%
EC_betas = zeros(9,9);
EC_betas_vfl = zeros(9,3);

% Training sample
num_iter = 9;
for sub = 1:num_iter
    % Extract matrices
    rdcm = ec_results(sub).rdcm;
    var = ec_results(sub).var;
    gc = ec_results(sub).gc;
    fask = ec_results(sub).fask;
    ccd = ec_results(sub).ccd;
    boss = ec_results(sub).boss;
    lingam = ec_results(sub).lingam;
    grasp = ec_results(sub).grasp;
    patel = ec_results(sub).patel;
    
    % Normalize matrices
    input_matrices = cell(9,1);
    input_matrices{1} = rdcm/max(rdcm(:));
    input_matrices{2} = var/max(var(:));
    input_matrices{3} = gc/max(gc(:));
    input_matrices{4} = fask/max(fask(:));
    input_matrices{5} = ccd/max(ccd(:));
    input_matrices{6} = boss/max(boss(:));
    input_matrices{7} = lingam/max(lingam(:));
    input_matrices{8} = grasp/max(grasp(:));
    input_matrices{9} = patel/max(patel(:));

    % Estimate beta values 
    Betas = findBetas_FLN(input_matrices, FLN);
    EC_betas(sub,:) = Betas;
    
    % Set up matrices for VFL
    vfl_matrices = cell(3,1);
    vfl_matrices{1} = var/max(var(:));
    vfl_matrices{2} = fask/max(fask(:));
    vfl_matrices{3} = lingam/max(lingam(:));

    % Estimate beta values for VFL
    Betas_vfl = findBetas_FLN(vfl_matrices, FLN);
    EC_betas_vfl(sub,:) = Betas_vfl;
end

% Save EC betas
save(fullfile(result_dir, 'EC_betas.mat'), 'EC_betas', 'EC_betas_vfl');

%% Create box plot for beta values
close all

% Create figure
fig = figure('Position', [100 100 700 200]);

% Define algorithm names
algorithms = {'rdcm', 'var', 'gc', 'fask', 'ccd', 'boss', 'lingam', 'grasp', 'patel'};

% Create box plot with individual points
h = boxplot(EC_betas, 'Colors', 'k', 'Symbol', '', 'Labels', algorithms);
hold on

% Customize box plot appearance
set(h, 'LineWidth', 1.5);
set(findobj(gca, 'type', 'line', 'Tag', 'Median'), 'Color', 'k', 'LineWidth', 2);

% Define pastel colors for each group
color_var_gc = [0.6 0.8 0.9];    % Pastel blue for rdcm, var and gc
color_others = [0.9 0.8 0.6];    % Pastel yellow/orange for the rest

% Fill the boxes with different colors based on groups
for i = 1:length(algorithms)
    if i == 1 || i == 2 || i == 3  % rdcm, var and gc
        box_color = color_var_gc;
    else
        box_color = color_others;
    end
    patch(get(h(5,i), 'XData'), get(h(5,i), 'YData'), box_color, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

% Add individual data points
for i = 1:length(algorithms)
    % Get data points for current algorithm
    curr_data = EC_betas(:,i);
    
    % Add jitter to x-coordinates
    x = i + (randn(length(curr_data),1)*0.1);
    
    % Plot scattered points
    scatter(x, curr_data, 40, [0.3 0.3 0.3], 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
end

% Customize plot appearance
box on
grid on
% set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])

% Statistical comparison of beta values
row_means = mean(EC_betas, 2);
col_diff = EC_betas - row_means;

% Perform one-tailed t-test for each column vs zero
nCols = size(EC_betas, 2);
pvals = zeros(1, nCols);

for j = 1:nCols
    [~, pvals(j)] = ttest(col_diff(:, j), 0, 'Tail', 'right');
end

% Apply FDR correction
alpha = 0.05;
[~, ~, ~, pvals_corrected] = fdr_bh(pvals, alpha, 'pdep', 'yes');
significant_idx = find(pvals_corrected < alpha);

% Add asterisks to significant boxes
for i = 1:length(significant_idx)
    idx = significant_idx(i);
    alg_y = max(EC_betas(:,idx)) + 0.03;
    text(idx, alg_y, '*', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
end

% Adjust y-axis limits to show data clearly
ylim_current = ylim;
ylim([ylim_current(1)-0.05, ylim_current(2)+0.07]);
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 5))
set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])

% Save figure with transparent background
% exportgraphics(fig, fullfile(main_dir, 'results/beta_boxplot.png'), 'Resolution', 1200);

