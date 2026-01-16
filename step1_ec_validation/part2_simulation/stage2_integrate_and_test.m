%% stage2_integrate_and_test.m
% This script integrates effective connectivity results from multiple algorithms
% and validates the iEC framework against ground truth directed structural connectivity.
% 
% Main workflow:
% 1. Beta estimation: Uses training set (50 iterations) to estimate integration weights
%    for all 8 algorithms and subset of 3 algorithms (VAR, FASK, LiNGAM)
% 2. Visualization: Creates box plots showing beta coefficient distributions
% 3. Validation testing: Tests integrated EC against ground truth using test set (50 iterations)
% 4. Statistical analysis: Compares iEC performance vs individual algorithms using ranksum tests
% 5. Results visualization: Creates correlation box plots with significance markers
%
% Author: Younghyun Oh
% Date: 2025-09-12
% Version: 1.0

% Clear workspace
clear; close all; clc

%-------------------------------------------------------------------------%
%% Integrate EC results
%-------------------------------------------------------------------------%

% Load directed SC and EC results
main_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy';
result_dir = fullfile(main_dir, 'step1_ec_validation/part2_simulation/results');
Js = importdata(fullfile(result_dir, 'directed_SC.mat'));
ec_results = importdata(fullfile(result_dir, 'ec_results.mat'));

% List of atlases
atlases = {'Schaefer', 'MMP'};

% Define algorithm names
algorithms = {'rdcm', 'var', 'fask', 'ccd', 'boss', 'lingam', 'grasp', 'patel'};
num_algorithms = length(algorithms);

% Bayes algorithm names
bayes_algs = {'fask', 'ccd', 'boss', 'lingam', 'grasp', 'patel'};

% Training sample
n_train = 50;
n_test = 50;

% Integrate EC results
EC_betas = zeros(num_algorithms,n_train,length(atlases));
EC_betas_subset = zeros(3,n_train,length(atlases)); % Array for var,fask,lingam betas

for atlas = 1:length(atlases)
    for iter = 1:n_train
        % Extract and normalize matrices
        matrices = cell(num_algorithms, 1);
        matrices_subset = cell(3, 1); % New cell array for subset
        
        for alg_idx = 1:num_algorithms
            matrices{alg_idx} = ec_results(atlas,iter).(algorithms{alg_idx});
            matrices{alg_idx} = matrices{alg_idx}/max(matrices{alg_idx}(:));
        end
        
        % Extract subset matrices (var,fask,lingam)
        matrices_subset{1} = matrices{2}; % var
        matrices_subset{2} = matrices{3}; % fask
        matrices_subset{3} = matrices{6}; % lingam

        % log-linear transform of target matrix 
        target = Js{atlas,iter};
        
        % Estimate beta values for all algorithms
        Betas = findBetas_FLN(matrices, target);
        EC_betas(:,iter,atlas) = Betas;
        
        % Estimate beta values for subset
        Betas_subset = findBetas_FLN(matrices_subset, target);
        EC_betas_subset(:,iter,atlas) = Betas_subset;
    end
end

% % Save results
% save(fullfile(result_dir, 'EC_betas.mat'), 'EC_betas', 'EC_betas_subset');

%-------------------------------------------------------------------------%
%% Create box plot for all beta values
%-------------------------------------------------------------------------%
close all
% Create figure
fig = figure('Position', [100 100 700 200]);

% Choose which atlas to plot (1 for Schaefer, 2 for MMP)
atlas_to_plot = 2;

% Define custom x positions for the boxes
x_positions = 1:num_algorithms;

% Get data for the selected atlas
plot_data_all = squeeze(EC_betas(:,:,atlas_to_plot));

% Create box plot
h = boxplot(plot_data_all', 'Colors', 'k', 'Symbol', '', ...
    'Labels', algorithms);

% Customize box plot appearance
set(h, 'LineWidth', 1.5);
set(findobj(gca, 'type', 'line', 'Tag', 'Median'), 'Color', 'k', 'LineWidth', 2);

% Define colors for each algorithm
colors_all = zeros(num_algorithms, 3);

% Assign pastel blue for rdcm and var
colors_all(1:2,:) = repmat([0.6 0.8 0.9], 2, 1);  % Pastel blue for rdcm, var

% Assign pastel yellow for all others
colors_all(3:end,:) = repmat([0.9 0.8 0.6], num_algorithms-2, 1);  % Pastel yellow for others

% Fill the boxes with colors
hGroups = findobj(gca, 'Tag', 'Box');
for j = 1:length(hGroups)
    box_idx = length(hGroups) - j + 1;  % Reverse the index
    patch(get(hGroups(j), 'XData'), get(hGroups(j), 'YData'), colors_all(box_idx,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

hold on;

% Add scatter plot points with jitter
jitterAmount = 0.15;
for i = 1:num_algorithms
    curr_data = plot_data_all(i,:);
    x_pos = repmat(x_positions(i), length(curr_data), 1) + ...
            (rand(length(curr_data), 1) * jitterAmount - jitterAmount/2);
    scatter(x_pos, curr_data, 40, [0.3 0.3 0.3], 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
end

% Statistical comparison of beta values
row_means = mean(plot_data_all, 1)';
col_diff = plot_data_all - row_means';

% Perform one-tailed t-test for each column vs zero
nCols = size(plot_data_all, 1);
pvals = zeros(1, nCols);

for j = 1:nCols
    [~, pvals(j)] = ttest(col_diff(j, :), 0, 'Tail', 'right');
end

% Apply FDR correction
alpha = 0.05;
[~, ~, ~, pvals_corrected] = fdr_bh(pvals, alpha, 'pdep', 'yes');
significant_idx = find(pvals_corrected < alpha);

% Add asterisks to significant boxes
for i = 1:length(significant_idx)
    idx = significant_idx(i);
    alg_y = max(plot_data_all(idx,:)) + 0.03;
    text(idx, alg_y, '*', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
end

% Customize plot appearance
box on
grid on

% Adjust y-axis limits to show data clearly
ylim_current = ylim;
ylim([ylim_current(1)-0.05, ylim_current(2)+0.07]);
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 5))
% set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])

switch atlas_to_plot
    case 1
        parcellation = 'Schaefer100';
    case 2
        parcellation = 'MMP360';
end

% exportgraphics(fig, fullfile(result_dir, [parcellation, '_all_algorithms_boxplot.png']), 'Resolution', 1200);

%-------------------------------------------------------------------------%
%% Integrate EC results and test against ground truth
%-------------------------------------------------------------------------%

% Initialize arrays to store correlation values
EC_corr_values = zeros(n_test, 10, 2);  % Changed from 11 to 10 columns after removing gc

% Calculate mean betas
Betas = squeeze(mean(EC_betas,2));

% Calculate mean betas for subset
Betas_subset = squeeze(mean(EC_betas_subset,2));

% Iterate over the test set
for atlas = 1:2
    for iter = 1:n_test
        % Extract and normalize matrices for test set
        matrices_test = cell(num_algorithms, 1);
        for alg_idx = 1:num_algorithms
            matrices_test{alg_idx} = ec_results(atlas,iter+n_train).(algorithms{alg_idx});
            matrices_test{alg_idx} = matrices_test{alg_idx}/max(matrices_test{alg_idx}(:));
        end
        
        target = 1.2*Js{atlas,iter+n_train}.^0.3;
        nonzero_idx = target ~= 0;

        % Compute correlations
        corr_values = zeros(1, num_algorithms + 2); % +2 for iec and iec2
        
        % Direct correlations for non-Bayes algorithms
        for i = 1:2 % rdcm, var (changed from 3 to 2 after removing gc)
            corr_values(i) = corr(target(nonzero_idx), matrices_test{i}(nonzero_idx));
        end
        
        % Correlations for Bayes algorithms using non-zero elements
        for i = 3:num_algorithms % fask through patel (indices shifted after removing gc)
            mask = matrices_test{i} ~= 0;
            corr_values(i) = corr(target(mask), matrices_test{i}(mask));
        end
        
        % Integrate EC results using all algorithms
        iec = zeros(size(target));
        for i = 1:num_algorithms
            iec = iec + matrices_test{i}*Betas(i,atlas);
        end
        iec = iec/num_algorithms;
        corr_values(num_algorithms + 1) = corr(target(nonzero_idx), iec(nonzero_idx));

        % Integrate EC results using subset (var, fask, lingam)
        iec2 = matrices_test{2}*Betas_subset(1,atlas) + ...  % var
               matrices_test{3}*Betas_subset(2,atlas) + ...  % fask (index changed from 4 to 3)
               matrices_test{6}*Betas_subset(3,atlas);       % lingam (index changed from 7 to 6)
        iec2 = iec2/3;
        corr_values(num_algorithms + 2) = corr(target(nonzero_idx), iec2(nonzero_idx));

        EC_corr_values(iter,:,atlas) = corr_values;
    end
end

% Save results
save(fullfile(result_dir, 'EC_corr_values.mat'), 'EC_corr_values');

%-------------------------------------------------------------------------%
%% Correlation box plot (for single atlas using all algorithms)
%-------------------------------------------------------------------------%
close all

% Configure figure size
figure('Position', [100 100 900 250]);
atlas = 2;

% Define algorithm names including integrated EC
all_algorithm_names = [algorithms, {'iec','iec2'}];

% Create box plot
h = boxplot(EC_corr_values(:,1:end,atlas),'Colors','k','Symbol','','Labels',all_algorithm_names);
hold on

% Customize box plot appearance
set(h, 'LineWidth', 1.5);
set(findobj(gca, 'type', 'line', 'Tag', 'Median'), 'Color', 'k', 'LineWidth', 2);

% Define colors for each group (same as beta values plot plus new color for iEC)
color_var = [0.6 0.8 0.9];       % Pastel blue for var
color_others = [0.9 0.8 0.6];    % Pastel yellow/orange for others
color_iec = [0.8 0.6 0.9];       % Pastel purple for iEC

% Fill the boxes with different colors based on groups
for i = 1:length(all_algorithm_names)
    if i == 1 || i == 2  % 
        box_color = color_var;
    elseif i == 9 || i == 10 
        box_color = color_iec;
    else  % rest of algorithms
        box_color = color_others;
    end
    patch(get(h(5,i), 'XData'), get(h(5,i), 'YData'), box_color, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

% Add individual data points
for i = 1:length(all_algorithm_names)
    % Get data points for current algorithm
    curr_data = EC_corr_values(:,i,atlas);
    
    % Add jitter to x-coordinates
    x = i + (randn(length(curr_data),1)*0.1);
    
    % Plot scattered points
    scatter(x, curr_data, 5, [0.3 0.3 0.3], 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
end

% Find best performing individual algorithm
[~, best_alg_idx] = max(median(EC_corr_values(:,1:8,atlas)));

% Get median values for significance testing
iec_median = median(EC_corr_values(:,9,atlas));
iec2_median = median(EC_corr_values(:,10,atlas));
best_alg_median = median(EC_corr_values(:,best_alg_idx,atlas));

% Perform statistical tests
p_iec = ranksum(EC_corr_values(:,9,atlas), EC_corr_values(:,best_alg_idx,atlas), 'tail', 'both');
p_iec2 = ranksum(EC_corr_values(:,10,atlas), EC_corr_values(:,best_alg_idx,atlas), 'tail', 'both');

% Add significance brackets and markers
y_max = max([iec_median, iec2_median, best_alg_median]) + 0.1;
best_alg_max = max(EC_corr_values(:,best_alg_idx,atlas));

% Draw horizontal line for iEC variants
line([9 10], [y_max y_max], 'Color', 'k', 'LineWidth', 1.5);

% Draw main horizontal line for the best algorithm
line([best_alg_idx 9.5], [y_max+0.03 y_max+0.03], 'Color', 'k', 'LineWidth', 1.5);
% Draw vertical lines connecting to best algorithm and iEC bracket
line([best_alg_idx best_alg_idx], [best_alg_max+0.03 y_max+0.03], 'Color', 'k', 'LineWidth', 1.5);
line([mean([9 10]) mean([9 10])], [y_max y_max+0.03], 'Color', 'k', 'LineWidth', 1.5);

% Add significance markers
if p_iec < 0.001 || p_iec2 < 0.001
    text(mean([best_alg_idx 10]), y_max+0.05, '***', 'HorizontalAlignment', 'center', 'FontSize', 12);
elseif p_iec < 0.01 || p_iec2 < 0.01
    text(mean([best_alg_idx 10]), y_max+0.05, '**', 'HorizontalAlignment', 'center', 'FontSize', 12);
elseif p_iec < 0.05 || p_iec2 < 0.05
    text(mean([best_alg_idx 10]), y_max+0.05, '*', 'HorizontalAlignment', 'center', 'FontSize', 12);
end

% Customize plot appearance
box on
grid on

% Adjust y-axis limits to show data and annotations clearly
ylim_current = ylim;
ylim([ylim_current(1)-0.05, y_max+0.15]);  % Increased upper limit to show annotations

% set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 4))

% Save figure
switch atlas
    case 1
        parcellation = 'Schaefer100';
    case 2
        parcellation = 'MMP360';
end

% exportgraphics(gcf, fullfile(result_dir, [parcellation, '_correlation_boxplot.png']), 'Resolution', 1200);