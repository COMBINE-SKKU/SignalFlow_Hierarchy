%% Stage 3: Main analysis and validation
% This script performs comprehensive validation of the iEC framework against
% macaque ground truth data and generates Figure 2 components.
%
% Main analyses:
% 1. iEC validation: Tests integrated EC against FLN ground truth using test subjects
% 2. Statistical comparison: Compares iEC vs individual algorithms with significance testing
% 3. F1 score analysis: Evaluates binary classification performance at different thresholds
% 4. SLN validation: Tests iEC against structural laminar network (feedforward/feedback)
% 5. Deconvolution comparison: Analyzes effects of hemodynamic deconvolution
% 6. Algorithm relationships: Examines correlations and contributions between algorithms
% 7. Supplementary analyses: Edge distributions, permutation tests, detailed visualizations
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
FLN_dir = fullfile(main_dir, 'data/FLN40.mat');
FLN = importdata(FLN_dir);
FLN = 1.2*FLN.^0.3; % log-linear transformation

% Load EC results
result_dir = fullfile(main_dir, 'results');
ec_results = load(fullfile(result_dir, 'ec_results.mat')).ec_results;
ec_results_deconv = load(fullfile(result_dir, 'ec_results_deconv.mat')).ec_results;
EC_betas = load(fullfile(result_dir, 'EC_betas.mat')).EC_betas;
EC_betas_vfl = load(fullfile(result_dir, 'EC_betas.mat')).EC_betas_vfl;

%--------------------------------------------------------------------------%
%% Integrate EC results and test against ground truth
%--------------------------------------------------------------------------%
% Initialize arrays to store correlation values
EC_correlation_values_FLN = zeros(11, 11);
EC_correlation_values_FLN_deconv = zeros(11, 11);

% Calculate mean betas
Betas = median(EC_betas);
Betas_vfl = median(EC_betas_vfl);

% Initialize array to store iEC values
iter = 1; 

% Iterate over the test set for both original and deconvolved results
for sub = 9:19
    % Extract matrices for original
    rdcm = ec_results(sub).rdcm;
    var = ec_results(sub).var;
    gc = ec_results(sub).gc;
    fask = ec_results(sub).fask;
    ccd = ec_results(sub).ccd;
    boss = ec_results(sub).boss;
    lingam = ec_results(sub).lingam;
    grasp = ec_results(sub).grasp;
    patel = ec_results(sub).patel;

    % Normalize matrices for original
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

    % Constructed iEC for original
    ec_results(sub).iec = (input_matrices{1}*Betas(1) + ...
          input_matrices{2}*Betas(2) + input_matrices{3}*Betas(3) + ...
          input_matrices{4}*Betas(4) + input_matrices{5}*Betas(5) + ...
          input_matrices{6}*Betas(6) + input_matrices{7}*Betas(7) + ...
          input_matrices{8}*Betas(8) + input_matrices{9}*Betas(9))/length(Betas);

    % VFL for original
    ec_results(sub).vfl = (input_matrices{2}*Betas_vfl(1) + input_matrices{4}*Betas_vfl(2) + input_matrices{7}*Betas_vfl(3))/3;

    % Calculate correlation with FLN for original
    corr_values = zeros(1, 11);
    
    % Loop through all individual algorithms
    for i = 1:9
        corr_values(i) = corr(FLN(:), input_matrices{i}(:));
    end    

    % Add correlations for integrated methods
    corr_values(10) = corr(FLN(:), ec_results(sub).iec(:));
    corr_values(11) = corr(FLN(:), ec_results(sub).vfl(:));
    
    EC_correlation_values_FLN(iter, :) = corr_values;

    % Repeat the process for deconvolved results
    rdcm_deconv = ec_results_deconv(sub).rdcm;
    var_deconv = ec_results_deconv(sub).var;
    gc_deconv = ec_results_deconv(sub).gc;
    fask_deconv = ec_results_deconv(sub).fask;
    ccd_deconv = ec_results_deconv(sub).ccd;
    boss_deconv = ec_results_deconv(sub).boss;
    lingam_deconv = ec_results_deconv(sub).lingam;
    grasp_deconv = ec_results_deconv(sub).grasp;
    patel_deconv = ec_results_deconv(sub).patel;

    % Normalize matrices for deconvolved
    input_matrices_deconv = cell(9,1);
    input_matrices_deconv{1} = rdcm_deconv/max(rdcm_deconv(:));
    input_matrices_deconv{2} = var_deconv/max(var_deconv(:));
    input_matrices_deconv{3} = gc_deconv/max(gc_deconv(:));
    input_matrices_deconv{4} = fask_deconv/max(fask_deconv(:));
    input_matrices_deconv{5} = ccd_deconv/max(ccd_deconv(:));
    input_matrices_deconv{6} = boss_deconv/max(boss_deconv(:));
    input_matrices_deconv{7} = lingam_deconv/max(lingam_deconv(:));
    input_matrices_deconv{8} = grasp_deconv/max(grasp_deconv(:));
    input_matrices_deconv{9} = patel_deconv/max(patel_deconv(:));

    % Constructed iEC for deconvolved
    ec_results_deconv(sub).iec = (input_matrices_deconv{1}*Betas(1) + ...
          input_matrices_deconv{2}*Betas(2) + input_matrices_deconv{3}*Betas(3) + ...
          input_matrices_deconv{4}*Betas(4) + input_matrices_deconv{5}*Betas(5) + ...
          input_matrices_deconv{6}*Betas(6) + input_matrices_deconv{7}*Betas(7) + ...
          input_matrices_deconv{8}*Betas(8) + input_matrices_deconv{9}*Betas(9))/length(Betas);

    % VFL for deconvolved
    ec_results_deconv(sub).vfl = (input_matrices_deconv{2}*Betas_vfl(1) + input_matrices_deconv{4}*Betas_vfl(2) + input_matrices_deconv{7}*Betas_vfl(3))/3;

    % Calculate correlation with FLN for deconvolved
    corr_values_deconv = zeros(1, 11);

    % Loop through all individual algorithms
    for i = 1:9
        corr_values_deconv(i) = corr(FLN(:), input_matrices_deconv{i}(:));
    end
    
    corr_values_deconv(10) = corr(FLN(:), ec_results_deconv(sub).iec(:));
    corr_values_deconv(11) = corr(FLN(:), ec_results_deconv(sub).vfl(:));
    EC_correlation_values_FLN_deconv(iter, :) = corr_values_deconv;

    iter = iter + 1;
end

%% Statistical comparison of iEC (top 3) vs best individual algorithm
% Define algorithm names including integrated EC
algorithms = {'rdcm', 'var', 'gc', 'fask', 'ccd', 'boss', 'lingam', 'grasp', 'patel','iec', 'vfl'};

% Calculate median correlation for each algorithm
median_corr = median(EC_correlation_values_FLN);

% Find the best individual algorithm (excluding integrated methods)
[best_ind_val, best_ind_idx] = max(median_corr(1:9));
best_ind_alg = algorithms{best_ind_idx};

% Get correlation values for VFL (top 3 integration)
vfl_corr = EC_correlation_values_FLN(:, 11);
best_ind_corr = EC_correlation_values_FLN(:, best_ind_idx);

% Perform Wilcoxon rank-sum test
% Testing if VFL is significantly better than the best individual algorithm
[p_val, h_val] = signrank(vfl_corr, best_ind_corr, 'tail', 'both');

%% Visualize correlation values
% close all

% Configure figure size
fig2 = figure('Position', [100 100 950 300]);

% Define algorithm names including integrated EC
algorithms = {'rdcm', 'var', 'gc', 'fask', 'ccd', 'boss', 'lingam', 'grasp', 'patel','iec','vfl'};

data = EC_correlation_values_FLN;
% Create box plot with individual points
h = boxplot(data, 'Colors', 'k', 'Symbol', '', 'Labels', algorithms);
hold on

% Customize box plot appearance
set(h, 'LineWidth', 1.5);
set(findobj(gca, 'type', 'line', 'Tag', 'Median'), 'Color', 'k', 'LineWidth', 2);

% Define pastel colors for each group
color_var_gc = [0.6 0.8 0.9];    % Pastel blue for rdcm, var and gc
color_others = [0.9 0.8 0.6];    % Pastel yellow/orange for the rest
color_iec = [0.8 0.6 0.9];       % Pastel purple for integrated versions

% Fill the boxes with different colors based on groups
for i = 1:length(algorithms)
    if i == 1 || i == 2 || i == 3  % rdcm, var and gc
        box_color = color_var_gc;
    elseif i == 10 || i == 11  % integrated versions
        box_color = color_iec;
    else  % rest of algorithms
        box_color = color_others;
    end
    patch(get(h(5,i), 'XData'), get(h(5,i), 'YData'), box_color, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

% Add individual data points
for i = 1:length(algorithms)
    % Get data points for current algorithm
    curr_data = data(:,i);
    
    % Add jitter to x-coordinates
    x = i + (randn(length(curr_data),1)*0.1);
    
    % Plot scattered points
    scatter(x, curr_data, 40, [0.3 0.3 0.3], 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
end

% Find best performing individual algorithm
[~, best_alg_idx] = max(median(data(:,1:9)));

% Get median values for significance testing
iec_median = median(data(:,10));
vfl_median = median(data(:,11));
best_alg_median = median(data(:,best_alg_idx));

% Perform statistical tests
p_iec = ranksum(data(:,10), data(:,best_alg_idx), 'tail', 'both');
p_vfl = ranksum(data(:,11), data(:,best_alg_idx), 'tail', 'both');

% Add significance brackets and markers
y_max = max([iec_median, vfl_median, best_alg_median]) + 0.1;
best_alg_max = max(data(:,best_alg_idx));

% Draw horizontal line for iEC variants
line([10 11], [y_max y_max], 'Color', 'k', 'LineWidth', 1.5);

% Draw main horizontal line for the best algorithm
line([best_alg_idx 10.5], [y_max+0.03 y_max+0.03], 'Color', 'k', 'LineWidth', 1.5);
% Draw vertical lines connecting to best algorithm and iEC bracket
line([best_alg_idx best_alg_idx], [best_alg_max+0.03 y_max+0.03], 'Color', 'k', 'LineWidth', 1.5);
line([mean([10 11]) mean([10 11])], [y_max y_max+0.03], 'Color', 'k', 'LineWidth', 1.5);

% Add significance markers
if p_iec < 0.001 || p_vfl < 0.001
    text(mean([best_alg_idx 11]), y_max+0.05, '***', 'HorizontalAlignment', 'center', 'FontSize', 12);
elseif p_iec < 0.01 || p_vfl < 0.01
    text(mean([best_alg_idx 11]), y_max+0.05, '**', 'HorizontalAlignment', 'center', 'FontSize', 12);
elseif p_iec < 0.05 || p_vfl < 0.05
    text(mean([best_alg_idx 11]), y_max+0.05, '*', 'HorizontalAlignment', 'center', 'FontSize', 12);
end

% Customize plot appearance
box on
grid on
set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])

% Adjust y-axis limits to show data clearly
ylim_current = ylim;
ylim([ylim_current(1)-0.05, y_max+0.11]);
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 4))

% Save figure with transparent background
% exportgraphics(fig2, fullfile(main_dir, 'results', 'ec_correlation_boxplot.png'), 'Resolution', 1200);

%--------------------------------------------------------------------------%
%% Compute F1 scores with specific thresholds
%--------------------------------------------------------------------------%
close all

% Define specific thresholds
threshold_values = [0.15, 0.3];
algorithms = {'rdcm', 'var', 'gc', 'fask', 'ccd', 'boss', 'lingam', 'grasp', 'patel', 'vfl'};

% Initialize storage for F1 scores
% Structure: thresholds x subjects x algorithms
F1_scores = zeros(length(threshold_values), 10, length(algorithms));

% Compute F1 scores for each threshold
for i = 1:length(threshold_values)
    binary_FLN = threshold_proportional(FLN, threshold_values(i)) > 0;
    
    iter = 1;  % Initialize counter for proper indexing
    for j = 9:19
        for k = 1:length(algorithms)-1  % Regular algorithms
            alg = algorithms{k};
            binary_pred = threshold_proportional(ec_results(j).(alg), threshold_values(i)) > 0;
            
            % Compute metrics
            TP = sum(binary_pred & binary_FLN, 'all');
            FP = sum(binary_pred & ~binary_FLN, 'all');
            FN = sum(~binary_pred & binary_FLN, 'all');
            
            % Compute precision and recall
            precision = TP / (TP + FP + eps);
            recall = TP / (TP + FN + eps);
            
            % Compute F1 score
            F1_scores(i, iter, k) = 2 * (precision * recall) / (precision + recall + eps);
        end
        
        % Compute F1 score for VFL
        binary_pred = threshold_proportional(ec_results(j).vfl, threshold_values(i)) > 0;
        TP = sum(binary_pred & binary_FLN, 'all');
        FP = sum(binary_pred & ~binary_FLN, 'all');
        FN = sum(~binary_pred & binary_FLN, 'all');
        
        precision = TP / (TP + FP + eps);
        recall = TP / (TP + FN + eps);
        F1_scores(i, iter, length(algorithms)) = 2 * (precision * recall) / (precision + recall + eps);
        
        iter = iter + 1;  % Increment counter
    end
end

%% Visualize thresholded FLN and F1 scores for each threshold
close all

for i = 1:length(threshold_values)
    % Create thresholded (not binarized) FLN
    thresholded_FLN = threshold_proportional(FLN, threshold_values(i));
    
    % Create FLN visualization figure
    fig_fln = figure('Position', [100 100 250 250]);
    imagesc(thresholded_FLN);
    colormap(generateColorMap(thresholded_FLN(:),100));
    set(gca, 'XTick', [], 'YTick', []);
    axis square
    exportgraphics(fig_fln, fullfile(main_dir, sprintf('results/thresholded_FLN_%g.png', threshold_values(i))), 'Resolution', 1200);
    close(fig_fln);
    
    % Create F1 score box plot
    F1_data_for_plot = squeeze(F1_scores(i, :, :));
    
    % Find algorithm with highest median F1 score (excluding VFL)
    median_F1 = median(F1_data_for_plot(:,1:end-1));
    [~, best_alg_idx] = max(median_F1);
    vfl_idx = size(F1_data_for_plot, 2);
    
    % Perform statistical test (Wilcoxon signed rank test, one-sided)
    [p_val, ~] = signrank(F1_data_for_plot(:,best_alg_idx), F1_data_for_plot(:,vfl_idx), 'tail', 'both');
    
    % Create figure
    fig_f1 = figure('Position', [100 100 900 300]);
    
    % Create box plot with individual points
    h = boxplot(F1_data_for_plot, 'Colors', 'k', 'Symbol', '', 'Labels', algorithms);
    hold on
    
    % Customize box plot appearance
    set(h, 'LineWidth', 1.5);
    set(findobj(gca, 'type', 'line', 'Tag', 'Median'), 'Color', 'k', 'LineWidth', 2);
    
    % Define colors for each group
    color_var_gc = [0.6 0.8 0.9];    % Pastel blue for rdcm, var and gc
    color_others = [0.9 0.8 0.6];    % Pastel yellow/orange for others
    color_iec = [0.8 0.6 0.9];       % Pastel purple for VFL
    
    % Fill the boxes with different colors
    for j = 1:length(algorithms)
        if j <= 3  % rdcm, var and gc
            box_color = color_var_gc;
        elseif j == length(algorithms)  % VFL
            box_color = color_iec;
        else  % rest of algorithms
            box_color = color_others;
        end
        patch(get(h(5,j), 'XData'), get(h(5,j), 'YData'), box_color, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    end
    
    % Add individual data points
    for j = 1:length(algorithms)
        curr_data = F1_data_for_plot(:,j);
        x = j + (randn(length(curr_data),1)*0.1);
        scatter(x, curr_data, 40, [0.3 0.3 0.3], 'filled', ...
            'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
    end
    
    % Add horizontal line and significance indicator
    y1 = median(F1_data_for_plot(:,best_alg_idx));
    y2 = median(F1_data_for_plot(:,vfl_idx));
    y_line = max([y1, y2]) + 0.1;
    plot([best_alg_idx, vfl_idx], [y_line y_line], 'k-', 'LineWidth', 1.5);
    
    x_text = mean([best_alg_idx, vfl_idx]);
    if p_val < 0.001
        text(x_text, y_line + 0.05, '***', 'HorizontalAlignment', 'center', 'FontSize', 14);
    elseif p_val < 0.01
        text(x_text, y_line + 0.05, '**', 'HorizontalAlignment', 'center', 'FontSize', 14);
    elseif p_val < 0.05
        text(x_text, y_line + 0.05, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    else
        text(x_text, y_line + 0.05, 'n.s.', 'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    
    % Customize plot appearance
    box on
    grid on
    set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
    ylim_current = ylim;
    ylim([ylim_current(1)-0.05, y_line + 0.15]);
    set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 4))
    
    % Save F1 score figure
    exportgraphics(fig_f1, fullfile(main_dir, sprintf('results/F1_scores_threshold_%g.png', threshold_values(i))), 'Resolution', 1200);
    close(fig_f1);
end

%--------------------------------------------------------------------------%
%% Create scatter plot for VFL-based iEC vs FLN
%--------------------------------------------------------------------------%
close all

% Find the median performance for VFL
corr_vals = EC_correlation_values_FLN(:,end); % Get VFL correlations
[~, median_idx] = min(abs(corr_vals - median(corr_vals)));
median_sub = median_idx + 8; % Convert to actual subject index

var_median = ec_results(median_sub).var;
fask_median = ec_results(median_sub).fask;
lingam_median = ec_results(median_sub).lingam;

var_median = var_median/max(var_median(:));
fask_median = fask_median/max(fask_median(:));
lingam_median = lingam_median/max(lingam_median(:));

vfl_median = (var_median*Betas_vfl(1) + ...
              fask_median*Betas_vfl(2) + ...
              lingam_median*Betas_vfl(3))/3;

iec_vfl = vfl_median/max(vfl_median(:));

% Create figure
figure('Position', [100 100 400 270]);

% Create scatter plot
scatter(iec_vfl(:), FLN(:), 40, [0.5 0.5 0.5], 'filled', ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
hold on

% Add regression line
p = polyfit(iec_vfl(:), FLN(:), 1);
x_range = [min(iec_vfl(:)) max(iec_vfl(:))];
y_fit = polyval(p, x_range);
plot(x_range, y_fit, '-r', 'LineWidth', 1.5);

% Calculate correlation
corr_val = corr(FLN(:), iec_vfl(:));

% Customize appearance
box on
grid on
% set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
set(gca, 'XTick', linspace(min(get(gca, 'XTick')), max(get(gca, 'XTick')), 5))
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 5))

% Set axis limits
xlim([min(iec_vfl(:)) max(iec_vfl(:))])
ylim([min(FLN(:)) max(FLN(:))])

% Save figure
exportgraphics(gcf, fullfile(main_dir, 'results/bei_recovery.png'), 'Resolution', 1200);

%% Create histogram showing VAR, FASK, LiNGAM contributions
close all
fig4 = figure('Position', [100 100 500 200]);

% Create histogram with custom binning
vfl_median_values = iec_vfl(:);
num_bins = 30;
histogram_edges = linspace(min(vfl_median_values)-0.1, max(vfl_median_values)+0.1, num_bins+1);
histogram_centers = (histogram_edges(1:end-1) + histogram_edges(2:end))/2;
bin_width = histogram_edges(2) - histogram_edges(1);

% Count occurrences in each bin and determine dominant contribution
bin_counts = zeros(num_bins, 1);
bin_types = zeros(num_bins, 1); % 1 for VAR, 2 for FASK, 3 for LiNGAM

% Get total number of samples
total_samples = length(vfl_median_values);

for i = 1:num_bins
    bin_mask = vfl_median_values >= histogram_edges(i) & vfl_median_values < histogram_edges(i+1);
    bin_counts(i) = sum(bin_mask) / total_samples; % Convert to probability
    
    if bin_counts(i) > 0
        % Get absolute contributions in this bin
        var_values = abs(var_median(bin_mask));
        fask_values = abs(fask_median(bin_mask));
        lingam_values = abs(lingam_median(bin_mask));
        
        % Calculate relative contributions (normalize by sum)
        total_contrib = mean(var_values) + mean(fask_values) + mean(lingam_values);
        var_contribution = mean(var_values) / total_contrib;
        fask_contribution = mean(fask_values) / total_contrib;
        lingam_contribution = mean(lingam_values) / total_contrib;
        
        % Find dominant algorithm (highest contribution)
        [~, max_idx] = max([var_contribution, fask_contribution, lingam_contribution]);
        bin_types(i) = max_idx;
    else
        % Skip bins with no samples
        bin_types(i) = NaN;  % Use NaN instead of 0 for empty bins
    end
end

% Create separate bars for each type
hold on
% Colors from the provided hex codes
colors = {[178/255 211/255 154/255],    % 99BC85 - Light green for VAR
          [251/255 228/255 214/255],    % FBE4D6 - Light peach for FASK
          [255/255 183/255 130/255]};   % FFC1B4 - Light coral for LiNGAM

% Set custom bar width (as a fraction of bin width)
bar_width =  bin_width;

% Plot only non-empty bins
for type = 1:3
    type_mask = bin_types == type;  % NaN values will be excluded automatically
    if any(type_mask)
        % Calculate bar positions
        x = histogram_centers(type_mask);
        y = bin_counts(type_mask);
        
        % Create bars with custom width
        for j = 1:length(x)
            rectangle('Position', [x(j)-bar_width/2, 1e-6, bar_width, y(j)], ...
                     'FaceColor', colors{type}, ...
                     'EdgeColor', 'k', ...
                     'LineWidth', 0.5);
        end
    end
end

% Customize appearance
box off
set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
set(gca, 'XTick', [0.1 0.4 0.7 1])
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 4))

% Adjust axis limits
% xlim([0.395 1.015])
% ylim([0 0.01])

% Save figure
% exportgraphics(fig4, fullfile(main_dir, 'results/contribution_histogram_inset.png'), 'Resolution', 1200);

%% Visualize the distribution of non-zero FLN values using a histogram
non_zero_FLN = FLN(FLN > 0);
figure('Position', [100, 100, 270, 70]);
histogram(non_zero_FLN, 'Normalization', 'probability', 'DisplayStyle', 'bar', 'FaceColor', [0.2, 0.6, 0.8]);
set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'XTick', 0:0.3:max(non_zero_FLN), 'YTick', 0:0.1:1);
set(gca, 'XTickLabel', [], 'YTickLabel', []);
box off;
exportgraphics(gcf, fullfile(main_dir, 'results', 'fln_distribution.png'), 'Resolution', 1200);


%--------------------------------------------------------------------------%
%% iEC and SLN
%--------------------------------------------------------------------------%

% Load SLN
SLN_dir = fullfile(main_dir, 'data/SLN40.mat');
SLN = importdata(SLN_dir);

% Define indices for feedback and feedforward
FBindex = SLN < 0.5 ;
FFindex = SLN >= 0.5 ;

% Initialize results storage
ratio_array = zeros(11, 4);

for sub = 9:num_subjects

    % Extract EC values using FF and FB indices
    iec_FF = ec_results(sub).vfl(FFindex);
    iec_FB = ec_results(sub).vfl(FBindex);

    % Calculate ratio of negative and positive values
    ratio_neg_FF = sum(iec_FF < 0) / length(iec_FF);
    ratio_pos_FF = sum(iec_FF > 0) / length(iec_FF);
    ratio_neg_FB = sum(iec_FB < 0) / length(iec_FB);
    ratio_pos_FB = sum(iec_FB > 0) / length(iec_FB);

    ratio_array(sub-8, :) = [ratio_neg_FF, ratio_pos_FF, ratio_neg_FB, ratio_pos_FB];
end

% Generate box plot for FF and FB with neg_VAR, pos_VAR
% Reshape the data vectors for each measure
data = reshape(ratio_array, [], 1);

% Generate group and legend labels
group = repmat({'neg', 'pos', 'neg', 'pos'}, 11, 1);
group = reshape(group, [], 1);
legend = repmat({'FF', 'FF', 'FB', 'FB'}, 11, 1);
legend = reshape(legend, [], 1);
combined_group = strcat(group, '_', legend);

%% Create the boxplot and set box colors
close all

% Configure figure size
figure('Position', [100, 100, 250, 350]);

% Define custom x positions for the boxes
x_positions = [0.7 1.3 2.7 3.3];  % Closer spacing within FF/FB pairs
box_width = 0.5;  % Control box width directly

% Create boxplot with custom positions and width
h = boxplot(data, combined_group, 'Colors', 'k', 'Symbol', '', ...
    'positions', x_positions, ...
    'widths', box_width);  % Set uniform width for all boxes

% Customize box plot appearance
set(h, 'LineWidth', 1.5);
set(findobj(gca, 'type', 'line', 'Tag', 'Median'), 'Color', 'k', 'LineWidth', 2);

% Define colors for each box
col_vals = [...
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf'];
col_vals = reshape(hex2dec(col_vals), [3 numel(col_vals)/6])' ./ 255;
P = size(col_vals, 1);
colors = interp1(1:size(col_vals, 1), col_vals, linspace(1, P, 4), 'linear');

% Set the colors for each box
hGroups = findobj(gca, 'Tag', 'Box');
for j = 1:length(hGroups)
    patch(get(hGroups(j), 'XData'), get(hGroups(j), 'YData'), colors(j, :), 'FaceAlpha', 0.8);
end

hold on;

% Add scatter plot points with adjusted positions
unique_groups = {'neg_FF', 'pos_FF', 'neg_FB', 'pos_FB'};
jitterAmount = 0.15;

for i = 1:length(unique_groups)
    % Use custom x positions for scatter points
    x_pos = repmat(x_positions(i), sum(strcmp(combined_group, unique_groups{i})), 1) + ...
            (rand(sum(strcmp(combined_group, unique_groups{i})), 1) * jitterAmount);
    y_pos = data(strcmp(combined_group, unique_groups{i}));
    scatter(x_pos, y_pos, 40, [0.5, 0.5, 0.5], 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
end

% Customize plot appearance
box on
grid on
set(gca, 'YTick', linspace(0, 1.15, 5))

% Set custom x-ticks in the middle of pairs
set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
set(gca, 'XTick', [mean(x_positions(1:2)), mean(x_positions(3:4))])
% set(gca, 'XTickLabel', {'FF', 'FB'})


% Adjust axis limits
xlim([0 4])
ylim([0 1.15]);

% Save figure
exportgraphics(gcf, fullfile(main_dir, 'results/ratio_boxplot.png'), 'Resolution', 1200);

%--------------------------------------------------------------------------%
%% Permutation test using iec_vfl
%--------------------------------------------------------------------------%
close all

% Calculate actual correlation
iec_emp = iec_vfl(:);
lingam_emp = lingam_median(:);
target = FLN(:);

actual_corr_iec = corr(iec_emp, target);
actual_corr_lingam = corr(lingam_emp, target);
corr_diff = actual_corr_iec - actual_corr_lingam;

% Generate null distribution
n_perms = 1000;
null_dist = zeros(n_perms, 1);

% Use parfor for parallel computation
null_corr_diff = zeros(n_perms,1);
parfor i = 1:n_perms
    var_null = null_model_dir_sign(var_median,5,1);
    fask_null = null_model_dir_sign(fask_median,5,1);
    lingam_null = null_model_dir_sign(lingam_median);
    iec_null = var_null*Betas_vfl(1) + fask_null*Betas_vfl(2) + lingam_null*Betas_vfl(3);
    null_corr_iec = corr(iec_null(:), target);
    null_corr_lingam = corr(lingam_null(:), target);
    null_corr_diff(i) = null_corr_iec - null_corr_lingam;
end

% One-sided p-value: fraction of permutations with a difference >= observed
p_value = sum(null_corr_diff >= corr_diff) / n_perms;

% Display results
fprintf('Observed difference in correlation (iEC - LiNGAM) = %.4f\n', corr_diff);
fprintf('Permutation test p-value = %.4f (n=%d)\n', p_value, n_perms);

% Optional: visualize the null distribution
figure('Position', [100 100 500 300]);
histogram(null_corr_diff, 30, 'Normalization', 'probability', ...
    'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k', 'LineWidth', 1);
hold on;
xline(corr_diff,  'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');

% Customize appearance
box on
grid on
% set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
set(gca, 'XTick', linspace(min(get(gca, 'XTick')), max(get(gca, 'XTick')), 5))
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 4))

% Save figure
% exportgraphics(gcf, fullfile(main_dir, 'supplementary/null_increament_test.png'), 'Resolution', 1200);

%% Create figure
figure('Position', [100 100 500 300]);

% Create histogram
histogram(null_dist, 30, 'Normalization', 'probability', ...
    'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k', 'LineWidth', 1);
hold on

% Add vertical line for actual correlation
line([actual_corr actual_corr], ylim, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');

% Calculate p-value
p_value = sum(null_dist >= actual_corr) / length(null_dist);

% Customize appearance
box on
grid on
% set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
set(gca, 'XTick', linspace(min(get(gca, 'XTick')), max(get(gca, 'XTick')), 5))
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 4))

% Adjust axis limits
xlim([min(null_dist)-0.05 max([max(null_dist) actual_corr])+0.05])

% Save figure
% exportgraphics(gcf, fullfile(main_dir, 'supplementary/permutation_test.png'), 'Resolution', 1200);


%--------------------------------------------------------------------------%
%% Compare original and deconvolved EC results
%--------------------------------------------------------------------------%

% Initialize storage for group averages and correlations
group_avg_orig = struct();
group_avg_deconv = struct();
correlations = struct();

% Compute group averages and correlations for each algorithm
for alg = algorithms
    alg_name = alg{1};
    
    % Extract matrices for all subjects
    orig_matrices = cat(3, ec_results.(alg_name));
    deconv_matrices = cat(3, ec_results_deconv.(alg_name));
    
    % Compute group averages
    group_avg_orig.(alg_name) = mean(orig_matrices, 3);
    group_avg_deconv.(alg_name) = mean(deconv_matrices, 3);
    
    % Compute correlation between original and deconvolved group averages
    orig_vec = group_avg_orig.(alg_name)(:);
    deconv_vec = group_avg_deconv.(alg_name)(:);
    correlations.(alg_name) = corr(orig_vec, deconv_vec);
end

% Display correlations
fprintf('\nCorrelations between original and deconvolved group averages:\n');
fprintf('--------------------------------------------------------\n');
for alg = algorithms
    fprintf('%s: %.3f\n', alg{1}, correlations.(alg{1}));
end

% Create bar plot of correlations
figure('Position', [100 100 950 300]);

% Extract correlation values in order
corr_values = zeros(1, length(algorithms));
for i = 1:length(algorithms)
    corr_values(i) = correlations.(algorithms{i});
end

% Create bar plot
b = bar(corr_values, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none', 'BarWidth', 0.6);
hold on;

% Customize plot appearance
box on;
grid on;
set(gca, 'XTick', 1:length(algorithms), ...
    'XTickLabel', [], ...  % Remove x tick labels
    'YTickLabel', [], ...  % Remove y tick labels
    'FontSize', 12, ...
    'LineWidth', 1.5, ...
    'GridAlpha', 0.15, ...
    'YTick', linspace(min(corr_values)*0.9, max(corr_values)*1.1, 6)); % Sparser y ticks with 4 points

% Set y-axis limits with some padding
ylim([min(corr_values)*0.9, max(corr_values)*1.1]);

% Save figure
% exportgraphics(gcf, fullfile(main_dir, 'supplementary/deconv_correlations.png'), ...
%     'Resolution', 1200, 'BackgroundColor', 'white');

%% Scatter plot of iEC before and after deconvolution
close all

% Compute average of VAR, FASK, LiNGAM before and after deconvolution   
avg_VAR_orig = mean(cat(3, ec_results.var), 3);
avg_VAR_deconv = mean(cat(3, ec_results_deconv.var), 3);
avg_FASK_orig = mean(cat(3, ec_results.fask), 3);
avg_FASK_deconv = mean(cat(3, ec_results_deconv.fask), 3);
avg_LINGAM_orig = mean(cat(3, ec_results.lingam), 3);
avg_LINGAM_deconv = mean(cat(3, ec_results_deconv.lingam), 3);

% Compute iEC before and after deconvolution
iEC_orig = avg_VAR_orig*Betas_vfl(1) + avg_FASK_orig*Betas_vfl(2) + avg_LINGAM_orig*Betas_vfl(3);
iEC_deconv = avg_VAR_deconv*Betas_vfl(1) + avg_FASK_deconv*Betas_vfl(2) + avg_LINGAM_deconv*Betas_vfl(3);

iEC_orig = iEC_orig/max(iEC_orig(:));
iEC_deconv = iEC_deconv/max(iEC_deconv(:));

%% Create scatter plot
figure('Position', [100 100 400 300]);

% Convert matrices to vectors for plotting
iEC_orig_vec = iEC_orig(:);
iEC_deconv_vec = iEC_deconv(:);

% Create scatter plot
scatter(iEC_orig_vec, iEC_deconv_vec, 40, [0.4 0.4 0.4], 'filled', ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);

% Add regression line
hold on;
p = polyfit(iEC_orig_vec, iEC_deconv_vec, 1);
x_fit = linspace(min(iEC_orig_vec), max(iEC_orig_vec), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, '--', 'Color', 'r', 'LineWidth', 1.5);

% Calculate correlation
r = corr(iEC_orig_vec, iEC_deconv_vec);

% Customize plot
% axis square;
box on;
grid on;

% Set axis limits with some padding separately for x and y
x_min = min(iEC_orig_vec) * 1.1;
x_max = max(iEC_orig_vec) * 1.1;
y_min = min(iEC_deconv_vec) * 1.1;
y_max = max(iEC_deconv_vec) * 1.1;
xlim([x_min x_max]);
ylim([y_min y_max]);

% Customize appearance
set(gca, 'LineWidth', 1.5, ...
    'XTick', linspace(x_min, x_max, 5), ...
    'YTick', linspace(y_min, y_max, 5), ...
    'TickLength', [0.02 0.02]);
% 
% % Save figure
% exportgraphics(gcf, fullfile(main_dir, 'supplementary/iEC_deconv_comparison.png'), ...
%     'Resolution', 1200, 'BackgroundColor', 'white');

%% Calculate median correlation values for each algorithm
median_corr_original = median(EC_correlation_values_FLN);
median_corr_deconv = median(EC_correlation_values_FLN_deconv);

% Calculate SEM for each algorithm
sem_corr_original = std(EC_correlation_values_FLN) / sqrt(size(EC_correlation_values_FLN, 1));
sem_corr_deconv = std(EC_correlation_values_FLN_deconv) / sqrt(size(EC_correlation_values_FLN_deconv, 1));

% Calculate differences and SEM of differences
diff_corr = median_corr_original - median_corr_deconv;
sem_diff = sqrt(sem_corr_original.^2 + sem_corr_deconv.^2);

% Create bar plot with SEM
% Define algorithm names including integrated EC
algorithms = {'rdcm', 'var', 'gc', 'fask', 'ccd', 'boss', 'lingam', 'grasp', 'patel'};
all_algorithm_names = [algorithms, {'iec', 'vfl'}];

figure('Position', [100 100 950 300]);
bar(diff_corr, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none', 'BarWidth', 0.6);
hold on;
errorbar(1:length(diff_corr), diff_corr, sem_diff, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
set(gca, 'LineWidth', 1.5);
set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 5))
grid on;

% Save figure with transparent background
exportgraphics(gca, fullfile(main_dir, 'supplementary', 'difference_barplot.png'), 'Resolution', 1200);

%--------------------------------------------------------------------------%
%% Compute pairwise correlations between all algorithms
%--------------------------------------------------------------------------%

% Algorithm names
alg_names = {'rdcm','var','gc','fask','ccd','boss','lingam','grasp','patel'};
n_algs = length(alg_names);

% Compute correlations for each subject
all_correlations = zeros(num_subjects, n_algs, n_algs);

for subj = 1:num_subjects
    for i = 1:n_algs
        for j = 1:n_algs
            alg1 = alg_names{i};
            alg2 = alg_names{j};
            
            % Get EC matrices for both algorithms
            ec1 = ec_results(subj).(alg1)(:);
            ec2 = ec_results(subj).(alg2)(:);
            
            % Compute correlation
            all_correlations(subj, i, j) = corr(ec1, ec2);
        end
    end
end

% Compute mean correlation across subjects
mean_corr = squeeze(mean(all_correlations, 1));
sem_corr = squeeze(std(all_correlations, 0, 1)) / sqrt(num_subjects);

%% Create heatmap figure
close all
figure('Position', [100 100 500 500]);

% Create copy of correlation matrix with NaN diagonal and upper triangle for plotting
plot_corr = mean_corr;
plot_corr(eye(n_algs) == 1) = NaN;  % Set diagonal to NaN
plot_corr(triu(ones(n_algs), 1) == 1) = NaN;  % Set upper triangle to NaN

% Create heatmap
h = imagesc(plot_corr);
colormap('viridis');
caxis([min(plot_corr(:)) max(plot_corr(:))]);
c = colorbar;
set(c, 'TickLength', 0, 'TickLabels', []);  % Remove colorbar ticks and labels

% Customize appearance
axis square
set(gca, 'XTick', 1:n_algs, 'XTickLabel', [], ...
    'YTick', 1:n_algs, 'YTickLabel', [], ...
    'TickLength', [0 0], ...  % Remove tick marks
    'Box', 'off', ...         % Turn off box
    'XColor', 'k', ...        % Set x-axis color to black
    'YColor', 'k');           % Set y-axis color to black

% Add only bottom and left spines
ax = gca;
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;

% Create mask for upper triangle and diagonal
mask = tril(ones(n_algs), -1);  % Create lower triangular mask
set(h, 'AlphaData', mask);  % Apply mask to make upper triangle transparent

% Set figure background to transparent
set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');

% Save figure with transparency
exportgraphics(gcf, fullfile(main_dir, 'supplementary/algorithm_correlation.png'), 'Resolution', 1200, 'BackgroundColor', 'white');

%--------------------------------------------------------------------------%
%% Create histograms of edge distributions for VAR, FASK, and LiNGAM
%--------------------------------------------------------------------------%
close all

% Extract the group mean results for VAR, FASK, and LiNGAM
% First, ensure group_mean_results is calculated
if ~exist('group_mean_results', 'var')
    % Calculate group mean results if not already done
    group_mean_results = struct();
    algorithms = {'var', 'fask', 'lingam'};
    
    for k = 1:length(algorithms)
        alg = algorithms{k};
        all_subjects_data = arrayfun(@(x) x.(alg), ec_results, 'UniformOutput', false);
        all_subjects_data = cat(3, all_subjects_data{:}); % Concatenate along the third dimension
        all_subjects_data = all_subjects_data/max(all_subjects_data(:));
        
        % Calculate the mean across the third dimension (subjects)
        group_mean_results.(alg) = mean(all_subjects_data, 3);
    end
end

% Create a single figure for the histogrami
figure('Position', [100 100 400 400]);

% Define colors for each algorithm using the specified hex colors
colors = {[178/255 211/255 154/255],    % 99BC85 - Light green for VAR
          [251/255 228/255 214/255],    % FBE4D6 - Light peach for FASK
          [255/255 183/255 130/255]};   % FFC1B4 - Light coral for LiNGAM

% Define edge colors (darker versions of fill colors)
edge_colors = {[128/255 161/255 104/255],    % Darker green for VAR edges
               [201/255 178/255 164/255],    % Darker peach for FASK edges
               [205/255 133/255 80/255]};    % Darker coral for LiNGAM edges

alg_labels = {'VAR', 'FASK', 'LiNGAM'};
algorithms = {'var', 'fask', 'lingam'};

% Set common number of bins and range for all histograms
num_bins = 30;


% Hold on to plot multiple datasets
hold on;

% Create histograms
for i = 1:3
    % Get data and handle sparsity for FASK and LiNGAM
    data = group_mean_results.(algorithms{i})(:);
    if strcmp(algorithms{i}, 'fask') || strcmp(algorithms{i}, 'lingam')
        data = data(data > 0); % Only include non-zero values for sparse algorithms
    end
    
    % Normalize data to [0,1] range
    data = data / max(abs(data));
    
    % Create histogram with custom appearance and transparency
    h = histogram(data, 'NumBins', num_bins, 'Normalization', 'probability', ...
        'FaceColor', colors{i}, 'EdgeColor', edge_colors{i}, ...
        'LineWidth', 1.5, 'FaceAlpha', 0.7);
end

% Customize appearance
% box on;
% grid on;
set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', []);
set(gca, 'XTick', [-0.2 0.1 0.4 0.7 1]);
set(gca, 'YTick', [0.17 0.33 0.5]);
% set(gca,'YLim',[0 0.05])
% set(gca,'XLim',[0.4 1.05])

% Save figure
exportgraphics(gcf, fullfile(main_dir, 'supplementary/edge_distribution_histograms.png'), 'Resolution', 1200, 'BackgroundColor', 'white');

%--------------------------------------------------------------------------%
%% Investigate relationship between FASK predictions and FLN values
%--------------------------------------------------------------------------%
close all

% Initialize arrays to store median FLN values
fln_medians = zeros(num_subjects, 2); % Column 1 for mask=0, Column 2 for mask=1

% Loop through subjects
for subj = 1:num_subjects
    % Create binary mask from LiNGAM results
    lingam_mask = ec_results(subj).fask > 0;
    
    % Only consider non-zero FLN values
    nonzero_mask = FLN > 0;
    
    % Sample FLN using both masks
    fln_positive = FLN(lingam_mask & nonzero_mask);  % FLN values where LiNGAM > 0 and FLN > 0
    fln_zero = FLN(~lingam_mask & nonzero_mask); % FLN values where LiNGAM == 0 and FLN > 0
    
    % Store median values
    fln_medians(subj, 1) = mean(fln_zero);
    fln_medians(subj, 2) = mean(fln_positive);
end

% Create box plot
figure('Position', [100 100 300 400]);

% Create box plot with individual points
h = boxplot(fln_medians, {'FASK = 0', 'FASK > 0'}, 'Colors', 'k', 'Symbol', '');
hold on

% Customize box plot appearance
set(h, 'LineWidth', 1.5);
set(findobj(gca, 'type', 'line', 'Tag', 'Median'), 'Color', 'k', 'LineWidth', 2);

% Fill boxes with colors
colors = [0.8 0.8 0.8; 0.4 0.4 0.4]; % Light gray and dark gray
for i = 1:2
    patch(get(h(5,i), 'XData'), get(h(5,i), 'YData'), colors(i,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

% Add individual data points
for i = 1:2
    % Add jitter to x-coordinates
    x = i + (randn(num_subjects,1)*0.1);
    
    % Plot scattered points
    scatter(x, fln_medians(:,i), 40, [0.5 0.5 0.5], 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
end

% Customize plot appearance
box on
grid on

% Perform statistical test
[h, p] = ttest2(fln_medians(:,1), fln_medians(:,2));

% Add significance marker if applicable
if p < 0.05
    y_max = max(fln_medians(:)) +0.03;
    line([1 2], [y_max y_max], 'Color', 'k', 'LineWidth', 1.5);
    if p < 0.001
        text(1.5, y_max+0.01, '***', 'HorizontalAlignment', 'center', 'FontSize', 12);
    elseif p < 0.01
        text(1.5, y_max+0.01, '**', 'HorizontalAlignment', 'center', 'FontSize', 12);
    else
        text(1.5, y_max+0.01, '*', 'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    ylim([min(fln_medians(:))*0.9 y_max*1.1]);
end

set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 5))
% set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
% Save figure
% exportgraphics(gcf, fullfile(main_dir, 'supplementary/fask_fln_relationship.png'), 'Resolution', 1200);
