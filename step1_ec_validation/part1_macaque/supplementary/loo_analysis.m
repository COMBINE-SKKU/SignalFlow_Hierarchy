%% Stage 4: Leave-One-Out (LOO) Analysis
% This script performs post-hoc analysis comparing three models:
% 1. VFL_full: VAR + FASK + LiNGAM (original top-3 model)
% 2. VFL_noFASK: VAR + LiNGAM only
% 3. VFL_noLiNGAM: VAR + FASK only
%
% Main workflow:
% 1. Recompute optimal weights for LOO models using training data (subjects 1-8)
% 2. Apply models to test subjects (9-19) with proper normalization
% 3. Compute correlations with FLN ground truth
% 4. Statistical comparison using two-tailed non-parametric tests
% 5. Visualization with boxplots and significance markers
%
% Author: Younghyun Oh
% Date: 2025-11-03
%
%--------------------------------------------------------------------------%
% Clear workspace
clear; close all; clc

% Set up directories
base_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy';
main_dir = fullfile(base_dir, 'step1_ec_validation/part1_macaque');

% Load FLN ground truth
FLN_dir = fullfile(main_dir, 'data/FLN40.mat');
FLN = importdata(FLN_dir);
FLN = 1.2*FLN.^0.3; % Log-linear transformation

% Load EC results
result_dir = fullfile(main_dir, 'results');
ec_results = load(fullfile(result_dir, 'ec_results.mat')).ec_results;

% Load existing VFL betas for comparison
EC_betas_existing = load(fullfile(result_dir, 'EC_betas.mat')).EC_betas_vfl;

fprintf('Loaded data successfully. Starting LOO analysis...\n');

%--------------------------------------------------------------------------%
%% Recompute weights for LOO models using training data (subjects 1-8)
%--------------------------------------------------------------------------%

% Initialize storage for new beta weights
EC_betas_noFASK = zeros(8, 2);   % VAR + LiNGAM
EC_betas_noLiNGAM = zeros(8, 2); % VAR + FASK

fprintf('Recomputing weights for LOO models...\n');

% Training on subjects 1-8
for sub = 1:8
    fprintf('Processing training subject %d/8...\n', sub);

    % Extract and normalize matrices
    var_norm = ec_results(sub).var / max(ec_results(sub).var(:));
    fask_norm = ec_results(sub).fask / max(ec_results(sub).fask(:));
    lingam_norm = ec_results(sub).lingam / max(ec_results(sub).lingam(:));

    % LOO Model 1: VAR + LiNGAM (no FASK)
    loo_matrices_noFASK = {var_norm, lingam_norm};
    Betas_noFASK = findBetas_FLN(loo_matrices_noFASK, FLN);
    EC_betas_noFASK(sub, :) = Betas_noFASK;

    % LOO Model 2: VAR + FASK (no LiNGAM)
    loo_matrices_noLiNGAM = {var_norm, fask_norm};
    Betas_noLiNGAM = findBetas_FLN(loo_matrices_noLiNGAM, FLN);
    EC_betas_noLiNGAM(sub, :) = Betas_noLiNGAM;
end

% Calculate median weights for test application
median_betas_full = median(EC_betas_existing);     % VAR, FASK, LiNGAM
median_betas_noFASK = median(EC_betas_noFASK);     % VAR, LiNGAM
median_betas_noLiNGAM = median(EC_betas_noLiNGAM); % VAR, FASK

fprintf('\nMedian weights computed:\n');
fprintf('Full model (VAR, FASK, LiNGAM): [%.4f, %.4f, %.4f]\n', median_betas_full);
fprintf('No FASK (VAR, LiNGAM): [%.4f, %.4f]\n', median_betas_noFASK);
fprintf('No LiNGAM (VAR, FASK): [%.4f, %.4f]\n', median_betas_noLiNGAM);

%--------------------------------------------------------------------------%
%% Apply models to test subjects (9-19) and compute correlations
%--------------------------------------------------------------------------%

% Initialize correlation storage
correlations = zeros(11, 3); % 11 test subjects x 3 models
test_subjects = 9:19;

fprintf('\nApplying models to test subjects...\n');

for i = 1:length(test_subjects)
    sub = test_subjects(i);
    fprintf('Processing test subject %d (%d/11)...\n', sub, i);

    % Extract and normalize matrices
    var_norm = ec_results(sub).var / max(ec_results(sub).var(:));
    fask_norm = ec_results(sub).fask / max(ec_results(sub).fask(:));
    lingam_norm = ec_results(sub).lingam / max(ec_results(sub).lingam(:));

    % Model 1: Full VFL (VAR + FASK + LiNGAM)
    vfl_full = (var_norm * median_betas_full(1) + ...
                fask_norm * median_betas_full(2) + ...
                lingam_norm * median_betas_full(3)) / 3;
    correlations(i, 1) = corr(FLN(:), vfl_full(:));

    % Model 2: VFL without FASK (VAR + LiNGAM)
    vfl_noFASK = (var_norm * median_betas_noFASK(1) + ...
                  lingam_norm * median_betas_noFASK(2)) / 2;
    correlations(i, 2) = corr(FLN(:), vfl_noFASK(:));

    % Model 3: VFL without LiNGAM (VAR + FASK)
    vfl_noLiNGAM = (var_norm * median_betas_noLiNGAM(1) + ...
                    fask_norm * median_betas_noLiNGAM(2)) / 2;
    correlations(i, 3) = corr(FLN(:), vfl_noLiNGAM(:));
end

% Display correlation statistics
fprintf('\nCorrelation Statistics:\n');
fprintf('======================\n');
model_names = {'VFL_full', 'VFL_noFASK', 'VFL_noLiNGAM'};
for i = 1:3
    fprintf('%s: Mean=%.4f, Std=%.4f, Median=%.4f\n', ...
        model_names{i}, mean(correlations(:,i)), std(correlations(:,i)), median(correlations(:,i)));
end

%--------------------------------------------------------------------------%
%% Statistical Testing
%--------------------------------------------------------------------------%

fprintf('\nPerforming statistical tests...\n');

% Pairwise comparisons using Wilcoxon signed-rank test (two-tailed)
% Comparison 1: Full vs noFASK
[p_full_vs_noFASK, h1] = signrank(correlations(:,1), correlations(:,2), 'tail', 'both');

% Comparison 2: Full vs noLiNGAM
[p_full_vs_noLiNGAM, h2] = signrank(correlations(:,1), correlations(:,3), 'tail', 'both');

% Comparison 3: noFASK vs noLiNGAM
[p_noFASK_vs_noLiNGAM, h3] = signrank(correlations(:,2), correlations(:,3), 'tail', 'both');

% Apply multiple comparisons correction (FDR)
p_values_uncorrected = [p_full_vs_noFASK, p_full_vs_noLiNGAM, p_noFASK_vs_noLiNGAM];
alpha = 0.05;
[h_corrected, crit_p, adj_ci_cvrg, p_values_corrected] = fdr_bh(p_values_uncorrected, alpha, 'pdep', 'yes');

fprintf('\nStatistical Test Results:\n');
fprintf('========================\n');
fprintf('Full vs noFASK: p = %.4f, p_corrected = %.4f (h = %d)\n', p_full_vs_noFASK, p_values_corrected(1), h_corrected(1));
fprintf('Full vs noLiNGAM: p = %.4f, p_corrected = %.4f (h = %d)\n', p_full_vs_noLiNGAM, p_values_corrected(2), h_corrected(2));
fprintf('noFASK vs noLiNGAM: p = %.4f, p_corrected = %.4f (h = %d)\n', p_noFASK_vs_noLiNGAM, p_values_corrected(3), h_corrected(3));

%--------------------------------------------------------------------------%
%% Visualization
%--------------------------------------------------------------------------%

close all

% Create figure
figure('Position', [100 100 400 350]);

% Define model labels
model_labels = {'Full', 'No FASK', 'No LiNGAM'};

% Create box plot
h = boxplot(correlations, 'Colors', 'k', 'Symbol', '', 'Labels', model_labels);
hold on

% Customize box plot appearance
set(h, 'LineWidth', 1.5);
set(findobj(gca, 'type', 'line', 'Tag', 'Median'), 'Color', 'k', 'LineWidth', 2);

% Define colors for boxes
colors = [0.8 0.6 0.9;     % Pastel purple for Full
          0.9 0.8 0.6;     % Pastel orange for No FASK
          0.6 0.8 0.9];    % Pastel blue for No LiNGAM

% Fill boxes with colors
for i = 1:3
    patch(get(h(5,i), 'XData'), get(h(5,i), 'YData'), colors(i,:), ...
        'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

% Add individual data points with jitter
for i = 1:3
    curr_data = correlations(:,i);
    x = i + (randn(length(curr_data),1)*0.1);
    scatter(x, curr_data, 40, [0.3 0.3 0.3], 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
end

% Add significance brackets and markers
y_max = max(correlations(:));
y_pad = 0.005;        % Small vertical padding above box for first line
y_step = 0.01;       % Small step between lines
marker_pad = 0.003;   % Padding for marker text above line

% Test significance levels and add markers (using corrected p-values)
comparisons = {[1 2], [1 3], [2 3]};

for comp_idx = 1:3
    p_val = p_values_corrected(comp_idx);
    pos = comparisons{comp_idx};

    if p_val < 0.05
        % Draw horizontal line closer to boxes
        y_line = y_max + y_pad + (comp_idx-1)*y_step;
        line(pos, [y_line y_line], 'Color', 'k', 'LineWidth', 1.5);

        % Add significance marker just above the line
        x_text = mean(pos);
        if p_val < 0.001
            marker = '***';
        elseif p_val < 0.01
            marker = '**';
        else
            marker = '*';
        end
        text(x_text, y_line + marker_pad, marker, 'HorizontalAlignment', 'center', 'FontSize', 12);
    end
end

% Customize plot appearance
box on
grid on
set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
ylabel_pos = y_max + y_pad + 2*y_step + marker_pad + 0.02; % minimal space above highest bracket
ylim([min(correlations(:))*0.95, ylabel_pos]);
set(gca, 'YTick', linspace(min(correlations(:)), ylabel_pos, 5))

% Save figure
exportgraphics(gcf, fullfile(main_dir, 'supplementary/loo_analysis_boxplot.png'), 'Resolution', 1200);

%--------------------------------------------------------------------------%
%% Save results
%--------------------------------------------------------------------------%

% Save analysis results
save(fullfile(result_dir, 'loo_analysis_results.mat'), ...
    'correlations', 'EC_betas_noFASK', 'EC_betas_noLiNGAM', ...
    'median_betas_full', 'median_betas_noFASK', 'median_betas_noLiNGAM', ...
    'p_full_vs_noFASK', 'p_full_vs_noLiNGAM', 'p_noFASK_vs_noLiNGAM', ...
    'p_values_corrected', 'h_corrected', 'crit_p');

fprintf('\nLOO analysis completed successfully!\n');
fprintf('Results saved to: %s\n', fullfile(result_dir, 'loo_analysis_results.mat'));
fprintf('Figure saved to: %s\n', fullfile(main_dir, 'results/loo_analysis_boxplot.png'));