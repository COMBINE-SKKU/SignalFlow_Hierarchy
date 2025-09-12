%% Supplementary analysis: Signal characteristics and algorithm variability
% This script performs supplementary analyses investigating signal characteristics
% and algorithm parameter sensitivity for the macaque EC validation study.
%
% Main analyses:
% 1. Non-Gaussian structure: Compares skewness of empirical vs random BOLD signals
% 2. FLN bidirectionality: Measures proportion of bidirectional connections in ground truth
% 3. VAR lambda sensitivity: Tests correlation performance across different lambda values
% 4. Visualization: Creates distribution plots and parameter sensitivity curves
%
% Author: Younghyun Oh
% Date: 2025-02-26
%
%--------------------------------------------------------------------------%
%% First, check the non-Gaussian structure of the BOLD signals (Macaque)
%--------------------------------------------------------------------------%
% Clear workspace
clear; close all; clc

% Set constant
species = 'macaque';
base_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy';
main_dir = fullfile(base_dir, 'step1_ec_validation/part1_macaque');
switch species
    case 'macaque'
        % Load BOLD
        BOLD_dir = fullfile(base_dir, 'data/BOLD(individual)_macaque.mat');
        BOLDs = importdata(BOLD_dir);
    case 'human'
        % Load BOLD from random 100 subjects
        BOLDs = cell(100, 1);
        fmri_dir = '/combinelab/03_user/younghyun/03_data/fmri_data';
        sublist = load('/combinelab/02_data/01_HCP/HCP_S1200/sublist_HCP.txt');
        sublist = sublist(1:440);
        sub_index = randperm(length(sublist), 100);
        for subject = 1:100
            subject_number = sublist(sub_index(subject));
            
            sub_dir = fullfile(fmri_dir, sprintf('sub%d', subject_number));
            fmri_file = fullfile(sub_dir, sprintf('rfMRI_REST1_LR_%s_TianLv1_hp2000_clean.csv', 'MMP360'));
            BOLDs{subject} = importdata(fmri_file)';
        end
end

%% Compute skewness for each region for each subject
num_subjects = length(BOLDs);
num_regions = size(BOLDs{1}, 2);
skewness_empirical = zeros(num_subjects, num_regions);

for subj = 1:num_subjects
    for region = 1:num_regions
        % Take absolute value of skewness
        skewness_empirical(subj, region) = skewness(BOLDs{subj}(:,region));
    end
end

% Compute skewness for random Gaussian variables
skewness_random = zeros(num_subjects, num_regions);
for subj = 1:num_subjects
    for region = 1:num_regions
        random_data = randn(size(BOLDs{subj}(:,region)));
        % Take absolute value of skewness
        skewness_random(subj, region) = skewness(random_data);
    end
end

%% Visualize the distribution of skewness_empirical and skewness_random
close all
figure('Position', [100, 100, 550, 450]);
hold on;
histogram(skewness_empirical(:), 'Normalization', 'probability', 'EdgeColor', 'b', 'FaceColor', [0.2, 0.6, 0.8], 'LineWidth', 1.5, 'NumBins', 20);
histogram(skewness_random(:), 'Normalization', 'probability', 'EdgeColor', 'r', 'FaceColor', [0.8, 0.2, 0.2], 'LineWidth', 1.5, 'NumBins', 20);
hold off;
grid on;
box on;
set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 4))
set(gca, 'XTick', linspace(min(get(gca, 'XTick')), max(get(gca, 'XTick')), 5))
exportgraphics(gcf, fullfile(main_dir, 'supplementary/skewness_distribution(macaque).png'), 'Resolution', 1200);

% Aggregate skewness values over all subjects and regions
mean_skewness_empirical = mean(skewness_empirical(:));
mean_skewness_random = mean(skewness_random(:));

% Calculate the standard error
std_error_empirical = std(skewness_empirical(:)) / sqrt(num_subjects * num_regions);
std_error_random = std(skewness_random(:)) / sqrt(num_subjects * num_regions);

% Two-tailed Z-test for aggregated skewness
z_score = (mean_skewness_empirical - mean_skewness_random) / sqrt(std_error_empirical^2 + std_error_random^2);
p_value = 2 * (1 - normcdf(abs(z_score)));

%--------------------------------------------------------------------------%
%% Second, Measure the proportion of bi-directional connections
%--------------------------------------------------------------------------%
% Load FLN
FLN_dir = fullfile(main_dir, 'data/FLN40.mat');
FLN = importdata(FLN_dir);
FLN = 1.2*FLN.^0.3; % Log-linear transformation

% Calculate median of non-zero FLN values
nonzero_FLN = FLN(FLN > 0);
median_FLN = mean(nonzero_FLN);

% Calculate proportion of edges below median
below_median_count = sum(FLN > 0 & FLN < median_FLN, 'all');
total_connections = sum(FLN > 0, 'all');
proportion_below_median = below_median_count / total_connections;

fprintf('Median FLN value: %.4f\n', median_FLN);
fprintf('Proportion of edges below median: %.4f\n', proportion_below_median);

num_regions = size(FLN, 1);
bidirectional_count = 0;
total_possible_pairs = 0;

for i = 1:num_regions
    for j = i+1:num_regions  % Only upper triangle to avoid counting pairs twice
        if i == j
            continue;
        end
        % Count as bidirectional if both directions exist
        if FLN(i, j) > 0 && FLN(j, i) > 0
            bidirectional_count = bidirectional_count + 1;
        end
        % Count total number of pairs that have at least one connection
        if FLN(i, j) > 0 || FLN(j, i) > 0
            total_possible_pairs = total_possible_pairs + 1;
        end
    end
end

proportion_bidirectional = bidirectional_count / total_possible_pairs;
fprintf('Number of bidirectional connections: %d\n', bidirectional_count);
fprintf('Total number of connected pairs: %d\n', total_possible_pairs);
fprintf('Proportion of bi-directional connections: %.3f\n', proportion_bidirectional);

%--------------------------------------------------------------------------%
%% Third, Analyze the impact of VAR's lambda
%--------------------------------------------------------------------------%
num_lambdas = 100;
lambda_range = linspace(1e+1,1e+10,100);
corr_FLN = zeros(num_subjects, num_lambdas);

% Loop over subjects
for subj = 1:num_subjects
    BOLD = BOLDs{subj};

    % Inner loop over lambda (from 0 to 1e+10 with 100 steps)
    for lambda_idx = 1:num_lambdas
        lambda = lambda_range(lambda_idx);

        % Run VAR
        VAR = run_var(BOLD', lambda);
        VAR = VAR/max(VAR(:)); % normalize

        % Compute correlation with FLN
        corr_FLN(subj, lambda_idx) = corr(FLN(:), VAR(:));
    end
end

% Save results
% save('supplementary/var_lambda_correlation.mat', 'corr_FLN');

%% Plot the results
close all

% First plot: Lambda vs FLN correlation
figure('Position', [100 100 400 400]);

sem_FLN = std(corr_FLN, [], 1) / sqrt(size(corr_FLN, 1));

% Plot shaded error bars for FLN correlation
fill([lambda_range fliplr(lambda_range)], ...
    [mean(corr_FLN, 1) + sem_FLN fliplr(mean(corr_FLN, 1) - sem_FLN)], ...
    [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;

% Plot mean line for FLN correlation
plot(lambda_range, mean(corr_FLN, 1), 'Color', [0 0 0.8], ...
    'LineWidth', 2);

% Customize plot
set(gca, 'XScale', 'log');
box on;
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 5))
set(gca, 'XTick', [1e+3, 2e+3, 3e+3,4e+3])
set(gca, 'LineWidth', 1.5)
set(gca,'YTickLabel',[], 'XTickLabel',[])
set(gca,'title',[])

% Save figure
exportgraphics(gcf, fullfile(main_dir, 'supplementary/corr_FC_vs_lambda_Schaefer100.png'), 'Resolution', 1200);