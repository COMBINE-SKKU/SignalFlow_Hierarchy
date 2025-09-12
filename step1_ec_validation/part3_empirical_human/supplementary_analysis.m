%% Supplementary analysis: VAR algorithm parameter optimization
% This script performs supplementary analysis to determine optimal VAR lambda
% parameters for human empirical fMRI data using Hopf model simulations.
%
% Main workflow:
% 1. Parameter sweep: Tests VAR lambda values from 10 to 1000 on subset of HCP subjects
% 2. Statistical testing: Performs t-tests with FDR correction on VAR matrices
% 3. Hopf optimization: Finds optimal G parameter for each lambda using simulation
% 4. Goodness-of-fit evaluation: Computes FC correlation and FCD KS statistics
% 5. Edge analysis: Calculates proportion of significant edges across lambda values
% 6. Visualization: Creates plots showing GOF and edge proportions vs lambda
%
% Output: Optimal lambda parameter and comprehensive parameter sensitivity analysis
%
% Author: Younghyun Oh
% Date: 2025-03-04
%
%--------------------------------------------------------------------------%
% Clear workspace and add required paths
clear; close all; clc;
addpath(genpath('/local_raid1/01_software/toolboxes/toolboxes/cifti-matlab/'));

% Define paths
fmri_dir = '/combinelab/03_user/younghyun/03_data/fmri_data';
main_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy';
sub_dir = '/combinelab/02_data/01_HCP/HCP_S1200/';
ec_dir = fullfile(main_dir, 'step1_ec_validation', 'part3_empirical_human');
sublist = load(fullfile(sub_dir, 'sublist_HCP.txt'));

% Parameters
total_subjects = 440;
subset_percentage = 0.5;
subset_size = round(total_subjects * subset_percentage);
num_lambdas = 10;
lambda_range = linspace(10, 1000, num_lambdas);
TR = 0.72;
Tmax = 2400;
num_sims = 10;
alpha = 0.05;  % FDR significance level

% Set parcellation
parcellation = 'Schaefer100';
rois = strcmp(parcellation, 'MMP360') * 360 + strcmp(parcellation, 'Schaefer100') * 100;
f_diff = importdata(fullfile(main_dir,'/data/',[parcellation,'_peak_freq.mat']));

% Load group data
group_data = importdata(fullfile(ec_dir, 'data', [parcellation,'_fc_fcd_td_train.mat']));
FC_emp_group = group_data.fc_emp_train;
FCD_emp_group = group_data.fcd_emp_train;
FC_emp_triu = FC_emp_group(triu(true(rois,rois),1));

% Initialize storage
triu_ind_fc = triu(true(rois,rois),1);
group_VAR = zeros(rois, rois, num_lambdas);
group_corr_FC = zeros(1, num_lambdas);
group_fcd_ks = zeros(1, num_lambdas);
group_gof = zeros(1, num_lambdas);
G_values = zeros(num_lambdas, 1);
significant_edges = cell(num_lambdas, 1);
all_gof_values = cell(num_lambdas, 1);
all_fc_corrs = cell(num_lambdas,1);
all_fcd_ks = cell(num_lambdas,1);

% Select random subset of subjects
subset_indices = 1:subset_size;
subject_numbers = sublist(subset_indices);

%% Process each lambda value
for lambda_idx = 1:num_lambdas
    lambda = lambda_range(lambda_idx);
    fprintf('Processing lambda %.2f (%d/%d)\n', lambda, lambda_idx, num_lambdas);
    
    % Initialize VAR matrices storage
    var_matrices = zeros(rois, rois, subset_size);

    % Compute VAR for each subject in parallel
    parfor i = 1:subset_size
        subject_number = sublist(i);
        fprintf('  Processing subject %d/%d (ID: %d)\n', i, subset_size, subject_number);
        
        % Load and process subject data
        fmri_data = importdata(fullfile(fmri_dir, sprintf('sub%d/rfMRI_REST1_LRRL_%s_TianLv1_hp2000_clean.csv', subject_number, parcellation)));
        var_result = run_var(fmri_data, 600);
        var_matrices(:,:,i) = var_result(33:end,33:end);
    end

    % Statistical testing
    p_values = zeros(rois, rois);
    for i = 1:rois
        for j = 1:rois
            [~, p_values(i,j)] = ttest(squeeze(var_matrices(i,j,:)));
        end
    end

    % FDR correction
    [~, ~, ~, fdr_p_values] = fdr_bh(p_values(:));
    adjusted_p_values_matrix = reshape(fdr_p_values, rois, rois);

    % Compute and normalize group average VAR
    VAR_group = mean(var_matrices, 3);
    VAR_group = VAR_group/max(VAR_group(:))*0.1;

    % Store results
    significant_edges{lambda_idx} = adjusted_p_values_matrix <= alpha;
    group_VAR(:,:,lambda_idx) = VAR_group;
    
    % Find optimal G and run simulations
    G = findOptimalParams(VAR_group, f_diff, FC_emp_group, TR, Tmax, num_sims, triu_ind_fc, FCD_emp_group, true, [], false, [0.1 10]);
    G_values(lambda_idx) = G;
    
    % Parallel simulation runs
    fc_corrs = zeros(num_sims, 1);
    fcd_ks_values = zeros(num_sims, 1);
    gof_values = zeros(num_sims, 1);
    
    parfor iter = 1:num_sims
        bold_sim = run_simulation(VAR_group, G, f_diff, TR, Tmax, -0.01)';
        fc_sim = corr(bold_sim);
        fcd_sim = calculate_fcd(bold_sim, TR, ceil(60/TR), 1);
        
        fc_corrs(iter) = corr(fc_sim(triu_ind_fc), FC_emp_triu);
        fcd_ks_values(iter) = max(abs(FCD_emp_group - fcd_sim));
        gof_values(iter) = fc_corrs(iter) - fcd_ks_values(iter);
    end
    
    % Store metrics
    all_gof_values{lambda_idx} = gof_values;
    all_fc_corrs{lambda_idx} = fc_corrs;
    all_fcd_ks{lambda_idx} = fcd_ks_values;
    group_corr_FC(lambda_idx) = mean(fc_corrs);
    group_fcd_ks(lambda_idx) = mean(fcd_ks_values);
    group_gof(lambda_idx) = mean(gof_values);
end

% Find optimal lambda
[max_gof, best_lambda_idx] = max(group_gof);
optimal_lambda = lambda_range(best_lambda_idx);

% Calculate edge statistics
edge_counts = cellfun(@(x) sum(x(:)), significant_edges);
total_possible_edges = rois*(rois-1);
edge_proportions = edge_counts / total_possible_edges;

% Save results
results_dir = fullfile(main_dir, 'step1_ec_validation', 'part3_empirical_human', 'results');
save(fullfile(results_dir, sprintf('optimal_lambda_results_%s.mat', parcellation)), ...
    'optimal_lambda', 'lambda_range', 'group_gof', 'all_gof_values', 'group_fcd_ks', 'group_corr_FC', ...
    'edge_counts', 'edge_proportions', 'G_values');

%% Visualization
figure;
% GOF plot with error bands
subplot(2,1,1);
hold on;
gof_sem = cellfun(@(x) std(x)/sqrt(num_sims), all_gof_values);
upper = (group_gof + gof_sem');
lower = (group_gof - gof_sem');
fill([lambda_range, fliplr(lambda_range)], ...
     [upper, fliplr(lower)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(lambda_range, group_gof, 'b-', 'LineWidth', 2);
grid on;

% Edge proportions plot
subplot(2,1,2);
plot(lambda_range, edge_proportions, 'g-o', 'LineWidth', 2);
grid on; box off;

% Set consistent ticks
for i = 1:2
    subplot(2,1,i);
    set(gca, 'XTick', linspace(min(lambda_range), max(lambda_range), 4));
    set(gca, 'YTick', linspace(min(get(gca, 'YLim')), max(get(gca, 'YLim')), 4));
    % set(gca, 'YTickLabel', [], 'XTickLabel', [])
    set(gca, 'Xlim', [0 780])
end

% Save figure
set(gcf, 'Position', [100, 100, 400, 350]);
savefig(fullfile(results_dir, sprintf('metrics_vs_lambda_%s_subset_parallel_fdr.fig', parcellation)));
exportgraphics(gcf, fullfile(results_dir, sprintf('metrics_vs_lambda_%s_subset_parallel_fdr.png', parcellation)), 'Resolution', 1200);
