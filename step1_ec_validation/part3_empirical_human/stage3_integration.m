%% Stage 2: iEC integration and weight optimization
% This script performs beta optimization to find optimal integration weights
% for combining multiple EC algorithms into integrated effective connectivity (iEC).
%
% Main workflow:
% 1. Data loading: Loads preprocessed group-averaged EC matrices from stage1
% 2. Beta optimization: Uses findBetas_new with Hopf simulations on training data
% 3. Algorithm selection: Supports both full (8 algorithms) and reduced (VAR+FASK) combinations
% 4. Statistical analysis: Performs significance testing on integration weights with FDR correction
% 5. Visualization: Creates box plots showing algorithm contribution distributions
% 6. iEC construction: Builds integrated connectivity matrices using optimal weights
%
% Output: Integration weights, iEC matrices, and weight distribution visualizations
%
% Author: Younghyun Oh
% Date: 2025-06-10
%--------------------------------------------------------------------------%
%% Set up
%--------------------------------------------------------------------------%
% Clear workspace
clear; close all; clc;

% Set constants
parcellation = 'MMP360';
n_subjects = 440;
n_train = n_subjects/2;
n_test = n_subjects - n_train;
Tmax = 2400; TR = 0.72; num_sims = 10;
rois = strcmp(parcellation, 'MMP360') * 360 + strcmp(parcellation, 'Schaefer100') * 100;

% Set directories
main_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy';
ec_dir = fullfile(main_dir, 'step1_ec_validation', 'part3_empirical_human');
f_diff = importdata(fullfile(main_dir,'/data/',[parcellation,'_peak_freq.mat']));
algorithm_names = {'rdcm','var','fask','ccd','boss','lingam','grasp','patel'};

%--------------------------------------------------------------------------%
%% Load empirical data
%--------------------------------------------------------------------------%
% Load empirical FC and FCD
fc_emp_train = importdata(fullfile(ec_dir, 'data', [parcellation,'_fc_fcd_td_train.mat'])).fc_emp_train;
fcd_emp_train = importdata(fullfile(ec_dir, 'data', [parcellation,'_fc_fcd_td_train.mat'])).fcd_emp_train;

% Load empirical EC algorithms
ec_alg_groups_train = importdata(fullfile(ec_dir, 'data', [parcellation,'_ec_alg_groups_train.mat']));
ec_alg_groups_test = importdata(fullfile(ec_dir, 'data', [parcellation,'_ec_alg_groups_test.mat']));

%--------------------------------------------------------------------------%
%% Find optimal beta values
%--------------------------------------------------------------------------%

% Flags for the findBetas function
use_all_algorithms = true;  % Set to true to use all algorithms
use_fcd = true;

% Get upper triangular indices for FC correlation calculation
triu_ind_fc = find(triu(ones(size(fc_emp_train)), 1));  % Upper triangular indices

% Define algorithms to use
if use_all_algorithms
    selected_algorithms = algorithm_names;
else
    selected_algorithms = {'var', 'fask'};  % Default selected algorithms
end

% Set up input matrices
n_algs = length(selected_algorithms);
input_matrices = cell(n_algs, 1);

% Extract EC matrices for selected algorithms
for a = 1:n_algs
    alg_idx = find(strcmp(algorithm_names, selected_algorithms{a}));
    if isempty(alg_idx)
        error('Algorithm %s not found in algorithm_names', selected_algorithms{a});
    end
    
    ec_temp = ec_alg_groups_train{alg_idx};
    
    % Normalize and store
    input_matrices{a} = ec_temp/max(abs(ec_temp(:)))*0.1; 
end

% Find optimal betas
if use_all_algorithms
    num_iter = 60;
    betas = zeros(num_iter, n_algs);
    gofs = zeros(num_iter, 1);
    G_opt = zeros(num_iter, 1);
    w_opt = zeros(num_iter, n_algs);
    for iter = 1:num_iter
        [betas(iter,:), G_opt(iter), w_opt(iter,:), gofs(iter)] = findBetas_new(input_matrices, f_diff, fc_emp_train, TR, Tmax, num_sims, triu_ind_fc, fcd_emp_train);
    end

    % Statistical analysis on normalized weights
    group2_w = w_opt(11:end,:);
    group2_w_norm = group2_w ./ sum(group2_w,2);
    group2_centered = group2_w_norm - 1/8;
    pvals = zeros(1, size(group2_centered,2));
    for i = 1:size(group2_centered,2)
        [~,pvals(i)] = ttest(group2_centered(:,i), 0, 'Tail', 'right');
    end

    % Multiple comparison correction
    [~,~,~,fdr_pvals] = fdr_bh(pvals);
else
    [betas, G_opt, w_opt, gofs] = findBetas_new(input_matrices, f_diff, fc_emp_train, TR, Tmax, num_sims, triu_ind_fc, fcd_emp_train);
    fdr_pvals = [];
end

% Save results
switch use_all_algorithms
    case true
        save(fullfile(ec_dir, 'data', [parcellation,'_betas_group.mat']), 'betas', 'G_opt', 'w_opt', 'gofs', 'selected_algorithms');
    case false
        save(fullfile(ec_dir, 'data', [parcellation,'_betas_group_vfl.mat']), 'betas', 'G_opt', 'w_opt', 'gofs', 'selected_algorithms');
end

%--------------------------------------------------------------------------%
%% Create box plot for w_opt values
%--------------------------------------------------------------------------%
close all

% Configure figure size
figure('Position', [100 100 650 300]);

% Create box plot with individual points
h = boxplot(w_opt, 'Colors', 'k', 'Symbol', '', 'Labels', selected_algorithms);
hold on

% Customize box plot appearance
set(h, 'LineWidth', 1.5);
set(findobj(gca, 'type', 'line', 'Tag', 'Median'), 'Color', 'k', 'LineWidth', 2);

% Define colors for each group
color_var = [0.6 0.8 0.9];    % Pastel blue for var
color_others = [0.9 0.8 0.6]; % Pastel yellow/orange for others

% Fill the boxes with different colors based on groups
for i = 1:length(selected_algorithms)
    if strcmp(selected_algorithms{i}, 'var') || strcmp(selected_algorithms{i}, 'rdcm')
        box_color = color_var;
    else
        box_color = color_others;
    end
    
    % Try to get box data
    try
        box_data = get(h(5,i), 'XData');
        box_data_y = get(h(5,i), 'YData');
        if ~isempty(box_data) && ~isempty(box_data_y)
            patch(box_data, box_data_y, box_color, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        end
    catch
        warning('Could not color box for algorithm %s', selected_algorithms{i});
    end
end

% Add individual data points
for i = 1:length(selected_algorithms)
    % Get data points for current algorithm
    curr_data = w_opt(:,i);
    
    % Add jitter to x-coordinates
    x = i + (randn(length(curr_data),1)*0.05);
    
    % Plot scattered points
    scatter(x, curr_data, 10, [0.3 0.3 0.3], 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
end

% Add significance markers for full algorithm case
if use_all_algorithms && ~isempty(fdr_pvals)
    % Get current y-axis limits for consistent marker placement
    ylim_current = ylim;
    y_range = ylim_current(2) - ylim_current(1);
    y_pos = ylim_current(2) + y_range * 0.02;
    
    % Add significance markers based on FDR-corrected p-values
    for i = 1:length(fdr_pvals)
        if fdr_pvals(i) < 0.05
            sig_str = '*';
            % Add significance stars at consistent y-position
            text(i, y_pos, sig_str, ...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold')
        end
    end
    
    % Update y-axis limits to show significance markers
    ylim([ylim_current(1), y_pos+y_range*0.1])
end

% Customize plot appearance
box on
grid on
set(gca, 'LineWidth', 1.5, 'TickDir', 'in')
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 5))
set(gca, 'YTickLabel', [], 'XTickLabel', [])

% Save the figure
suffix = use_all_algorithms * '_full' + ~use_all_algorithms * '_vfl';
% exportgraphics(gcf, fullfile(ec_dir, 'results', [parcellation, '_wopt_boxplot', suffix, '.png']), 'Resolution', 1200);
% savefig(gcf, fullfile(ec_dir, 'results', [parcellation, '_wopt_boxplot', suffix, '.fig']));

%--------------------------------------------------------------------------%
%% Construct iEC matrices
%--------------------------------------------------------------------------%

% Load the optimization results for both full and VFL versions
data_all = importdata(fullfile(ec_dir, 'data', [parcellation,'_betas_group.mat']));

% For VFL version, run separate optimization if not done above
if ~use_all_algorithms
    data_vfl = importdata(fullfile(ec_dir, 'data', [parcellation,'_betas_group_vfl.mat']));
else
    % If we ran full algorithm optimization, also create VFL version
    vfl_algorithms = {'var', 'fask'};
    n_vfl = length(vfl_algorithms);
    input_matrices_vfl = cell(n_vfl, 1);
    
    for a = 1:n_vfl
        alg_idx = find(strcmp(algorithm_names, vfl_algorithms{a}));
        ec_temp = ec_alg_groups_train{alg_idx};
        input_matrices_vfl{a} = ec_temp/max(abs(ec_temp(:)))*0.1;
    end
    
    [betas_vfl, G_opt_vfl, w_opt_vfl, gofs_vfl] = findBetas_new(input_matrices_vfl, f_diff, fc_emp_train, TR, Tmax, num_sims, triu_ind_fc, fcd_emp_train);
    data_vfl = struct('betas', betas_vfl, 'G_opt', G_opt_vfl, 'w_opt', w_opt_vfl, 'gofs', gofs_vfl, 'selected_algorithms', {vfl_algorithms});
    save(fullfile(ec_dir, 'data', [parcellation,'_betas_group_vfl.mat']), '-struct', 'data_vfl');
end

% Find the index of the maximum gof for each algorithm set
[~, opt_idx_all] = max(data_all.gofs(:,1));
opt_beta_all = data_all.betas(opt_idx_all,:);
opt_beta_vfl = data_vfl.betas;

% Construct iEC_full using all algorithms
n_algs_full = length(data_all.selected_algorithms);
iEC_Full = zeros(rois,rois);
for a = 1:n_algs_full
    alg_idx = find(strcmp(algorithm_names, data_all.selected_algorithms{a}));
    ec_temp = ec_alg_groups_test{alg_idx};
    ec_temp = ec_temp/max(abs(ec_temp(:)))*0.1;
    iEC_Full = iEC_Full + opt_beta_all(a)*ec_temp;
end
save(fullfile(ec_dir, 'data', [parcellation,'_resting_iEC_full.mat']), 'iEC_Full');

% Construct iEC_vfl using selected algorithms
n_algs_vfl = length(data_vfl.selected_algorithms);
iEC_vfl = zeros(rois,rois);
for a = 1:n_algs_vfl
    alg_idx = find(strcmp(algorithm_names, data_vfl.selected_algorithms{a}));
    ec_temp = ec_alg_groups_test{alg_idx};
    ec_temp = ec_temp/max(abs(ec_temp(:)))*0.1;
    iEC_vfl = iEC_vfl + opt_beta_vfl(a)*ec_temp;
end
save(fullfile(ec_dir, 'data', [parcellation,'_resting_iEC_vfl.mat']), 'iEC_vfl');
