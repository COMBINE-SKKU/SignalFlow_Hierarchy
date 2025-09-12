%% Stage 3: Simulation validation and visualization
% This script performs Hopf model simulations to validate iEC performance
% against individual EC algorithms using goodness-of-fit metrics.
%
% IMPORTANT NOTE FOR FIGURE REPLICATION:
% Hopf simulations are inherently stochastic. Running this code will produce
% somewhat different results each time, though overall patterns remain consistent.
% 
% For exact figure replication, load pre-computed results using:
%   load([parcellation,'_gof_results.mat'])
% This file contains stored simulation results with the same random seed.
%
% Main workflow:
% 1. Data loading: Loads iEC matrices and preprocessed empirical data from previous stages
% 2. Simulation validation: Runs Hopf simulations for all algorithms + iEC variants
% 3. Performance evaluation: Calculates FC correlation, FCD KS, and goodness-of-fit metrics
% 4. Statistical analysis: Compares iEC vs best individual algorithm performance
% 5. Visualization: Creates comprehensive validation plots (GOF, FC, FCD, time-delay)
%
% Performance metric: Goodness-of-fit = FC correlation - FCD KS statistic
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

% Set flag for using stored results vs new simulations
use_stored_results = true;  % Set to false to run new simulations

%--------------------------------------------------------------------------%
%% Load empirical data and iEC matrices
%--------------------------------------------------------------------------%
% Load empirical FC and FCD for test set
fc_emp_test = importdata(fullfile(ec_dir, 'data', [parcellation,'_fc_fcd_td_test.mat'])).fc_emp_test;
fcd_emp_test = importdata(fullfile(ec_dir, 'data', [parcellation,'_fc_fcd_td_test.mat'])).fcd_emp_test;

% Load empirical EC algorithms for test set
ec_alg_groups_test = importdata(fullfile(ec_dir, 'data', [parcellation,'_ec_alg_groups_test.mat']));

% Load iEC matrices
iEC_Full = importdata(fullfile(ec_dir, 'data', [parcellation,'_resting_iEC_full.mat']));
iEC_vfl = importdata(fullfile(ec_dir, 'data', [parcellation,'_resting_iEC_vfl.mat']));

%--------------------------------------------------------------------------%
%% Load or compute simulation results
%--------------------------------------------------------------------------%

if use_stored_results
    % Load pre-computed results for exact figure replication
    fprintf('Loading pre-computed simulation results for figure replication...\n');
    try
        results_data = load(fullfile(ec_dir, 'data', [parcellation,'_gof_results.mat']));
        gof = results_data.gof;
        fc_corr = results_data.fc_corr;
        fcd_ks = results_data.fcd_ks;
        G_all = results_data.G_all;
        fprintf('Successfully loaded stored results.\n');
    catch
        warning('Could not load stored results. Running new simulations...');
        use_stored_results = false;
    end
end

if ~use_stored_results
    fprintf('Running new Hopf model simulations...\n');
    
    % Initialize parameters
    num_iter = 50;
    n_algs = length(algorithm_names);
    triu_ind_fc = find(triu(ones(size(fc_emp_test)), 1));
    fc_emp_triu = fc_emp_test(triu_ind_fc);
    
    % Initialize result matrices
    fc_corr = zeros(num_iter, n_algs+2);
    fcd_ks = zeros(num_iter, n_algs+2);
    gof = zeros(num_iter, n_algs+2);
    G_all = zeros(n_algs+2, 1);
    
    % Run simulations for individual algorithms
    for a = 1:n_algs
        fprintf('Processing algorithm %d/%d: %s\n', a, n_algs, algorithm_names{a});
        
        % Prepare EC matrix
        ec_temp = ec_alg_groups_test{a};
        ec_alg = ec_temp/max(ec_temp(:))*0.1;   
        
        % Set G range based on algorithm
        G_range = (a <= 2) * [1, 6] + (a > 2) * [3, 50];
        
        % Find optimal G and run simulations
        G_all(a) = findOptimalParams(ec_alg, f_diff, fc_emp_test, TR, Tmax, num_sims, triu_ind_fc, fcd_emp_test, true, [], false, G_range);
        
        % Pre-allocate temporary variables for parfor
        temp_fc_corr = zeros(num_iter, 1);
        temp_fcd_ks = zeros(num_iter, 1);
        temp_gof = zeros(num_iter, 1);
        G = G_all(a);
        
        parfor i = 1:num_iter
            % Run simulation
            bold_sim = run_simulation(ec_alg, G, f_diff, TR, Tmax, -0.01)';
            fc_sim = corr(bold_sim);
            fcd_sim = calculate_fcd(bold_sim, TR, ceil(60/TR), 1);
            
            % Calculate performance metrics
            temp_fc_corr(i) = corr(fc_sim(triu_ind_fc), fc_emp_triu);
            temp_fcd_ks(i) = max(abs(fcd_emp_test - fcd_sim));
            temp_gof(i) = temp_fc_corr(i) - temp_fcd_ks(i);
        end
        
        % Store results
        fc_corr(:,a) = temp_fc_corr;
        fcd_ks(:,a) = temp_fcd_ks;
        gof(:,a) = temp_gof;
    end
    
    % Run simulations for iEC variants
    iec_variants = {iEC_Full, iEC_vfl};
    iec_names = {'iEC_full', 'iEC_vfl'};
    
    for v = 1:2
        fprintf('Processing %s...\n', iec_names{v});
        
        % Find optimal G and run simulations
        G_all(n_algs+v) = findOptimalParams(iec_variants{v}, f_diff, fc_emp_test, TR, Tmax, num_sims, triu_ind_fc, fcd_emp_test, true, [], false, [0.7, 1.5]);
        
        % Pre-allocate temporary variables for parfor
        temp_fc_corr = zeros(num_iter, 1);
        temp_fcd_ks = zeros(num_iter, 1);
        temp_gof = zeros(num_iter, 1);
        G = G_all(n_algs+v);
        iec = iec_variants{v};
        
        parfor i = 1:num_iter
            % Run simulation
            bold_sim = run_simulation(iec, G, f_diff, TR, Tmax, -0.01)';
            fc_sim = corr(bold_sim);
            fcd_sim = calculate_fcd(bold_sim, TR, ceil(60/TR), 1);
            
            % Calculate performance metrics
            temp_fc_corr(i) = corr(fc_sim(triu_ind_fc), fc_emp_triu);
            temp_fcd_ks(i) = max(abs(fcd_emp_test - fcd_sim));
            temp_gof(i) = temp_fc_corr(i) - temp_fcd_ks(i);
        end
        
        % Store results
        fc_corr(:,n_algs+v) = temp_fc_corr;
        fcd_ks(:,n_algs+v) = temp_fcd_ks;
        gof(:,n_algs+v) = temp_gof;
    end
    
    % Save results for future replication
    save(fullfile(ec_dir, 'data', [parcellation,'_gof_results.mat']), 'gof', 'fc_corr', 'fcd_ks', 'G_all');
    fprintf('Simulation results saved for future replication.\n');
end

%--------------------------------------------------------------------------%
%% Create GOF boxplot
%--------------------------------------------------------------------------%

% Configure figure
figure('Position', [100 100 800 300]);
algorithm_labels = [algorithm_names, {'iEC_full', 'iEC_vfl'}];
num_iter = size(gof, 1);

% Create boxplot
h = boxplot(gof, 'Colors', 'k', 'Symbol', '', 'Labels', algorithm_labels);
hold on

% Customize boxplot appearance
set(h, 'LineWidth', 1.5);
set(findobj(gca, 'type', 'line', 'Tag', 'Median'), 'Color', 'k', 'LineWidth', 2);

% Define colors
colors = struct('var', [0.6 0.8 0.9], 'iec', [0.8 0.6 0.9], 'others', [0.9 0.8 0.6]);

% Color boxes
for i = 1:length(algorithm_labels)
    % Select color based on algorithm type
    if strcmp(algorithm_labels{i}, 'var') || strcmp(algorithm_labels{i}, 'rdcm')
        box_color = colors.var;
    elseif contains(algorithm_labels{i}, 'iEC')
        box_color = colors.iec;
    else
        box_color = colors.others;
    end
    
    % Apply color to box
    try
        box_data = get(h(5,i), 'XData');
        box_data_y = get(h(5,i), 'YData');
        if ~isempty(box_data) && ~isempty(box_data_y)
            patch(box_data, box_data_y, box_color, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        end
    catch
        warning('Could not color box for algorithm %s', algorithm_labels{i});
    end
end

% Add individual data points
for i = 1:length(algorithm_labels)
    x = i + (randn(num_iter,1)*0.05);
    scatter(x, gof(:,i), 10, [0.3 0.3 0.3], 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
end

% Statistical analysis
[~, best_idx] = max(median(gof(:,1:end-2), 1));
iec_data = gof(:,end);
best_alg_data = gof(:,best_idx);

% Perform statistical test
[p, ~, stats] = signrank(iec_data, best_alg_data, 'tail', 'right');

% Add significance markers if significant
if p < 0.05
    % Determine significance level
    if p < 0.001
        sig_str = '***';
    elseif p < 0.01
        sig_str = '**';
    else
        sig_str = '*';
    end
    
    % Get current y-axis limits
    ylim_current = ylim;
    y_range = ylim_current(2) - ylim_current(1);
    y_pos = ylim_current(2) + y_range * 0.02;
    x_pos = [best_idx, length(algorithm_labels)];
    
    % Draw significance bar
    plot(x_pos, [y_pos y_pos], 'k-', 'LineWidth', 1.5)
    plot([x_pos(1) x_pos(1)], [y_pos-y_range*0.01 y_pos], 'k-', 'LineWidth', 1.5)
    plot([x_pos(2) x_pos(2)], [y_pos-y_range*0.01 y_pos], 'k-', 'LineWidth', 1.5)
    
    % Add significance stars
    text(mean(x_pos), y_pos+y_range*0.02, sig_str, ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold')
    
    % Update y-axis limits
    ylim([ylim_current(1), y_pos+y_range*0.1])
end

% Finalize plot appearance
box on
grid on
set(gca, 'LineWidth', 1.5, 'TickDir', 'in')
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 5))
set(gca, 'YTickLabel', [], 'XTickLabel', [])

% Save figure
exportgraphics(gcf, fullfile(ec_dir, 'results', [parcellation, '_gof_comparison_boxplot.png']), 'Resolution', 1200);
savefig(gcf, fullfile(ec_dir, 'results', [parcellation, '_gof_comparison_boxplot.fig']));

%--------------------------------------------------------------------------%
%% Visualization of FC correlation and FCD comparison for iEC
%--------------------------------------------------------------------------%

% Run batch simulation for iEC_vfl
fprintf('Running batch simulation for FC/FCD visualization...\n');
[fc_sim_avg, fcd_sim_avg] = batch_simulation(iEC_vfl, 0.9, f_diff, TR, Tmax, num_sims);

% Prepare data for FC correlation scatter plot
triu_ind_fc = find(triu(ones(size(fc_emp_test)), 1));
fc_emp_triu = fc_emp_test(triu_ind_fc);
fc_sim_triu = fc_sim_avg(triu_ind_fc);

% Create FC correlation scatter plot
figure('Position', [100 100 350 350]);

% Plot scatter points
scatter(fc_emp_triu, fc_sim_triu, 10, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.7);
hold on;

% Add regression line
p = polyfit(fc_emp_triu, fc_sim_triu, 1);
x_range = linspace(min(fc_emp_triu), max(fc_emp_triu), 100);
y_fit = polyval(p, x_range);
plot(x_range, y_fit, 'r-', 'LineWidth', 2);

% Calculate correlation
r_value = corr(fc_emp_triu, fc_sim_triu);
fprintf('FC correlation: r = %.3f\n', r_value);

% Customize FC plot
grid on;
box on;
set(gca, 'XTick', linspace(min(fc_emp_triu), max(fc_emp_triu), 4));
set(gca, 'YTick', linspace(min(fc_sim_triu), max(fc_sim_triu), 4));
set(gca, 'XTickLabel', [], 'YTickLabel', []);
axis square;

% Save FC plot
% exportgraphics(gcf, fullfile(ec_dir, 'results', [parcellation, '_fc_correlation_scatter_iec.png']), 'Resolution', 1200);

% Create FCD CDF comparison plot
figure('Position', [500 100 350 350]);

% Plot CDFs
x_axis = linspace(-1, 1, length(fcd_emp_test));
plot(x_axis, fcd_emp_test, 'r-', 'LineWidth', 2, 'DisplayName', 'Empirical'); 
hold on;
plot(x_axis, fcd_sim_avg, 'b-', 'LineWidth', 2, 'DisplayName', 'Simulated');

% Find and plot maximum difference
[ks_dist, max_idx] = max(abs(fcd_emp_test - fcd_sim_avg));
plot([x_axis(max_idx) x_axis(max_idx)], [fcd_emp_test(max_idx), fcd_sim_avg(max_idx)], 'k--', 'LineWidth', 1.5);

fprintf('FCD KS distance: %.3f\n', ks_dist);

% Customize FCD plot
grid on;
box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 12, 'XTickLabel', [], 'YTickLabel', []);
axis square;

% Save FCD plot
% exportgraphics(gcf, fullfile(ec_dir, 'results', [parcellation, '_fcd_cdf_comparison_iec.png']), 'Resolution', 1200);

%--------------------------------------------------------------------------%
%% Time-delay analysis
%--------------------------------------------------------------------------%

% Load pre-computed empirical time delay data
fprintf('Running time-delay analysis...\n');
lag_data = importdata(fullfile(ec_dir, 'data', [parcellation,'_fc_fcd_td_test.mat'])).mean_lag_test;

% Find optimal G for iEC specifically for time delay analysis
iEC_vfl_norm = iEC_vfl/max(iEC_vfl(:))*0.1;
triu_ind_fc = find(triu(ones(size(fc_emp_test)), 1));
G = findOptimalParams(iEC_vfl_norm, f_diff, fc_emp_test, TR, Tmax, num_sims, triu_ind_fc, fcd_emp_test, false, lag_data, true, [0 10]);

% Run simulations with optimal G to calculate time delay matrices
td_matrix_iec = zeros(size(iEC_vfl_norm, 1), size(iEC_vfl_norm, 1), n_test);
fprintf('Computing time delay matrices (this may take a while)...\n');
parfor iter = 1:n_test
    ts_iec = run_simulation(iEC_vfl_norm, G, f_diff, TR, Tmax, -0.01)';
    td_matrix_iec(:,:,iter) = calculate_lag_threads(ts_iec, TR);
end
lag_sim = mean(td_matrix_iec, 3);

% Create scatter plot for mean time lag comparison
figure('Position', [100 100 350 350]);
lag_data_flipped = -lag_data;
lag_sim_flipped = -lag_sim; % flip the sign 
scatter(lag_data_flipped, lag_sim_flipped, 20, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.7);
hold on;

% Add linear regression line
p = polyfit(lag_data_flipped, lag_sim_flipped, 1);
x_range = linspace(min(lag_data_flipped), max(lag_data_flipped), 100);
y_fit = polyval(p, x_range);
plot(x_range, y_fit, 'r-', 'LineWidth', 2);

% Calculate correlation coefficient
r_value = corr(lag_data_flipped, lag_sim_flipped);
fprintf('Time delay correlation: r = %.3f\n', r_value);

% Customize appearance
grid on;
box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 12, 'XTickLabel', [], 'YTickLabel', []);
set(gca, 'XTick', linspace(min(lag_data_flipped), max(lag_data_flipped), 4));
set(gca, 'YTick', linspace(min(lag_sim_flipped), max(lag_sim_flipped), 4));
axis square;

% Save the figure
exportgraphics(gcf, fullfile(ec_dir, 'results', [parcellation, '_meanlag_correlation_scatter_iec.png']), 'Resolution', 1200);

% Create surface plots for time delay patterns
fprintf('Creating surface plots...\n');
surfaceplot(lag_sim_flipped, parcellation, 'both', 'BlRd')
exportgraphics(gcf, fullfile(ec_dir, 'results', [parcellation, 'iec_mean_lag.png']), 'Resolution', 1200);

surfaceplot(lag_data_flipped, parcellation, 'both', 'BlRd')
exportgraphics(gcf, fullfile(ec_dir, 'results', [parcellation, 'emp_mean_lag.png']), 'Resolution', 1200);
