%% stage0_base_model.m
% This script compares empirical and simulated BOLD dynamics using FCD analysis
%
% The script:
% 1. Loads empirical fMRI data and FLN
% 2. Simulates BOLD signals using Hopf model
% 3. Computes and compares FCD distributions between empirical and simulated data
%
% Measures:
% - FCD (Functional Connectivity Dynamics)
% - KS statistic for distribution comparison
%
% Dependencies:
% - calculate_fcd.m: Function to compute FCD matrices and vectors
% - run_simulation.m: Function to run Hopf model simulation
%
% Author: Younghyun Oh
% Date: 2025-01-14
% Version: 1.1

%% Initialize
clear; close all; clc

% Load FLN
base_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy';
main_dir = fullfile(base_dir, 'step1_ec_validation');
FLN_dir = fullfile(main_dir, 'part1_macaque/data/FLN40.mat');
FLN = importdata(FLN_dir);
FLN = FLN/mean(FLN(:))*0.02;
TR = 1.6; Tmax = 490;
num_sims = 50;

% Load BOLD
BOLD_dir = fullfile(base_dir, 'data/BOLD(individual)_macaque.mat');
BOLDs = importdata(BOLD_dir);
f_diff = importdata(fullfile(base_dir,'/data/Macaque_peak_freq.mat'));

% Compute FC and FCD
empFC = corr(BOLDs{1,1});
empFC = empFC - diag(diag(empFC));
empFCD = calculate_fcd(BOLDs{1,1}, TR, ceil(60/TR), 1);
triu_ind = triu(true(size(empFC)),1);

% Run Hopf model
G = findOptimalParams(FLN, f_diff, empFC, TR, Tmax, num_sims, triu_ind, empFCD, true, 0, false,[0 5]);
hopf_result = run_simulation(FLN, G, f_diff,TR,Tmax)';

% Load BEI simulation result
BEI_dir = fullfile('/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy/step1_ec_validation/part2_simulation/supplementary/BEI_results.mat');
BEI_temp = importdata(BEI_dir); 

BEI_result = zscore(squeeze(BEI_temp(1,31:end,:)), 0, 1);

%-----------------------------------------------------------%
%% Compare functional connectivity and show scatter plot
%-----------------------------------------------------------%
close all

% compute functional connectivity
hopfFC = corr(hopf_result);
BEIFC = corr(BEI_result);

% remove diagonal elements
hopfFC = hopfFC - diag(diag(hopfFC));
BEIFC = BEIFC - diag(diag(BEIFC));

% compute correlation
corr_emp_hopf = corr(empFC(triu_ind), hopfFC(triu_ind));
corr_emp_BEI = corr(empFC(triu_ind), BEIFC(triu_ind));

% Create scatter plots
figure('Position', [100 100 700 300])

% Get upper triangular values for plotting
emp_vals = empFC(triu_ind);
hopf_vals = hopfFC(triu_ind);
BEI_vals = BEIFC(triu_ind);

% First subplot: Empirical vs Hopf
subplot(1,2,1)
scatter(hopf_vals, emp_vals, 50, 'filled', ...
    'MarkerFaceColor', [0.2 0.2 0.8], ...
    'MarkerFaceAlpha', 0.5);
hold on

% Add regression line
p = polyfit(hopf_vals, emp_vals, 1);
x_fit = linspace(min(hopf_vals), max(hopf_vals), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

% Customize plot
box on
grid on
set(gca, 'LineWidth', 1.5);
set(gca, 'XTick', linspace(min(xlim), max(xlim), 5), 'YTick', linspace(min(ylim), max(ylim), 5));
% set(gca, 'XTickLabel', [], 'YTickLabel', []);
axis square

% Second subplot: Empirical vs BEI
subplot(1,2,2)
scatter(BEI_vals, emp_vals, 50, 'filled', ...
    'MarkerFaceColor', [0.2 0.6 0.2], ...
    'MarkerFaceAlpha', 0.5);
hold on

% Add regression line
p = polyfit(BEI_vals, emp_vals, 1);
x_fit = linspace(min(BEI_vals), max(BEI_vals), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

% Customize plot
box on
grid on
set(gca, 'LineWidth', 1.5);
set(gca, 'XTick', linspace(min(xlim), max(xlim), 5), 'YTick', linspace(min(ylim), max(ylim), 5));
% set(gca, 'XTickLabel', [], 'YTickLabel', []);
axis square

% Export graphics
% exportgraphics(gcf, fullfile(main_dir, 'part2_simulation/supplementary/fig1_fc_comparison.png'), 'Resolution', 1200);

%---------------------------------------------------------------%
%% Compare functional connectivity dynamics and show pdf of FCD
%---------------------------------------------------------------%
% Compute FCD
hopfFCD = calculate_fcd(hopf_result, TR, ceil(60/TR), 1);
BEIFCD = calculate_fcd(BEI_result, TR, ceil(60/TR), 1);

% Setup binning
edges = -1:0.0002:1;
x = edges(1:end-1) + diff(edges)/2;

% Create FCD distribution plots
figure('Position', [100 100 700 300])

% Get max KS points
[KS_emp_hopf, idxMax_emp_hopf] = max(abs(empFCD - hopfFCD));
xMax1 = x(idxMax_emp_hopf);
y1Max1 = empFCD(idxMax_emp_hopf);
y2Max1 = hopfFCD(idxMax_emp_hopf);

[KS_emp_BEI, idxMax_emp_BEI] = max(abs(empFCD - BEIFCD));
xMax2 = x(idxMax_emp_BEI);
y1Max2 = empFCD(idxMax_emp_BEI);
y2Max2 = BEIFCD(idxMax_emp_BEI);

subplot(1,2,1)
plot(x, empFCD, 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
hold on
plot(x, hopfFCD, 'Color', [0.2 0.2 0.8], 'LineWidth', 2);
plot([xMax1, xMax1], [y1Max1, y2Max1], '--k', 'LineWidth', 1.5);

% Customize plot
box on
grid on
set(gca, 'LineWidth', 1.5);
set(gca, 'XTick', linspace(min(xlim), max(xlim), 5), 'YTick', linspace(min(ylim), max(ylim), 5));
% set(gca, 'XTickLabel', [], 'YTickLabel', []);
axis square

% Second subplot: Empirical vs BEI FCD
subplot(1,2,2)
plot(x, empFCD, 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
hold on
plot(x, BEIFCD, 'Color', [0.2 0.6 0.2], 'LineWidth', 2);
plot([xMax2, xMax2], [y1Max2, y2Max2], '--k', 'LineWidth', 1.5);

% Customize plot
box on
grid on
set(gca, 'LineWidth', 1.5);
set(gca, 'XTick', linspace(min(xlim), max(xlim), 5), 'YTick', linspace(min(ylim), max(ylim), 5));
% set(gca, 'XTickLabel', [], 'YTickLabel', []);
axis square

% Export graphics
% exportgraphics(gcf, fullfile(main_dir, 'part2_simulation/supplementary/fig2_fcd_comparison.png'), 'Resolution', 1200);

%---------------------------------------------------------------%
%% Run EC algorithms on Hopf and BEI results
%---------------------------------------------------------------%

% Define hyperparameters
lambda = 1e+8;
threshold = 0.1;
iter = 50;
TR = 1.6;

% Define algorithms
num_algorithms = 9;
alg_names = {'var', 'fask', 'lingam'};

% Initialize results storage as arrays of structures
ec_results = struct('var', [], 'fask', [], 'lingam', []);
ec_results_hopf = repmat(ec_results, 1, 1); % Only one simulation
ec_results_BEI = repmat(ec_results, 1, 1); % Only one simulation

% Run EC algorithms
% Only one simulation, so no loop needed
% First, Hopf model
ec_results_hopf.var = run_var(hopf_result', lambda);
ec_results_hopf.fask = run_fask(hopf_result, 1e-5, threshold, iter);
ec_results_hopf.lingam = run_lingam(hopf_result, threshold, iter);

% Second, BEI model
ec_results_BEI.var = run_var(BEI_result', lambda);
ec_results_BEI.fask = run_fask(BEI_result, 1e-5, threshold, iter);
ec_results_BEI.lingam = run_lingam(BEI_result, threshold, iter);

% Save results
save(fullfile(main_dir, 'part2_simulation/supplementary/ec_macaque_results_hopf.mat'), 'ec_results_hopf');
save(fullfile(main_dir, 'part2_simulation/supplementary/ec_macaque_results_BEI.mat'), 'ec_results_BEI');

%---------------------------------------------------------------%
%% Integrate EC results using beta values from empirical data
%---------------------------------------------------------------%
% Load empirical beta values
beta_dir = fullfile(main_dir, 'part1_macaque/results/EC_betas.mat');
beta_data = importdata(beta_dir).EC_betas_vfl;
beta_vfl = mean(beta_data);
FLN = 1.2*FLN.^0.3;
% beta_vfl = Betas([2,4,7]);

% Integrate EC results using var, fask, and lingam
iec_vfl_hopf = (ec_results_hopf.var/max(ec_results_hopf.var(:)) * beta_vfl(1) + ...
                ec_results_hopf.fask/max(ec_results_hopf.fask(:)) * beta_vfl(2) + ...
                ec_results_hopf.lingam/max(ec_results_hopf.lingam(:)) * beta_vfl(3)) / 3;

iec_vfl_BEI = (ec_results_BEI.var/max(ec_results_BEI.var(:)) * beta_vfl(1) + ...
                ec_results_BEI.fask/max(ec_results_BEI.fask(:)) * beta_vfl(2) + ...
                ec_results_BEI.lingam/max(ec_results_BEI.lingam(:)) * beta_vfl(3)) / 3;

% Normalize integrated results
iec_vfl_hopf = iec_vfl_hopf / max(iec_vfl_hopf(:));
iec_vfl_BEI = iec_vfl_BEI / max(iec_vfl_BEI(:));

% Calculate correlation with FLN
corr_hopf(1) = corr(FLN(:), ec_results_hopf.var(:));
corr_hopf(2) = corr(FLN(:), ec_results_hopf.fask(:));
corr_hopf(3) = corr(FLN(:), ec_results_hopf.lingam(:));
corr_hopf(4) = corr(FLN(:), iec_vfl_hopf(:));

corr_BEI(1) = corr(FLN(:), ec_results_BEI.var(:));
corr_BEI(2) = corr(FLN(:), ec_results_BEI.fask(:));
corr_BEI(3) = corr(FLN(:), ec_results_BEI.lingam(:));
corr_BEI(4) = corr(FLN(:), iec_vfl_BEI(:));

%% Create scatter plots for Hopf and BEI models
figure('Position', [100 100 800 300]);

% Hopf model scatter plot
subplot(1, 2, 1);
scatter(iec_vfl_hopf, FLN, 40, [0.2 0.2 0.8], 'filled', ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
hold on;
p = polyfit(iec_vfl_hopf, FLN, 1);
x_range = [min(iec_vfl_hopf) max(iec_vfl_hopf)];
y_fit = polyval(p, x_range);
plot(x_range, y_fit, '--k', 'LineWidth', 1.5);
box on;
grid on;
set(gca, 'XTick', linspace(min(get(gca, 'XTick')), max(get(gca, 'XTick')), 5))
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 5))
xlim([min(iec_vfl_hopf) max(iec_vfl_hopf)])
ylim([min(FLN) max(FLN)])
% set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])

% BEI model scatter plot
subplot(1, 2, 2);
scatter(iec_vfl_BEI, FLN, 40, [0.2 0.6 0.2], 'filled', ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
hold on;
p = polyfit(iec_vfl_BEI, FLN, 1);
x_range = [min(iec_vfl_BEI) max(iec_vfl_BEI)];
y_fit = polyval(p, x_range);
plot(x_range, y_fit, '--k', 'LineWidth', 1.5);
box on;
grid on;
set(gca, 'XTick', linspace(min(get(gca, 'XTick')), max(get(gca, 'XTick')), 5))
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 5))
xlim([min(iec_vfl_BEI) max(iec_vfl_BEI)])
ylim([min(FLN) max(FLN)])
% set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])

% Export graphics
% exportgraphics(gcf, fullfile(main_dir, 'part2_simulation/supplementary/fig3_ec_comparison.png'), 'Resolution', 1200);

%---------------------------------------------------------------%
%% Human level analysis - Model comparison
%---------------------------------------------------------------%
fprintf('\nStarting human level analysis...\n');

%% Load data
% Load BEI simulation result (human)
BEI_human_dir = fullfile('/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy/step1_ec_validation/part2_simulation/supplementary/Schaefer100_results.mat');
BEI_temp = importdata(BEI_human_dir);
BEI_result_human = BEI_temp.sim_bold(:,467:end)';
sc = BEI_temp.SC;
TR_human = 0.72; 
Tmax_human = 1200;

% Load human BOLD
data_dir = '/combinelab/02_data/01_HCP/HCP_S1200/';
sublist_dir = fullfile(data_dir, 'sublist_HCP.txt');
sublist = load(sublist_dir);
f_diff_human = importdata(fullfile(base_dir, 'data/Schaefer100_peak_freq.mat'));
f_diff_human = f_diff_human(1:50);
parcellation = 'Schaefer100'; 
rois = 50;
bold_dir = '/combinelab/03_user/younghyun/03_data/fmri_data/';

% Load empirical BOLD data
sub_dir = fullfile(bold_dir, sprintf('sub%d', sublist(1)));
fmri_file = fullfile(sub_dir, sprintf('rfMRI_REST%d_LR_%s_TianLv1_hp2000_clean.csv', 1, parcellation));
data_temp = importdata(fmri_file)';
emp_bold_human = data_temp(:, 33:33+rois-1);

% Compute empirical FC and FCD
empFC_human = corr(emp_bold_human);
empFC_human = empFC_human - diag(diag(empFC_human));
empFCD_human = calculate_fcd(emp_bold_human, TR_human, ceil(60/TR_human), 1);
triu_ind_human = triu(true(size(empFC_human)), 1);

% Run Hopf model for human
% G = findOptimalParams(sc, f_diff_human, empFC_human, TR_human, Tmax_human, num_sims, triu_ind_human, empFCD_human, true, 0, false, [0 10]);
[hopfFC_human,hopfFCD_human] = batch_simulation(sc, 2.5, f_diff_human, TR_human, Tmax_human,10);

%-----------------------------------------------------------%
%% Compare functional connectivity and show scatter plot
%-----------------------------------------------------------%
close all

% Compute functional connectivity
% hopfFC_human = corr(hopf_result_human);
BEIFC_human = corr(BEI_result_human);

% Remove diagonal elements
% hopfFC_human = hopfFC_human - diag(diag(hopfFC_human));
BEIFC_human = BEIFC_human - diag(diag(BEIFC_human));

% Compute correlation
corr_emp_hopf_human = corr(empFC_human(triu_ind_human), hopfFC_human(triu_ind_human));
corr_emp_BEI_human = corr(empFC_human(triu_ind_human), BEIFC_human(triu_ind_human));

fprintf('Human FC correlation with Hopf model: %.4f\n', corr_emp_hopf_human);
fprintf('Human FC correlation with BEI model: %.4f\n', corr_emp_BEI_human);

% Create scatter plots
figure('Position', [100 100 700 300])

% Get upper triangular values for plotting
emp_vals_human = empFC_human(triu_ind_human);
hopf_vals_human = hopfFC_human(triu_ind_human);
BEI_vals_human = BEIFC_human(triu_ind_human);

% First subplot: Empirical vs Hopf
subplot(1,2,1)
scatter(hopf_vals_human, emp_vals_human, 50, 'filled', ...
    'MarkerFaceColor', [0.2 0.2 0.8], ...
    'MarkerFaceAlpha', 0.5);
hold on

% Add regression line
p = polyfit(hopf_vals_human, emp_vals_human, 1);
x_fit = linspace(min(hopf_vals_human), max(hopf_vals_human), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

% Customize plot
box on
grid on
set(gca, 'LineWidth', 1.5);
set(gca, 'XTick', linspace(min(xlim), max(xlim), 5), 'YTick', linspace(min(ylim), max(ylim), 5));
% set(gca, 'XTickLabel', [], 'YTickLabel', []);
axis square

% Second subplot: Empirical vs BEI
subplot(1,2,2)
scatter(BEI_vals_human, emp_vals_human, 50, 'filled', ...
    'MarkerFaceColor', [0.2 0.6 0.2], ...
    'MarkerFaceAlpha', 0.5);
hold on

% Add regression line
p = polyfit(BEI_vals_human, emp_vals_human, 1);
x_fit = linspace(min(BEI_vals_human), max(BEI_vals_human), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

% Customize plot
box on
grid on
set(gca, 'LineWidth', 1.5);
set(gca, 'XTick', linspace(min(xlim), max(xlim), 5), 'YTick', linspace(min(ylim), max(ylim), 5));
% set(gca, 'XTickLabel', [], 'YTickLabel', []);
axis square

% Export graphics
% exportgraphics(gcf, fullfile(main_dir, 'part2_simulation/supplementary/fig4_human_fc_comparison.png'), 'Resolution', 1200);

%---------------------------------------------------------------%
%% Compare functional connectivity dynamics and show pdf of FCD
%---------------------------------------------------------------%
% Compute FCD
% hopfFCD_human = calculate_fcd(hopf_result_human, TR_human, ceil(60/TR_human), 1);
BEIFCD_human = calculate_fcd(BEI_result_human, TR_human, ceil(60/TR_human), 1);

% Setup binning
edges = -1:0.0002:1;
x = edges(1:end-1) + diff(edges)/2;

% Calculate KS statistics
[KS_emp_hopf_human, idxMax_emp_hopf_human] = max(abs(empFCD_human - hopfFCD_human));
[KS_emp_BEI_human, idxMax_emp_BEI_human] = max(abs(empFCD_human - BEIFCD_human));

fprintf('Human FCD KS statistic with Hopf model: %.4f\n', KS_emp_hopf_human);
fprintf('Human FCD KS statistic with BEI model: %.4f\n', KS_emp_BEI_human);

% Get max KS points
xMax1_human = x(idxMax_emp_hopf_human);
y1Max1_human = empFCD_human(idxMax_emp_hopf_human);
y2Max1_human = hopfFCD_human(idxMax_emp_hopf_human);

xMax2_human = x(idxMax_emp_BEI_human);
y1Max2_human = empFCD_human(idxMax_emp_BEI_human);
y2Max2_human = BEIFCD_human(idxMax_emp_BEI_human);

% Create FCD distribution plots
figure('Position', [100 100 700 300])

subplot(1,2,1)
plot(x, empFCD_human, 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
hold on
plot(x, hopfFCD_human, 'Color', [0.2 0.2 0.8], 'LineWidth', 2);
plot([xMax1_human, xMax1_human], [y1Max1_human, y2Max1_human], '--k', 'LineWidth', 1.5);

% Customize plot
box on
grid on
set(gca, 'LineWidth', 1.5);
set(gca, 'XTick', linspace(min(xlim), max(xlim), 5), 'YTick', linspace(min(ylim), max(ylim), 5));
% set(gca, 'XTickLabel', [], 'YTickLabel', []);
axis square

% Second subplot: Empirical vs BEI FCD
subplot(1,2,2)
plot(x, empFCD_human, 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
hold on
plot(x, BEIFCD_human, 'Color', [0.2 0.6 0.2], 'LineWidth', 2);
plot([xMax2_human, xMax2_human], [y1Max2_human, y2Max2_human], '--k', 'LineWidth', 1.5);

% Customize plot
box on
grid on
set(gca, 'LineWidth', 1.5);
set(gca, 'XTick', linspace(min(xlim), max(xlim), 5), 'YTick', linspace(min(ylim), max(ylim), 5));
% set(gca, 'XTickLabel', [], 'YTickLabel', []);
axis square
% Export graphics
% exportgraphics(gcf, fullfile(main_dir, 'part2_simulation/supplementary/fig5_human_fcd_comparison.png'), 'Resolution', 1200);

%---------------------------------------------------------------%
%% Statistical tests for FC and FCD comparisons
%---------------------------------------------------------------%
% FC comparison using Fisher's z-transformation
z_hopf = atanh(corr_emp_hopf);
z_BEI = atanh(corr_emp_BEI);
z_diff = (z_hopf - z_BEI) / sqrt(1/(length(triu_ind)-3) + 1/(length(triu_ind)-3));
p_fc = 2*(1-normcdf(abs(z_diff)));

fprintf('\nFC Comparison Statistics:\n');
fprintf('Hopf FC correlation: %.4f\n', corr_emp_hopf);
fprintf('BEI FC correlation: %.4f\n', corr_emp_BEI);
fprintf('Z-score difference: %.4f\n', z_diff);
fprintf('p-value: %.4e\n', p_fc);

% FCD comparison using KS test
[~, p_fcd] = kstest2(hopfFCD, BEIFCD);

fprintf('\nFCD Comparison Statistics:\n');
fprintf('KS statistic: %.4f\n', KS_emp_hopf);
fprintf('p-value: %.4e\n', p_fcd);

% Human level statistics
z_hopf_human = atanh(corr_emp_hopf_human);
z_BEI_human = atanh(corr_emp_BEI_human);
z_diff_human = (z_hopf_human - z_BEI_human) / sqrt(1/(length(triu_ind_human)-3) + 1/(length(triu_ind_human)-3));
p_fc_human = 2*(1-normcdf(abs(z_diff_human)));

[~, p_fcd_human] = kstest2(hopfFCD_human, BEIFCD_human);

fprintf('\nHuman Level Statistics:\n');
fprintf('FC Z-score difference: %.4f\n', z_diff_human);
fprintf('FC p-value: %.4e\n', p_fc_human);
fprintf('FCD p-value: %.4e\n', p_fcd_human);
