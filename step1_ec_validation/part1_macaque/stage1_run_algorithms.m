%% Stage 1: Run EC algorithms
% This script runs effective connectivity estimation on macaque fMRI data using
% multiple algorithms on both original and deconvolved BOLD signals.
%
% Main workflow:
% 1. Loads FLN ground truth and individual macaque BOLD signals (19 subjects)
% 2. Performs blind deconvolution using rsHRF to remove hemodynamic response
% 3. Runs 9 EC algorithms: rDCM, VAR, GC, FASK, CCD, BOSS, LiNGAM, GRASP, Patel
% 4. Saves EC matrices for each subject and algorithm for subsequent integration
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
num_regions = size(FLN,1);
num_subjects = 19;

% Load BOLD signals
BOLD_dir = fullfile(main_dir, 'data/BOLD(individual)_macaque.mat');
BOLDs = importdata(BOLD_dir);
TR = 1.6;
Tmax = size(BOLDs{1,1},1);

%--------------------------------------------------------------------------%
%% Run individual EC algorithms
%--------------------------------------------------------------------------%
% Set parameters for blind deconvolution
para.TR = TR;
para.name = 'Canonical HRF (with time and dispersion derivatives)';
temporal_mask = [];
para.T = 3;
para.T0 = 1;
para.order = 3;
para.dt = para.TR/para.T;
para.AR_lag = 1;
para.thr = 1;
para.len = 24;
min_onset_search = 4;
max_onset_search = 8;
para.lag = fix(min_onset_search/para.dt):fix(max_onset_search/para.dt);
T = round(para.len/TR);
data_deconv = zeros(Tmax,num_regions,num_subjects);

% Loop over subjects
for subj = 1:num_subjects
    data = zscore(BOLDs{subj},0,1);
    
    sigma = std(data);
    
     [beta_hrf, bf, ~] = rsHRF_estimation_temporal_basis(data,para,temporal_mask);
     hrf = bf*beta_hrf(1:size(bf,2),:);
     hrf = resample(hrf,1,para.T);
     for roi = 1:num_regions
         data_deconv(:,roi,subj) = zscore(rsHRF_iterative_wiener_deconv(data(:,roi),hrf(:,roi),80));  
     end
end

%% Initialize constants
lambda = 1e+8;
threshold = 0.1;
iter = 50;

% Initialize results storage as arrays of structures
ec_results = struct('rdcm', [], 'var', [], 'gc', [], 'fask', [], ...
    'ccd', [], 'boss', [], 'lingam', [], 'grasp', [], 'patel', []);
ec_results = repmat(ec_results, num_subjects, 1);

% Determine whether to use deconvolved data or not
use_deconv = false;

% Loop over subjects
for subj = 1:num_subjects
    if use_deconv
        BOLD = data_deconv(:,:,subj);
    else
        BOLD = BOLDs{subj};
    end
    
    % run EC algorithms
    ec_results(subj).rdcm = run_rdcm(BOLD, TR);
    ec_results(subj).var = run_var(BOLD', lambda);
    ec_results(subj).gc = run_gc(BOLD', 0, lambda);
    ec_results(subj).fask = run_fask(BOLD, 1e-6, threshold, iter);
    ec_results(subj).ccd = run_ccd(BOLD, threshold, iter);
    ec_results(subj).boss = run_boss(BOLD, threshold, iter);
    ec_results(subj).lingam = run_lingam(BOLD, threshold, iter);
    ec_results(subj).grasp = run_grasp(BOLD, threshold, iter);
    ec_results(subj).patel = run_patel(BOLD, threshold, iter);

end

% Save results
if use_deconv
    save(fullfile(main_dir, 'results/ec_results_deconv.mat'), 'ec_results');
else
    save(fullfile(main_dir, 'results/ec_results.mat'), 'ec_results');
end

%--------------------------------------------------------------------------%
%% Save group average ECs
%--------------------------------------------------------------------------%
group_mean_results = struct();

algorithms = fieldnames(ec_results);
for k = 1:length(algorithms)
    alg = algorithms{k};
    all_subjects_data = arrayfun(@(x) x.(alg), ec_results, 'UniformOutput', false);
    all_subjects_data = cat(3, all_subjects_data{:}); % Concatenate along the third dimension
    all_subjects_data = all_subjects_data/max(all_subjects_data(:));
    
    % Calculate the mean across the third dimension (subjects)
    group_mean_results.(alg) = mean(all_subjects_data, 3);
end

% Create individual connectivity matrix figures
figure('Position', [100 100 200 200]);  % Square figure for square plots

for k = 1:length(algorithms)
    alg = algorithms{k};
    
    % Clear figure for each new plot
    clf;
    
    % Plot connectivity matrix
    imagesc(group_mean_results.(alg));
    colormap(gca, generateColorMap(reshape(group_mean_results.(alg),1,[]),100));
    
    % Customize plot
    axis square;
    set(gca, 'LineWidth', 1.5, 'Box', 'on', ...
        'XTick', [], 'YTick', [], 'XTickLabel', [], 'YTickLabel', []); % Remove ticks and labels
    
    % Adjust axes position to be centered in figure
    pos = get(gca, 'Position');
    pos(1) = (1-pos(3))/2;  % Center horizontally
    pos(2) = (1-pos(4))/2;  % Center vertically
    set(gca, 'Position', pos);
    
    % Save individual figure
    exportgraphics(gcf, fullfile(main_dir, 'supplementary', ['conn_matrix_', alg, '.png']), ...
        'Resolution', 300, 'BackgroundColor', 'white');
end
