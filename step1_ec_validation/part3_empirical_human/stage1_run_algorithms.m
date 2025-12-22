%% run_sch_algs.m
% This script runs various effective connectivity (EC) algorithms on fMRI data.
% It processes data from multiple subjects, applying selected EC algorithms and
% saving the results. The script supports optional HRF deconvolution and allows
% for selective execution of algorithms.
%
% Author: Younghyun Oh
% Date: 2023

%% Set up constants and configuration
%Clear workspace
% clear; close all; clc

% Data paths
fmri_dir = '/combinelab/03_user/younghyun/03_data/fmri_data';
result_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy/step1_ec_validation/part3_empirical_human/data';

% Algorithm selection (set to 1 to run, 0 to skip) - all enabled by default
run_rdcm_alg = 1;
run_var_alg = 1;
run_gc_alg = 1;
run_fask_alg = 1;
run_ccd_alg = 1;
run_boss_alg = 1;
run_lingam_alg = 1;
run_grasp_alg = 1;
run_patel_alg = 1;

% Subject range
start_subject = 1;
end_subject = 440;

% Timing measurement (1 to measure execution time, 0 to skip)
measure_timing = 1;

% Deconvolution setting (1 to perform deconvolution, 0 to skip)
do_deconv = 0;

% Algorithm parameters
threshold = 0.1;
iter = 100;
o_alpha = 1e-6;

% fMRI parameters
TR = 0.72;
Tmax = 2400;
ses = 1;
parcellation = 'MMP360';
switch parcellation
    case 'Schaefer100'
        rois = 100;
        var_lambda = 500;
    case 'MMP360'
        rois = 360;
        var_lambda = 1000;
end

% % Set up the parallel cluster with the optimal number of workers
% pc = parcluster('local');
% pc.NumWorkers = 50;
% saveProfile(pc);

% % Check if parallel pool is already open, if not open with the optimal worker count
% if isempty(gcp('nocreate'))
%     parpool(pc, 50);
%     fprintf('Opened parallel pool with %d workers\n', 50);
% else
%     currentPool = gcp('nocreate');
%     fprintf('Using existing parallel pool with %d workers\n', currentPool.NumWorkers);
%     % Optionally restart the pool if the worker count is significantly different
%     if currentPool.NumWorkers < (50 * 0.75)
%         fprintf('Current pool has significantly fewer workers than optimal. Restarting pool.\n');
%         delete(currentPool);
%         parpool(pc, 50);
%         fprintf('Opened new parallel pool with %d workers\n', 50);
%     end
% end

%% Load subject list and initialize structures
sublist = load('/combinelab/02_data/01_HCP/HCP_S1200/sublist_HCP.txt');
num_subjects = end_subject;

% Create field names for enabled algorithms
field_names = {};
if run_rdcm_alg, field_names{end+1} = 'rdcm'; end
if run_var_alg, field_names{end+1} = 'var'; end
if run_gc_alg, field_names{end+1} = 'gc'; end
if run_fask_alg, field_names{end+1} = 'fask'; end
if run_ccd_alg, field_names{end+1} = 'ccd'; end
if run_boss_alg, field_names{end+1} = 'boss'; end
if run_lingam_alg, field_names{end+1} = 'lingam'; end
if run_grasp_alg, field_names{end+1} = 'grasp'; end
if run_patel_alg, field_names{end+1} = 'patel'; end

% Initialize results structure with empty arrays for each algorithm
ec_results = struct();
for i = 1:length(field_names)
    ec_results.(field_names{i}) = [];
end
ec_results = repmat(ec_results, num_subjects, 1);

% Initialize timing structure if needed
if measure_timing
    timing_results = struct();
    for i = 1:length(field_names)
        timing_results.(field_names{i}) = [];
    end
    timing_results = repmat(timing_results, num_subjects, 1);
end

%% Set up HRF deconvolution parameters
para = struct('TR', TR, 'Tmax', Tmax, ...
              'name', 'Canonical HRF (with time and dispersion derivatives)', ...
              'T', 3, 'T0', 1, 'order', 3, 'AR_lag', 1, 'thr', 1, 'len', 24);
para.dt = para.TR/para.T;
min_onset_search = 4;
max_onset_search = 8;
para.lag = fix(min_onset_search/para.dt):fix(max_onset_search/para.dt);
temporal_mask = [];

%% Run algorithms for each subject
% csv_results_dir = fullfile(result_dir, sprintf('py_tetrad_results_%s', parcellation));
% Standard processing for other workstations
for subject = start_subject:end_subject
    subject_number = sublist(subject);
    fprintf('Processing subject %d (%d/%d)\n', subject_number, subject-start_subject+1, end_subject-start_subject+1);
    
    % Load fMRI data
    sub_dir = fullfile(fmri_dir, sprintf('sub%d', subject_number));
    fmri_file = fullfile(sub_dir, sprintf('rfMRI_REST%d_LRRL_%s_TianLv1_hp2000_clean.csv', ses, parcellation));
    fmri_data = importdata(fmri_file)';
    
   
    % Perform deconvolution if enabled
    if do_deconv
        data_deconv = zeros(Tmax, rois);
        [beta_hrf, bf, ~] = rsHRF_estimation_temporal_basis(fmri_data, para, temporal_mask);
        hrf = bf*beta_hrf(1:size(bf,2),:);
        hrf = resample(hrf, 1, para.T);
        for roi = 1:rois
            data_deconv(:,roi) = zscore(rsHRF_iterative_wiener_deconv(fmri_data(:,roi), hrf(:,roi), 100));  
        end
        analysis_data = data_deconv;
    else
        analysis_data =fmri_data;
    end

    % Run selected EC algorithms with optional timing
    % Create a structure to map algorithm flags to their execution functions
    alg_map = struct();
    if run_rdcm_alg, alg_map.rdcm = @() run_rdcm(analysis_data, TR); end
    if run_var_alg, alg_map.var = @() run_var(analysis_data', var_lambda); end
    if run_gc_alg, alg_map.gc = @() run_gc(analysis_data', 0, var_lambda); end
    if run_fask_alg, alg_map.fask = @() run_fask(analysis_data, o_alpha, threshold, iter); end
    if run_ccd_alg, alg_map.ccd = @() run_ccd(analysis_data, threshold, iter); end
    if run_boss_alg, alg_map.boss = @() run_boss(analysis_data, threshold, iter); end
    if run_lingam_alg, alg_map.lingam = @() run_lingam(analysis_data, threshold, iter); end
    if run_grasp_alg, alg_map.grasp = @() run_grasp(analysis_data, threshold, iter); end
    if run_patel_alg, alg_map.patel = @() run_patel(analysis_data, threshold, iter); end
    
    % Execute each algorithm
    alg_names = fieldnames(alg_map);
    for i = 1:length(alg_names)
        alg_name = alg_names{i};
        fprintf('  Running %s algorithm\n', alg_name);

        % % Check if csv file exists
        % csv_file = fullfile(csv_results_dir, sprintf('sub%d_%s.csv', subject_number, upper(alg_name)));
        % if exist(csv_file, 'file')
        %     fprintf('CSV file already exists for subject %d: %s\n', subject_number, csv_file);
        %     continue;
        % end
        
        if measure_timing
            tic;
            ec_results(subject).(alg_name) = alg_map.(alg_name)();
            timing_results(subject).(alg_name) = toc;
        else
            ec_results(subject).(alg_name) = alg_map.(alg_name)();
        end

        csv_file = fullfile(csv_results_dir, sprintf('sub%d_%s.csv', subject_number, upper(alg_name)));
        writematrix(ec_results(subject).(alg_name), csv_file);
        fprintf('    Saved CSV for %s: %s\n', alg_name, csv_file);
    end
end



%% Save results with appropriate filename based on deconvolution setting
if do_deconv
    save_filename = fullfile(csv_results_dir, sprintf('ec_results_deconv_%s.mat', parcellation));
else
    save_filename = fullfile(csv_results_dir, sprintf('%s_ec_results.mat', parcellation));
end

if measure_timing
    save(save_filename, 'ec_results', 'timing_results');
else
    save(save_filename, 'ec_results');
end

fprintf('Results saved to: %s\n', save_filename);
fprintf('Processed subjects %d to %d\n', start_subject, end_subject);
