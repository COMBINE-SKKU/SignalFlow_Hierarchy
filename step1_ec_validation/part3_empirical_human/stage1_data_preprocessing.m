%% Stage 1: Data preprocessing and group averaging
% This script preprocesses human empirical fMRI data from HCP dataset and
% creates group-averaged connectivity matrices for iEC framework validation.
%
% Main workflow:
% 1. EC preprocessing: Computes group-average EC matrices for train/test sets (8 algorithms)
% 2. BOLD data loading: Loads individual subject BOLD signals from HCP dataset
% 3. FC computation: Calculates functional connectivity matrices for all subjects
% 4. FCD computation: Calculates functional connectivity dynamics distributions
% 5. Time delay analysis: Computes lag matrices and principal components
% 6. Data splitting: Separates data into training (220) and testing (220) sets
% 7. Group averaging: Creates group-level statistics for empirical validation
%
% Output: Preprocessed group data files for subsequent iEC analysis
%
% Author: Younghyun Oh
% Date: 2025-02-06
%--------------------------------------------------------------------------%
%% Set up and preprocess data
%--------------------------------------------------------------------------%
% Clear workspace
clear; close all; clc;

% Set constants that will be used throughout the script
parcellation = 'MMP360';
n_subjects = 440;
n_train = n_subjects/2;
n_test = n_subjects - n_train;
Tmax = 2400; TR = 0.72; num_sims = 10;
alpha = 0.05;

% Set directories
main_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy';
data_dir = '/combinelab/02_data/01_HCP/HCP_S1200/';
ec_dir = fullfile(main_dir, 'step1_ec_validation', 'part3_empirical_human');
f_diff = importdata(fullfile(main_dir,'/data/',[parcellation,'_peak_freq.mat']));
algorithm_names = {'rdcm','var','fask','ccd','boss','lingam','grasp','patel'};
switch parcellation
    case 'Schaefer100'
        rois = 100;
    case 'MMP360'
        rois = 360;
end
ec_results = importdata(fullfile(ec_dir, 'data', [parcellation,'_ec_results.mat']));

%% Compute group average EC matrices for training and testing sets

% Initialize cell arrays to store results
n_alg = length(algorithm_names);
ec_alg_groups_train = cell(1,n_alg);
ec_alg_groups_test = cell(1,n_alg);

% Loop over algorithms
for algorithm = 1:n_alg
    % Initialize matrices for both sets
    ec_alg_test = zeros(rois,rois,n_train);  
    % Create a matrix for all subjects to compute p0 once
    ec_alg_whole = zeros(rois,rois,n_subjects);
    
    % Extract data for all subjects first
    for sub = 1:n_subjects
        ec_alg_whole(:,:,sub) = ec_results(sub).(algorithm_names{algorithm})(33:end,33:end);
    end
    
    % Then separate into training and test sets
    ec_alg_train = ec_alg_whole(:,:,1:n_train);
    for sub = 1:n_train
        test_idx = sub + n_train;
        ec_alg_test(:,:,sub) = ec_alg_whole(:,:,test_idx);
    end
    
    % Compute group average
    ec_alg_group_train = mean(ec_alg_train, 3);
    ec_alg_group_test = mean(ec_alg_test, 3);
    
    % Apply thresholding for non-rDCM and non-VAR algorithms
    % For other algorithms, remove edges detected in less than 10% of subjects
    current_alg_name = algorithm_names{algorithm};
    if ~strcmp(current_alg_name, 'rdcm') && ~strcmp(current_alg_name, 'var')
        % Calculate proportion of subjects with non-zero edges
        prop_train = mean(ec_alg_train ~= 0, 3);
        prop_test = mean(ec_alg_test ~= 0, 3);
        
        % Apply 0.1 threshold (remove edges detected in <10% of subjects)
        ec_alg_group_train(prop_train < 0.1) = 0;
        ec_alg_group_test(prop_test < 0.1) = 0;
    end
    
    % Normalize the resulting matrices
    ec_alg_group_train = ec_alg_group_train/max(ec_alg_group_train(:))*0.1;
    ec_alg_group_test = ec_alg_group_test/max(ec_alg_group_test(:))*0.1;
    
    % Store the results in the cell arrays
    ec_alg_groups_train{algorithm} = ec_alg_group_train;
    ec_alg_groups_test{algorithm} = ec_alg_group_test;
end

% Save results
save(fullfile(ec_dir, 'data', [parcellation,'_ec_alg_groups_train.mat']), 'ec_alg_groups_train');
save(fullfile(ec_dir, 'data', [parcellation,'_ec_alg_groups_test.mat']), 'ec_alg_groups_test');

%% Load individual BOLD data and compute group average FC and FCD for train and test sets
% Load empirical BOLD data
sublist_dir = fullfile(data_dir, 'sublist_HCP.txt');
sublist = load(sublist_dir);
bold_dir = '/combinelab/03_user/younghyun/03_data/fmri_data/';

% Load BOLD data for each subject
fmri_data_all = cell(1, n_subjects);
fc_all = zeros(rois,rois,n_subjects);
fcd_all = cell(n_subjects,1);
fcd_pdf_all = cell(n_subjects,1);
td_matrix_all = zeros(rois, rois, n_subjects); % Initialize time delay matrix array

parfor i = 1:n_subjects
    sub_dir = fullfile(bold_dir, sprintf('sub%d', sublist(i)));
    fmri_file = fullfile(sub_dir, sprintf('rfMRI_REST%d_LRRL_%s_TianLv1_hp2000_clean.csv', 1, parcellation));
    data_temp = importdata(fmri_file)';
    fmri_data_all{i} = data_temp(:,33:33+rois-1);

    % Compute FC and FCD
    fc_temp = corr(fmri_data_all{i});
    fc_temp(eye(size(fc_temp)) == 1) = 0;
    fc_all(:,:,i) = atanh(fc_temp);
    [~,~,fcd_pdf] = calculate_fcd(fmri_data_all{i}, TR, ceil(60/TR), 1);
    fcd_pdf_all{i} = fcd_pdf;
    
    % Calculate time delay for this subject
    td_matrix_all(:,:,i) = calculate_lag_threads(fmri_data_all{i}, TR);
end

% Separate into training and testing sets
fc_train = fc_all(:,:,1:n_train);
fc_test = fc_all(:,:,n_train+1:end);
fcd_pdf_train = fcd_pdf_all(1:n_train);
fcd_pdf_test = fcd_pdf_all(n_train+1:end);
td_matrix_train = td_matrix_all(:,:,1:n_train);
td_matrix_test = td_matrix_all(:,:,n_train+1:end);

% Compute group average FC for training set
fc_emp_train = tanh(nanmean(fc_train,3));

% Compute group average FC for testing set
fc_emp_test = tanh(nanmean(fc_test,3));

% Calculate average time delay matrices for training and test sets
mean_td_train = mean(td_matrix_train, 3);
mean_td_test = mean(td_matrix_test, 3);
mean_lag_train = mean(mean_td_train, 1);
mean_lag_test = mean(mean_td_test, 1);

% Calculate principal components for time delay matrices
[coeff_train, ~, ~] = svd(mean_td_train);
PC1_lag_train = coeff_train(:,1);
PC2_lag_train = coeff_train(:,2);

[coeff_test, ~, ~] = svd(mean_td_test);
PC1_lag_test = coeff_test(:,1);
PC2_lag_test = coeff_test(:,2);

% Compute group average FCD for training set
fcd_train_concat = [];
for i = 1:n_train
    fcd_train_concat = [fcd_train_concat; fcd_pdf_train{i}];
end

edges = -1:0.0002:1;
pdf_group = histcounts(fcd_train_concat, edges, 'Normalization', 'pdf');
cdf_group = cumsum(pdf_group);
fcd_emp_train = cdf_group / cdf_group(end);

% Compute group average FCD for testing set
fcd_test_concat = [];
for i = 1:n_train
    fcd_test_concat = [fcd_test_concat; fcd_pdf_test{i}];
end
edges = -1:0.0002:1;
pdf_group = histcounts(fcd_test_concat, edges, 'Normalization', 'pdf');
cdf_group = cumsum(pdf_group);
fcd_emp_test = cdf_group / cdf_group(end);

% Save the empirical FC, FCD, and time delay data for both sets
save(fullfile(ec_dir, 'data', [parcellation,'_fc_fcd_td_train.mat']), 'fc_emp_train', 'fcd_emp_train', 'mean_td_train', 'mean_lag_train', 'PC1_lag_train', 'PC2_lag_train');
save(fullfile(ec_dir, 'data', [parcellation,'_fc_fcd_td_test.mat']), 'fc_emp_test', 'fcd_emp_test', 'mean_td_test', 'mean_lag_test', 'PC1_lag_test', 'PC2_lag_test');