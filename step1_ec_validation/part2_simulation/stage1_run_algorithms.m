%% stage1_run_algorithms.m
% This script performs Hopf model simulations and effective connectivity estimation
% for validation of the iEC framework using directed structural connectivity matrices.
% 
% Main workflow:
% 1. Lambda optimization: Tests VAR lambda parameters (1e1 to 1e8) across 10 iterations
%    to find optimal values for each atlas (Schaefer100, MMP360)
% 2. Directed SC generation: Creates 100 directed SC variants using randmio_dir_connected
% 3. Hopf simulation: Generates synthetic BOLD signals from directed SCs
% 4. EC estimation: Applies 8 algorithms (rDCM, VAR, FASK, CCD, BOSS, LiNGAM, GRASP, Patel)
%    to estimate effective connectivity from synthetic BOLD data
% 5. Results storage: Saves EC estimates and directed SC matrices for validation
%
% Author: Younghyun Oh
% Date: 2025-09-12
% Version: 1.0

% Clear workspace
clear; close all; clc

% Load structural connectivity matrix and frequency data
main_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy';
result_dir = fullfile(main_dir, 'step1_ec_validation/part2_simulation/results');

% List of atlases to be processed
atlases = {'Schaefer100', 'MMP360'};

% Modified parameters
num_iterations = 10;  % Reduced to 10 for lambda optimization
lambda_range = logspace(1, 8, 20);  % Create range of lambda values from 1e1 to 1e8

% Initialize matrices to store optimal lambdas
optimal_lambdas = cell(length(atlases), 1);

% Simulation parameters
TR = 0.72;
Tmax = 1200; 

% Loop over each atlas
for atlas = 1:length(atlases)
    parcellation = atlases{atlas};
    SC = importdata(fullfile(main_dir, 'data', [parcellation, '_SC.mat']));
    freq = importdata(fullfile(main_dir, 'data', [parcellation, '_peak_freq.mat']));

    % Use single cortex
    N = size(SC, 1);
    freq = freq(1:N);


    % Store optimal lambda for each iteration
    atlas_optimal_lambdas = zeros(num_iterations, 1);
    
    % Pre-generate all Js using parfor
    Js = cell(num_iterations, 1);
    parfor iter = 1:num_iterations
        Js{iter} = randmio_dir_connected(SC, 1);
    end

    % Now use the pre-generated Js for lambda optimization
    for iter = 1:num_iterations
        J = Js{iter};
        
        % Run simulation
        BOLD = run_simulation(J, 1, freq, TR, Tmax)';

        % Test different lambda values
        best_corr = -inf;
        best_lambda = nan;

        for lambda = lambda_range
            var_result = run_var(BOLD', lambda);
            corr_val = corr(J(:), var_result(:));
            
            if corr_val > best_corr
                best_corr = corr_val;
                best_lambda = lambda;
            end
        end

        atlas_optimal_lambdas(iter) = best_lambda;
    end

    % Calculate median and round to first digit of 1
    median_lambda = median(atlas_optimal_lambdas);
    exp_val = round(log10(median_lambda));
    optimal_lambdas{atlas} = 10^exp_val;
end

% Save optimal lambda values
save(fullfile(result_dir, 'optimal_lambdas.mat'), 'optimal_lambdas');

%-------------------------------------------------------------------------%
%% Use optimal lambda values for final analysis
%-------------------------------------------------------------------------%

% Number of iterations for the experiment
num_iterations = 100;
num_algorithms = 8;  


% Initialize matrices to store beta values and beta histories for each atlas
ec_results = struct('rdcm', [], 'var', [], 'fask', [], ...
    'ccd', [], 'boss', [], 'lingam', [], 'grasp', [], 'patel', []);  % Removed 'gc' field
ec_results = repmat(ec_results, length(atlases), num_iterations);
Js = cell(length(atlases), num_iterations);

% Hyperparameters of algorithms
threshold = 0.1;

% Loop over each atlas
for atlas = 1:length(atlases)
    parcellation = atlases{atlas};
    SC = importdata(fullfile(main_dir, 'data', [parcellation, '_SC.mat']));
    freq = importdata(fullfile(main_dir, 'data', [parcellation, '_peak_freq.mat']));

    % Use single cortex
    N = size(SC, 1);
    freq = freq(1:N);

    switch atlas
        case 1
            var_lambda = 1e+3;
        case 2
            var_lambda = 1e+4;
    end

    % Generate directional SC
    parfor iter = 1:num_iterations
        Js{atlas, iter} = randmio_dir_connected(SC, 1);
    end

    % Loop over each network generation
    for iter = 1:num_iterations
        J = Js{atlas, iter};
        J = J/max(J(:));

        % Run the simulation to generate BOLD signals
        BOLD = run_simulation(J, 1, freq, TR, Tmax)';

        % Run EC algorithms
        ec_results(atlas, iter).rdcm = run_rdcm(BOLD, 0.72);
        ec_results(atlas, iter).var = run_var(BOLD', var_lambda);
        ec_results(atlas, iter).fask = run_fask(BOLD, 1e-6, threshold, 50);
        ec_results(atlas, iter).ccd = run_ccd(BOLD, threshold, 50);
        ec_results(atlas, iter).boss = run_boss(BOLD, threshold, 50);
        ec_results(atlas, iter).lingam = run_lingam(BOLD, threshold, 50);
        ec_results(atlas, iter).grasp = run_grasp(BOLD, threshold, 50);
        ec_results(atlas, iter).patel = run_patel(BOLD, threshold, 50);

    end 

end

save(fullfile(result_dir, 'ec_results.mat'), 'ec_results');
save(fullfile(result_dir, 'directed_SC.mat'), 'Js');
