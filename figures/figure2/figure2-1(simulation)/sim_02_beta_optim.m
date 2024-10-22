%% Beta Estimation for Effective Connectivity (EC) Algorithms
% This script estimates the beta values for different EC algorithms using 
% optimized parameters from previously saved files. The results are stored 
% and saved for further analysis.

% Clear workspace, close all figures, and clear command window
close all; clear; clc

% Load parameters of individual algorithms
Lambda = importdata('fig2_results/bestLambdas.mat', 'bestLambdas');
Alpha = importdata('fig2_results/bestO_Alphas.mat', 'bestO_Alphas');
O_Alpha = importdata('fig2_results/bestAlphas.mat', 'bestAlphas');
Threshold = importdata('fig2_results/bestThresholds.mat', 'bestThresholds');

% Number of iterations for the experiment
num_iterations = 100;

% List of atlases to be processed
atlases = {'Schaefer', 'MMP'};

% Initialize matrices to store beta values and beta histories for each atlas
EC_betas = zeros(num_iterations, 3, length(atlases));
betaHistories = cell(length(atlases), 1);

% Loop through each atlas
for atlas = 1:length(atlases)
    % Load structural connectivity matrix and frequency data for the current atlas
    atlas_name = atlases{atlas};
    SC = importdata(fullfile(sprintf('../../../data/%s_SC.mat', atlas_name)));
    f_diff = importdata(fullfile(sprintf('../../../data/%s_peak_freq.mat', atlas_name)));
    
    % Use only left hemisphere data
    N = size(SC, 1);
    f_diff = f_diff(1:N);
    
    % Set shuffle_iter based on atlas
    shuffle_iter = (atlas == 1) * 50 + (atlas == 2) * 10;
      
    % Generate network (randomize SC while preserving degree distribution)
    Js = cell(num_iterations, 1);
    parfor iter = 1:num_iterations
        Js{iter} = randmio_dir_connected(SC, shuffle_iter);
    end
    
    % Initialize matrix to store beta histories for the current atlas
    betaHistories_curr = [];
    
    % Loop over network generations
    for iter = 1:num_iterations
        J = Js{iter};
        J = J / max(J(:)) * 0.2;  % Normalize connectivity matrix

        % Calculate in-degree and out-degree for J
        in_degree_J = sum(J, 2);  % column sum
        out_degree_J = sum(J, 1);  % row sum
        
        % Run the simulation to generate BOLD signals
        BOLD = run_simulation(J, 2, f_diff);

        % Compute inferred EC using different algorithms
        rDCM = run_rdcm(BOLD', 1); 
        VAR = var_ols_ridge(BOLD, Lambda(atlas));   
        FASK = subsample_wrapper(BOLD', O_Alpha(atlas), Alpha(atlas), Threshold(atlas), 100);

        % Normalize EC matrices
        input_matrices = cell(3, 1);
        input_matrices{1} = rDCM / max(rDCM(:)) * 0.2;
        input_matrices{2} = VAR / max(VAR(:)) * 0.2;
        input_matrices{3} = FASK / max(FASK(:)) * 0.2;

        % Optimize EC to find beta values
        [Betas, betaHistory] = findBetas(input_matrices, 'Target', J);
        betaHistories_curr = cat(3, betaHistories_curr, betaHistory);                
        EC_betas(iter, :, atlas) = Betas;
    end
    
    % Store beta histories for the current atlas
    betaHistories{atlas} = betaHistories_curr;
end

% Save the estimated beta values and beta histories
save('fig2_results/EC_beta_values.mat', 'EC_betas')
save('fig2_results/EC_beta_histories.mat', 'betaHistories')
