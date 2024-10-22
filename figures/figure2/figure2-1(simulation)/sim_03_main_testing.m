%% Main Testing of iEC Framework
% This script tests the integrated Effective Connectivity (iEC) framework 
% by loading optimized parameters and simulating BOLD signals. The script 
% computes various correlation metrics and saves the results.

% Clear workspace, close all figures, and clear command window
clear; close all; clc

% Load parameters of individual algorithms
Lambda = importdata('fig2_results/bestLambdas.mat', 'bestLambdas');
Alpha = importdata('fig2_results/bestO_Alphas.mat', 'bestO_Alphas');
O_Alpha = importdata('fig2_results/bestAlphas.mat', 'bestAlphas');
Threshold = importdata('fig2_results/bestThresholds.mat', 'bestThresholds');
EC_betas = importdata('fig2_results/EC_beta_values.mat', 'EC_betas');

% Set the number of iterations
num_iterations = 100;

% List of atlases to be processed
atlases = {'Schaefer', 'MMP'};

% Allocate the size of the result matrices by adding threshold dimension
in_degree_corr_values = zeros(num_iterations, 4, length(atlases));
out_degree_corr_values = zeros(num_iterations, 4, length(atlases));
degree_corr_values = zeros(num_iterations, 4, length(atlases));
EC_correlation_values = zeros(num_iterations, 4, length(atlases));

% Compute mean beta values for each atlas
Betas_all = mean(EC_betas, 1);

% Loop through each atlas
for atlas = 1:length(atlases)
    % Load structural connectivity matrix and frequency data
    atlas_name = atlases{atlas};
    SC = importdata(fullfile(sprintf('../../../data/%s_SC.mat', atlas_name)));
    f_diff = importdata(fullfile(sprintf('../../../data/%s_peak_freq.mat', atlas_name)));
    Betas = Betas_all(:, :, atlas);
    
    % Use only left hemisphere data
    N = size(SC, 1);
    f_diff = f_diff(1:N);
    
    % Generate network (randomize SC while preserving degree distribution)
    Js = cell(num_iterations, 1);
    parfor iter = 1:num_iterations
        Js{iter} = randmio_dir_connected(SC, 50);
    end
    
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
        [~,GC] = StateSpaceGC(BOLD,0,Lambda(atlas));

        % Normalize EC matrices
        rDCM = rDCM / max(rDCM(:)) * 0.2;
        VAR = VAR / max(VAR(:)) * 0.2;
        FASK = FASK / max(FASK(:)) * 0.2;
        GC = GC / max(GC(:))*0.2;
        
        % Compute integrated EC (iEC)
        iEC = (Betas(1) * rDCM + Betas(2) * VAR + Betas(3) * FASK) / 3;
        
        % Calculate correlations between J and inferred EC matrices
        correlation_VAR = corr(J(:), VAR(:));
        correlation_GC = corr(J(:), GC(:));
        correlation_rDCM = corr(J(:), rDCM(:));
        correlation_FASK = corr(J(:), FASK(:));
        correlation_iEC = corr(J(:), iEC(:));

        % Calculate in-degree and out-degree for each method
        in_degree_VAR = sum(VAR, 2);  
        out_degree_VAR = sum(VAR, 1);  

        in_degree_rDCM = sum(rDCM, 2);  
        out_degree_rDCM = sum(rDCM, 1);  

        in_degree_FASK = sum(FASK, 2);  
        out_degree_FASK = sum(FASK, 1);  

        in_degree_iEC = sum(iEC, 2);  
        out_degree_iEC = sum(iEC, 1);  

        % Calculate correlation of in-degree and out-degree with J
        in_degree_corr_J_VAR = corr(in_degree_J, in_degree_VAR);
        in_degree_corr_J_rDCM = corr(in_degree_J, in_degree_rDCM);
        in_degree_corr_J_FASK = corr(in_degree_J, in_degree_FASK);
        in_degree_corr_J_iEC = corr(in_degree_J, in_degree_iEC);

        out_degree_corr_J_VAR = corr(out_degree_J', out_degree_VAR');
        out_degree_corr_J_rDCM = corr(out_degree_J', out_degree_rDCM');
        out_degree_corr_J_FASK = corr(out_degree_J', out_degree_FASK');
        out_degree_corr_J_iEC = corr(out_degree_J', out_degree_iEC');

        % Calculate average degree correlations
        degree_corr_J_VAR = mean([in_degree_corr_J_VAR, out_degree_corr_J_VAR], 2);
        degree_corr_J_rDCM = mean([in_degree_corr_J_rDCM, out_degree_corr_J_rDCM], 2);
        degree_corr_J_FASK = mean([in_degree_corr_J_FASK, out_degree_corr_J_FASK], 2);
        degree_corr_J_iEC = mean([in_degree_corr_J_iEC, out_degree_corr_J_iEC], 2);

        % Store the results
        in_degree_corr_values(iter, :, atlas) = [in_degree_corr_J_rDCM, in_degree_corr_J_VAR, in_degree_corr_J_FASK, in_degree_corr_J_iEC];
        out_degree_corr_values(iter, :, atlas) = [out_degree_corr_J_rDCM, out_degree_corr_J_VAR, out_degree_corr_J_FASK, out_degree_corr_J_iEC];
        degree_corr_values(iter, :, atlas) = [degree_corr_J_rDCM, degree_corr_J_VAR, degree_corr_J_FASK, degree_corr_J_iEC];
        EC_correlation_values(iter, :, atlas) = [correlation_rDCM, correlation_VAR, correlation_FASK, correlation_iEC];
    end
end

% Save the results
save('fig2_results/EC_correlation_values.mat', 'EC_correlation_values')
save('fig2_results/EC_in_degree_corr_values.mat', 'in_degree_corr_values')
save('fig2_results/EC_out_degree_corr_values.mat', 'out_degree_corr_values')
save('fig2_results/EC_degree_corr_values.mat', 'degree_corr_values')

