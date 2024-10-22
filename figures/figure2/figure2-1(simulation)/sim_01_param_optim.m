%% Parameter Optimization of Individual Algorithms
% This script optimizes parameters for individual algorithms using structural connectivity (SC) matrices
% and frequency data from different atlases. The results are stored and saved for further analysis.

% Clear workspace, close all figures, and clear command window
clear; close all; clc

% Number of optimization iterations
nIter = 10;

% List of atlases to be processed
atlases = {'Schaefer', 'MMP'};

% Initialize matrices to store best parameters for each atlas and iteration
bestLambdas = zeros(length(atlases), nIter);
bestAlphas = zeros(length(atlases), nIter);
bestO_Alphas = zeros(length(atlases), nIter);
bestThresholds = zeros(length(atlases), nIter);

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
    Js = cell(nIter, 1);
    parfor iter = 1:nIter
        Js{iter} = randmio_dir_connected(SC, shuffle_iter);
    end
    
    % Loop through iterations for optimization
    for iter = 1:nIter
        
        J = Js{iter};
        J = J / max(J(:)) * 0.2; % Normalize connectivity matrix
        
        % Run the simulation to generate BOLD signals
        BOLD = run_simulation(J, 2, f_diff);
        
        % Perform parameter optimization for VAR and FASK models
        bestLambda = VAR_optimization(BOLD, J);
        [bestAlpha, bestO_Alpha, bestThreshold] = FASK_optimization(BOLD', J);
        
        % Store the best parameters for this iteration
        bestLambdas(atlas, iter) = bestLambda;
        bestAlphas(atlas, iter) = bestAlpha;
        bestO_Alphas(atlas, iter) = bestO_Alpha;
        bestThresholds(atlas, iter) = bestThreshold;
    end
end

% Calculate median parameter values for each atlas
Lambda = median(bestLambdas, 2);
Alpha = median(bestAlphas, 2);
O_Alpha = median(bestO_Alphas, 2);
Threshold = median(bestThresholds, 2);

% Save the best parameters to files for further analysis
save('fig2_results/bestLambdas.mat', 'bestLambdas');
save('fig2_results/bestO_Alphas.mat', 'bestO_Alphas');
save('fig2_results/bestAlphas.mat', 'bestAlphas');
save('fig2_results/bestThresholds.mat', 'bestThresholds');
