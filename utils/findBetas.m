function [Betas, BetaHistory] = findBetas(inputMats, varargin)
% FINDBETAS finds optimal beta values to integrate algorithms
%   [BETAS, BETAHISTORY] = FINDBETAS(INPUTMATS, ...) returns
%   N-by-1 vector BETAS containing the estimated optimal beta values of
%   individual algorithms and N-by-Iter matrix containing history of beta
%   values.
%
%   Parameters:
%   'Target'    The ground truth (or its proxy) matrix to compute 
%               Pearson's linear correlation coefficient.
%   'Atlas'     The atlas resolution of the ground truth (e.g., 'MMP' or
%               'Schaefer'). This information is used to extract the
%               peak frequency information of each region of the atlas.
%   'State'     The state information of the fMRI data (e.g., 'resting' or
%               'pain'). This information is also used to set the right
%               file path for the peak frequency information obtained
%               from each state.
%
%   This function requires the 'Target' parameter to be provided as an
%   input argument.
%
%   Example usage:
%   [Betas, BetaHistory] = findBetas(inputMats, 'Target', targetMatrix, 'Atlas', 'Schaefer', 'State', 'resting');
%
%   Author: Younghyun Oh
%   Date: 2024-07-31

    % Ensure that 'Target' parameter is provided
    if nargin < 2
        error('The ''Target'' parameter must be provided as an input argument.');
    end

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'inputMats');
    addParameter(p, 'Target', []);

    % Only add optional parameters if nargin is greater than 2
    if nargin > 2
        addParameter(p, 'Atlas', []);
        addParameter(p, 'State', []);
    end

    parse(p, inputMats, varargin{:});

    % Check if 'Target' is not provided
    if isempty(p.Results.Target)
        error('The ''Target'' parameter must be provided and cannot be empty.');
    end

    % Assign parsed values to variables
    inputMats = p.Results.inputMats;
    target = p.Results.Target;

    % Only assign atlas and state if nargin is greater than 2
    if nargin > 2
        atlas = p.Results.Atlas;
        state = p.Results.State;
    else
        atlas = '';
        state = '';
    end

    num_matrices = length(inputMats); % Get the number of input matrices

    BetaHistory = [];

    % Define the hyperparameter search space
    hyperparameters = [];
    for i = 1:num_matrices
        hyperparameters = [hyperparameters; optimizableVariable(['Beta', num2str(i)], [0, 50], 'Type', 'integer')];
    end

    % Define convergence criteria
    epsilon = 1e-100;  % Small value for convergence criterion
    bestSoFar = Inf;  % To keep track of best result so far
    count = 0;  % To keep track of how many times the change has been smaller than epsilon

    % Custom output function for bayesopt
    function stop = customOutputFcn(results, state)
        stop = false;  % Don't stop by default
        if strcmp(state, 'iteration')  % state 0 is 'iteration'
            % Check convergence criterion
            currentEstimatedBetas = table2array(results.XAtMinEstimatedObjective);
            BetaHistory = [BetaHistory; currentEstimatedBetas];

            if abs(results.MinEstimatedObjective - bestSoFar) < epsilon
                count = count + 1;  % Increment count
            else
                count = 0;  % Reset count
            end
            bestSoFar = results.MinEstimatedObjective;  % Update best so far

            % If the change has been smaller than epsilon for 5 iterations, stop
            if count >= 50
                stop = true;
            end
        end
    end

    % Bayesian optimization
    results = bayesopt(@(params) evaluateObjective(params, inputMats, target, atlas, state), hyperparameters, ...
        'MaxObjectiveEvaluations', 50, ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'OutputFcn', @customOutputFcn, ...
        'PlotFcn', []);

    % Best hyperparameters
    Betas = [];
    for i = 1:num_matrices
        Betas(i) = results.XAtMinEstimatedObjective.(['Beta', num2str(i)]);
        fprintf('Best Beta%d: %f\n', i, Betas(i));
    end
end

function performanceValue = evaluateObjective(hyperparams, input_matrices, target, atlas, state)
    % Extract matrices from input_matrices cell array
    iEC = 0;
    for i = 1:length(input_matrices)
        Beta = hyperparams.(['Beta', num2str(i)]);
        iEC = iEC + input_matrices{i} * Beta;
    end
    iEC = iEC / length(input_matrices);

    % Handle the absence of atlas and state parameters
    if isempty(atlas) && isempty(state)
        performanceValue = 1 - corr(target(:), iEC(:));
    else
        performanceValue = 1 - calculateCC(iEC, target, 1, atlas, state);
    end
end
