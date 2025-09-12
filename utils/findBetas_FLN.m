function [Betas, corr_opt, w_opt] = findBetas_FLN(inputMats, target)
% FINDBETAS_FLN finds optimal beta values to integrate algorithms
%   [BETAS, CORR_OPT, W_OPT] = FINDBETAS_FLN(INPUTMATS, TARGET) returns 
%   an N-by-1 vector BETAS containing weights that maximize the correlation 
%   between the integrated EC and TARGET.
%
%   INPUTMATS   Cell array of connectivity matrices to be integrated
%   TARGET      Ground truth (or proxy) matrix to compute Pearson correlation
%
% OUTPUTS
%   Betas       1xM vector of weights
%   corr_opt    Correlation at the optimum
%   w_opt       1xM vector of optimal weights
%
%   Example usage:
%       [Betas, corr_opt, w_opt] = findBetas_FLN(inputMats, targetMatrix);

% Input validation
if isempty(target)
    error('The target matrix cannot be empty.');
end

M = length(inputMats);

% Define hyperparameters for Bayesian Optimization
vars = optimizableVariable.empty();
for i = 1:M
    vars(end+1) = optimizableVariable(sprintf('x%d', i), [-5, 5], 'Type', 'real');
end

% Nested objective function
function obj_val = evalObj(params)
    %---- Extract parameters ----
    xs = zeros(M, 1);
    for j = 1:M
        xs(j) = params.(['x', num2str(j)]);
    end
    
    %---- Softmax – weights w ----
    ex = exp(xs);
    w = ex ./ sum(ex);    % w sums to 1
    
    %---- Build integrated EC = sum_i w(i)*inputMats{i} ----
    iEC = zeros(size(inputMats{1}));
    for j = 1:M
        iEC = iEC + w(j) * inputMats{j};
    end
    
    %---- Compute correlation-based objective ----
    rVal = corr(iEC(:), target(:));
    obj_val = 1 - rVal;  % Minimizing => maximizing correlation
end

% Run Bayesian optimization
results = bayesopt(@evalObj, vars, ...
    'MaxObjectiveEvaluations', 100, ...
    'NumSeedPoints', 80, ...
    'IsObjectiveDeterministic', true, ...
    'PlotFcn', []);

% Extract optimal weights
opt = results.XAtMinObjective;
xs_opt = zeros(M, 1);
for j = 1:M
    xs_opt(j) = opt.(['x', num2str(j)]);
end
ex_opt = exp(xs_opt);
w_opt = ex_opt ./ sum(ex_opt);

% Extract best correlation
corr_opt = 1 - results.MinObjective;

% Return the optimal weights
Betas = w_opt';   % 1 x M vector

% Print final Betas
fprintf('\nOptimal Beta Values:\n');
for i = 1:M
    fprintf('   Beta%d = %.6f\n', i, Betas(i));
end

end 
