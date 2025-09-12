function [Betas, G_opt, w_opt, gof_opt] = findBetas_new(inputMats, f_diff, FC_emp_group, TR, Tmax, ...
    num_sims, triu_ind_fc, FCD_emp_group)
% FINDBETAS  Estimate relative ECmatrix weights and global scale G
%   Implements a two-stage parameterization:
%     w  = simplex (relative weights),
%     G = [0,10]       (absolute scale).
%   Also seeds RNG deterministically from the current params so that
%   repeated evaluations at the same params produce the same noise.
%
% INPUTS
%   inputMats       Cell array of NxN effective-connectivity matrices
%   f_diff, FC_emp_group, TR, Tmax, num_sims, triu_ind_fc, ...
%   FCD_emp_group
% OUTPUTS
%   Betas           1xM vector of Gxw_i
%   metrics         Struct with fields:
%                     .gof       = -MinObjective
%                     .fc_corr   = FC correlation at optimum
%                     .fcd_ks    = FCD KS statistic at optimum
%   G_opt            Global scale parameter at optimum
%   w_opt            1M vector of weights at optimum

M = numel(inputMats);

% Define hyperparameters for Bayesian Optimization
vars = optimizableVariable.empty();
for i = 1:M
    vars(end+1) = optimizableVariable(sprintf('x%d', i), [-5, 5], 'Type', 'real');
end
vars(end+1) = optimizableVariable('G', [0, 20], 'Type', 'real');

% Nested objective function
function neg_gof = evalObj(params)
    %----  Deterministic RNG seed from params ----
    pv = zeros(M+1, 1);
    for j = 1:M
        pv(j) = params.(['x', num2str(j)]);
    end
    pv(end) = params.G;
    str = sprintf('%.8f,', pv);           % string repr of vector
    seed = sum(double(str));             % simple hash
    rng(mod(seed, 2^32), 'twister');      % set RNG
    
    %---- 2b) Softmax – weights w ----
    xs = zeros(M, 1);
    for j = 1:M
        xs(j) = params.(['x', num2str(j)]);
    end
    ex = exp(xs);
    w  = ex ./ sum(ex);                   % w sums to 1
    
    %---- 2c) Build iEC = G * sum_i w(i)*inputMats{i} ----
    G_local = params.G;
    iEC = zeros(size(inputMats{1}));
    for j = 1:M
        iEC = iEC + w(j) * G_local * inputMats{j};
    end
    
    %---- 2d) Simulation & metric computation ----
 
    test_data = run_simulation(iEC, 1, f_diff, TR, Tmax, -0.01)';
    if any(isnan(test_data(:))) || any(isinf(test_data(:)))
        neg_gof = Inf; 
        return;
    end
    [FC_sim, FCD_sim] = batch_simulation(iEC, 1, f_diff, TR, Tmax, num_sims);
    if any(isnan(FC_sim(:))) || any(isinf(FC_sim(:)))
        neg_gof = Inf; 
        return;
    end
    
    %---- 2e) Extract upper-triangular FC and compute stats ----
    FC_sim_triu = FC_sim(triu_ind_fc);
    FC_emp_triu = FC_emp_group(triu_ind_fc);
    fc_corr = corr(FC_emp_triu, FC_sim_triu, 'Rows', 'complete');
    fcd_ks = max(abs(FCD_sim - FCD_emp_group));
    neg_gof = (1-fc_corr) + fcd_ks;
end

% Run Bayesian optimization
results = bayesopt(@evalObj, vars, ...
'MaxObjectiveEvaluations', 200, ...
'NumSeedPoints', 80, ...
'IsObjectiveDeterministic', false, ...
'AcquisitionFunctionName','expected-improvement-plus', ...
'PlotFcn', []);

% Extract optimal weights and global scale parameter
opt = results.XAtMinObjective;
xs_opt = zeros(M, 1);
for j = 1:M
    xs_opt(j) = opt.(['x', num2str(j)]);
end
ex_opt = exp(xs_opt);
w_opt = ex_opt ./ sum(ex_opt);
G_opt = opt.G;

% Extract best gof
gof_opt = -results.MinObjective;

Betas = (G_opt * w_opt)';   % 1 x M vector
end
