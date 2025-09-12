 function [best_G, gof] = findOptimalParams(ec_alg, f_diff, FC_emp_group, TR, Tmax, ...
    num_sims, triu_ind_fc, FCD_emp_group, optimizeFCD, emp_mean_lag, optimizeLag, G_range)
% FINOPTIMALPARAMS  Find optimal global coupling G (linear range)
%
%   [best_G, metrics] = findOptimalParams(ec_alg, f_diff, FC_emp_group, TR, ...
%                 Tmax, num_sims, triu_ind_fc, FCD_emp_group, ...
%                 optimizeFCD, emp_mean_lag, optimizeLag, G_range)
%
%   - ec_alg         : base EC matrix
%   - f_diff, TR, Tmax, num_sims, triu_ind_fc : as before
%   - FCD_emp_group  : empirical FCD CDF (vector)
%   - optimizeFCD    : boolean, include FCD in loss?
%   - emp_mean_lag   : empirical mean-lag vector
%   - optimizeLag    : boolean, optimize lag instead of FC/FCD?
%   - G_range        : two-element vector [Gmin, Gmax] (default [0.1,10])
%
% OUTPUTS:
%   best_G          : optimal global coupling G
%   metrics         : struct with fields
%                       .gof       = -MinObjective
%                       .fc_corr   = FC correlation at optimum
%                       .fcd_ks    = FCD KS statistic (or NaN)
%                       .lag_corr  = lag correlation (or NaN)

    if nargin < 9
        optimizeFCD = false;
    end
    
    if nargin < 11
        optimizeLag = false;
    end
    
    if nargin < 12 || isempty(G_range)
        G_range = [0.5, 5];
    end

    % Precompute constant vectors/matrices
    FC_emp_triu = FC_emp_group(triu_ind_fc);

    % Define G as a linear hyperparameter
    vars = optimizableVariable('G', G_range, 'Type', 'real');

    % Nested objective function
    function neg_gof = objFun(params)
        G = params.G;
        
        % Set a deterministic RNG seed from G
        s = sprintf('%.8f', G);
        seed = sum(double(s));
        rng(mod(seed,2^32), 'combRecursive');
        
        % If optimizing lag correlation
        if optimizeLag && ~isempty(emp_mean_lag)
            ts0 = run_simulation(ec_alg, G, f_diff, TR, Tmax, -0.01)';
            if any(isnan(ts0(:))) || any(isinf(ts0(:)))
                neg_gof = Inf;
                return;
            end
            
            rois = size(ec_alg, 1);
            td_accum = zeros(rois, rois, num_sims);
            parfor k = 1:num_sims
                ts = run_simulation(ec_alg, G, f_diff, TR, Tmax, -0.01)';
                td_accum(:,:,k) = calculate_lag_threads(ts, TR);
            end
            mean_td_sim = mean(td_accum, 3);
            sim_mean_lag = mean(mean_td_sim, 1);
            lag_corr = corr(emp_mean_lag', sim_mean_lag');
            neg_gof = -lag_corr;
            return;
        end
        
        % Else compute FC (and optionally FCD)
        ts0 = run_simulation(ec_alg, G, f_diff, TR, Tmax, -0.01)';
        if any(isnan(ts0(:))) || any(isinf(ts0(:)))
            neg_gof = Inf;
            return;
        end
        
        % Always use batch simulation since num_sims is never 1
        [FC_sim, FCD_sim] = batch_simulation(ec_alg, G, f_diff, TR, Tmax, num_sims);
        
        FC_sim_triu = FC_sim(triu_ind_fc);
        fc_corr = corr(FC_emp_triu, FC_sim_triu);
        ks = max(abs(FCD_sim - FCD_emp_group));
        neg_gof = (1-fc_corr) + ks;

    end

    % Run Bayesian optimization
    results = bayesopt(@objFun, vars, ...
        'MaxObjectiveEvaluations', 30, ...
        'IsObjectiveDeterministic', false, ...
        'UseParallel', false, ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'PlotFcn', []);

    % Extract optimum and recompute metrics
    best_G = results.XAtMinObjective.G;
    gof = 1-results.MinObjective;
end
