function [Betas, gof] = findBetas(inputMats, f_diff, FC_emp_group, TR, Tmax, ...
    num_sims, triu_ind_fc, FCD_emp_group, use_fcd_ks, range_betas)

    % Number of matrices
    num_matrices = length(inputMats);

    % Define optimizable hyperparameters (adjust Type as needed)
    hyperparameters = optimizableVariable.empty();
    for i = 1:num_matrices
        hyperparameters(end+1) = optimizableVariable(['Beta', num2str(i)], ...
            range_betas, 'Type', 'real');
        % If you want real-valued betas: 'Type','real'
    end

    %--- Convergence criteria
    epsilon = 1e-5;      % Tolerance for small improvement
    bestSoFar = Inf;     % Track the best (lowest) objective so far
    count = 0;           % Count how many consecutive times improvement < epsilon
    consecutiveItersNeeded = 30;  % Stop if improvement < epsilon for 5 iterations

    function stop = customOutputFcn(results, state)
        stop = false;
        if strcmp(state, 'iteration')
            currentObj = results.MinObjective;  % The best (lowest) observed so far

            if abs(currentObj - bestSoFar) < epsilon
                count = count + 1;
            else
                count = 0;
            end
            bestSoFar = currentObj;

            % If we've had 'consecutiveItersNeeded' iterations of tiny improvement, stop
            if count >= consecutiveItersNeeded
                stop = true;
            end
        end
    end

    %--- Nested objective wrapper so we can pass extra args cleanly
    function neg_gof = evaluateObjectiveNested(params)
        neg_gof = evaluateObjective(params, inputMats, f_diff, FC_emp_group, ...
            FCD_emp_group, TR, Tmax, num_sims, triu_ind_fc, use_fcd_ks);
    end

    %--- Run Bayesian optimization
    results = bayesopt(@evaluateObjectiveNested, hyperparameters, ...
        'MaxObjectiveEvaluations', 500, ...
        'OutputFcn', @customOutputFcn, ...
        'Verbose', 1, ...
        'PlotFcn', []);

    %--- Extract best hyperparameters and GOF
    Betas = zeros(1, num_matrices);
    for i = 1:num_matrices
        Betas(i) = results.XAtMinObjective.(['Beta', num2str(i)]);
        fprintf('Best Beta%d: %g\n', i, Betas(i));
    end
    gof = -results.MinObjective;

    %--- Recompute gof using optimal Betas
    % Construct iEC using optimal Betas
    iEC_optimal = 0;
    for i = 1:num_matrices
        iEC_optimal = iEC_optimal + inputMats{i} * Betas(i);
    end
    iEC_optimal = iEC_optimal / num_matrices;

    % Run batch simulation with optimal parameters
    [FC_sim_group_optimal, FCD_sim_group_optimal] = batch_simulation(iEC_optimal, 1, f_diff, TR, Tmax, num_sims);

    % Extract upper triangular elements for comparison
    FC_sim_triu = FC_sim_group_optimal(triu_ind_fc);
    FC_emp_triu = FC_emp_group(triu_ind_fc);
    
    % Correlation between upper-triangular elements
    fc_corr_optimal = corr(FC_emp_triu, FC_sim_triu, 'Rows', 'complete');

    cdf1 = FCD_sim_group_optimal;
    cdf2 = FCD_emp_group;
    fcd_ks_optimal = max(abs(cdf1 - cdf2));

    % Display the results
    fprintf('FC Correlation: %g\n', fc_corr_optimal);
    fprintf('FCD KS Statistic: %g\n', fcd_ks_optimal);

end

%######################################################################
function neg_gof = evaluateObjective(params, input_matrices, f_diff, ...
    FC_emp_group, FCD_emp_group, TR, Tmax, num_sims, triu_ind_fc, use_fcd_ks)

    % Combine connectivity matrices weighted by each Beta
    num_matrices = length(input_matrices);
    iEC = 0;
    for i = 1:num_matrices
        Beta = params.(['Beta', num2str(i)]);
        iEC = iEC + input_matrices{i} * Beta;
    end
    iEC = iEC / num_matrices;  % Average if that's your intended formula

    % Handle single simulation case separately
    if num_sims == 1
        % Quick stability check via single simulation
        test_data = run_simulation(iEC, 1, f_diff, TR, Tmax, -0.01)';
        if any(isnan(test_data(:))) || any(isinf(test_data(:)))
            neg_gof = Inf;
            return;
        end
        
        % Compute FC
        fc_temp = corr(test_data);
        fc_temp(eye(size(fc_temp))==1) = 0;  % optional removal of diagonal
        FC_sim_group = fc_temp;
        
        % Compute FCD if needed
        if use_fcd_ks && ~isempty(FCD_emp_group)
            [FCD_sim_group, ~] = calculate_fcd(test_data, TR, 83, 1);
        end
    else
        % First run a single simulation to check for NaN/Inf values
        test_data = run_simulation(iEC, 1, f_diff, TR, Tmax, -0.01)';
        if any(isnan(test_data(:))) || any(isinf(test_data(:)))
            neg_gof = Inf;
            return;
        end
        
        % If the test simulation is valid, use batch_simulation for multiple simulations
        [FC_sim_group, FCD_sim_group] = batch_simulation(iEC, 1, f_diff, TR, Tmax, num_sims);
        
        % Additional check for invalid FC (just in case)
        if any(isnan(FC_sim_group(:))) || any(isinf(FC_sim_group(:)))
            neg_gof = Inf;
            return;
        end
    end

    % Correlation between upper-triangular elements
    FC_sim_triu = FC_sim_group(triu_ind_fc);
    FC_emp_triu = FC_emp_group(triu_ind_fc);
    fc_corr = corr(FC_emp_triu, FC_sim_triu, 'Rows', 'complete');
    
    % Compute absolute difference of mean values
    mean_diff = abs(mean(FC_sim_triu) - mean(FC_emp_triu));
    
    % Compute normalized Euclidean distance
    euclidean = sqrt(sum((FC_sim_triu - FC_emp_triu).^2));
    max_euclidean = 2 * sqrt(length(FC_sim_triu)); % Maximum possible distance
    euclidean_norm = euclidean / max_euclidean; % Normalized to 0-1 range
    
    % Define weights for metric components
    w_fc = 1.0;   % Weight for FC correlation
    w_fcd = 1.0;  % Weight for FCD KS statistic
    w_mean = 1.0; % Weight for mean difference
    w_euc = 1.0;  % Weight for Euclidean distance

    if use_fcd_ks && ~isempty(FCD_emp_group)
        % KS statistic for difference in distribution
        cdf1 = FCD_sim_group;
        cdf2 = FCD_emp_group;
        fcd_ks = max(abs(cdf1 - cdf2));

        % Negative objective: we want to maximize fc_corr, minimize fcd_ks, mean_diff, and euclidean_norm
        neg_gof = -w_fc * fc_corr + w_fcd * fcd_ks;
    else
        neg_gof = -w_fc * fc_corr + w_mean * mean_diff + w_euc * euclidean_norm;
    end
end