function [bestO_Alpha, bestAlpha, bestThreshold] = FASK_optimization2(BOLD, atlas_name, state)
    % Define search space for hyperparameters
    alphas = [1e-7, 1e-6, 1e-5];  % The alpha search space (R^3)
    oAlphas = [1e-7, 1e-6, 1e-5]; % The oAlpha search space (R^3)
    thresholds = [0.05, 0.1, 0.15];  % The threshold search space (R^3)
    
    % Initialize the best performance variable and corresponding hyperparameters
    bestPerformance = Inf;
    bestAlpha = NaN;
    bestO_Alpha = NaN;
    bestThreshold = NaN;
    
    % Compute functional connectivity (FC)
    FC = corr(BOLD);
    FC = FC - diag(diag(FC));
    
    % Loop over the entire search space to find the optimal hyperparameters (O(27) complexity)
    for alpha = alphas
        for oAlpha = oAlphas
            for threshold = thresholds
                % Generate output matrix using the current set of hyperparameters and data
                FASK = subsample_wrapper(BOLD, oAlpha, alpha, threshold, 50);  % Assuming 100 for bootstrap sampling

                performanceValue = 1- calculateCC(FASK,FC,3,atlas_name,state);
                
                % Update the best hyperparameters if the current one is better
                if performanceValue < bestPerformance
                    bestPerformance = performanceValue;
                    bestAlpha = alpha;
                    bestO_Alpha = oAlpha;
                    bestThreshold = threshold;
                end
            end
        end
    end
    
    % Display the best hyperparameters
    fprintf('Best Alpha: %e\n', bestAlpha);
    fprintf('Best o-Alpha: %e\n', bestO_Alpha);
    fprintf('Best Threshold: %f\n', bestThreshold);
end