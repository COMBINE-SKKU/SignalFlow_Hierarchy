function [fc_sim, fcd_sim] = batch_simulation(ec, G, f_diff, TR, Tmax, num_sims)
% BATCH_SIMULATION Runs multiple Hopf model simulations and averages results
%
% Inputs:
%   ec       - Effective connectivity matrix
%   G        - Global coupling parameter
%   f_diff   - Frequency differences for each node
%   TR       - Repetition time (in seconds)
%   Tmax     - Maximum simulation time (in seconds)
%   num_sims - Number of simulations to run
%
% Outputs:
%   fc_sim  - Average functional connectivity matrix
%   fcd_sim - Average functional connectivity dynamics vector
%
% This function runs multiple simulations of the Hopf whole-brain model
% and returns the average FC and FCD measures across all simulations.

    % Get number of ROIs
    rois = size(ec, 1);
    
    % Prepare storage for FC and FCD from multiple simulations
    FC_sim_all = zeros(rois, rois, num_sims);
    FCD_all = cell(1, num_sims);
    
    % Precompute static parameters to avoid recalculation in each iteration
    scale_a = -0.01;
    
    % Run multiple simulations in parallel
    parfor sim_idx = 1:num_sims
        % Run simulation
        sim_data = run_simulation(ec, G, f_diff, TR, Tmax, scale_a)';
        
        % Compute FC - preallocate temporary arrays
        fc_temp = corr(sim_data);
        fc_temp(logical(eye(rois))) = 0; % Using logical indexing is faster
        FC_sim_all(:,:,sim_idx) = atanh(fc_temp);  % Fisher z-transform
        
        % Compute FCD
        [~, ~, fcd_hopf_vec] = calculate_fcd(sim_data, TR, ceil(60/TR), 1);
        FCD_all{sim_idx} = fcd_hopf_vec;
    end
    
    % Average FC across simulations (inverse Fisher transform)
    fc_sim = tanh(nanmean(FC_sim_all, 3));
    
    % Preallocate size for combined FCD vectors
    total_fcd_length = 0;
    for sim_idx = 1:num_sims
        total_fcd_length = total_fcd_length + numel(FCD_all{sim_idx});
    end
    hopf_FCD_concat = zeros(total_fcd_length, 1);
    
    % Combine all FCD vectors with preallocated array
    idx = 1;
    for sim_idx = 1:num_sims
        fcd_length = numel(FCD_all{sim_idx});
        hopf_FCD_concat(idx:idx+fcd_length-1) = FCD_all{sim_idx};
        idx = idx + fcd_length;
    end
    
    % Compute group CDF from the combined FCD vectors
    edges = -1:0.0002:1;
    pdf_group = histcounts(hopf_FCD_concat, edges, 'Normalization', 'pdf');
    cdf_group = cumsum(pdf_group);
    fcd_sim = cdf_group / cdf_group(end);
end 