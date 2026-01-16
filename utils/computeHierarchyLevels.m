function hierarchyLevels = computeHierarchyLevels(EC, thresholdValue)
% COMPUTEHIERARCHYLEVELS computes the hierarchy levels of a given square matrix
%   HIERARCHYLEVELS = COMPUTEHIERARCHYLEVELS(EC, THRESHOLDVALUE) returns
%   an N-by-1 vector HIERARCHYLEVELS containing the computed hierarchy levels
%   based on the input square matrix EC and a specified threshold value.
%
%   Parameters:
%   EC              A square matrix representing the connectivity data.
%   thresholdValue  A scalar value used to threshold the input matrix.
%
%   This function involves creating a binary adjacency matrix, constructing
%   a directed graph, computing an incidence matrix, solving for linear
%   regression estimates, and normalizing these estimates to determine
%   hierarchy levels.
%
%   Example usage:
%   hierarchyLevels = computeHierarchyLevels(EC, 0.5);
%
%   Author: Younghyun Oh
%   Date: 2024-07-31

    % Validate input dimensions
    if size(EC, 1) ~= size(EC, 2)
        error('EC must be a square matrix');
    end
    nNodes = size(EC,1);  
    
    % Step 1: Threshold the input matrix to construct binary adjacency matrix
    EC = threshold_absolute(EC,thresholdValue);
    Xtemp = abs(EC) > 0;

    % Step 2: Create a directed graph from the binary adjacency matrix
    G = digraph(Xtemp');

    % Step 3: Compute the incidence matrix and transpose it
    X = incidence(G)';

    % Define indices for which hierarchy levels are set to zero
    zero_indices = [1, 51];

    % Step 4: Remove specific columns to adjust degrees of freedom
    X(:, zero_indices) = [];

    % Step 5: Find indices of non-zero elements in the thresholded matrix
    [target, source] = find(Xtemp > 0);

    % Step 6: Extract corresponding values from EC using the target and source indices
    Y = arrayfun(@(i, j) EC(i, j), target, source);

    % Step 7: Solve the linear regression (OLS) equation using left division operator
    estimates = X\Y;
   
    % Step 8: Winsorize the estimates at the 99th percentile
    percentile_cutoff = 99;
    upper_limit = prctile(estimates, percentile_cutoff);
    estimates(estimates > upper_limit) = upper_limit;

    % Step 9: Normalize the estimated values to the range [0, 10]
    hier_vals = normalize(estimates, 'range', [0 10]);

    % Step 10: Initialize hierarchy levels with zeros
    hierarchyLevels = zeros(nNodes, 1);

    % Iterate through all positions, assigning computed hierarchy levels
    current_idx = 1;
    for i = 1:nNodes
        if ~ismember(i, zero_indices)
            hierarchyLevels(i) = hier_vals(current_idx);
            current_idx = current_idx + 1;
        end
    end
end
