function TD = calculate_lag_threads(BOLD, TR)
    % Calculate lag threads from BOLD time series
    %
    % Inputs:
    % BOLD - Time series matrix (T timepoints × N regions)
    % TR - Repetition time in seconds
    %
    % Outputs:
    % TD - Time delay matrix (N × N)

    % Get dimensions
    [~, N] = size(BOLD);
    
    % Maximum lag to consider (at least 2 seconds converted to TR units)
    max_lag = max(3, round(2/TR));
    
    % Pre-allocate temporary array for parallel processing
    TD_upper = zeros(N, N);
    
    % Use parallel processing for upper triangle
    parfor i = 1:N
        TD_row = zeros(1, N);
        for j = (i+1):N
            % Calculate cross-covariance
            [c, lags] = xcov(BOLD(:,i), BOLD(:,j), max_lag, 'coeff');
            
            % Find lag with maximum absolute correlation
            [~, idx] = max(abs(c));
            TD_row(j) = lags(idx) * TR;  % Convert to seconds
        end
        TD_upper(i,:) = TD_row;
    end
    
    % Fill in the full matrix using anti-symmetry
    TD = TD_upper - TD_upper';  % This automatically fills both upper and lower triangles
end