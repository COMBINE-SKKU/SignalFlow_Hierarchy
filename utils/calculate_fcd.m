function [FCD_CDF, FCD_matrix, FCD_vec] = calculate_fcd(TS, TR, window_size, window_step, method)
    % calculate_fcd - Computes the Functional Connectivity Dynamics (FCD) from BOLD time series data.
    %
    % Syntax: [FCD_CDF, FCD_matrix, FCD_PDF] = calculate_fcd(TS, TR, window_size, window_step, method)
    %
    % Inputs:
    %   TS - Time series data (time points x regions)
    %   TR - Repetition time in seconds
    %   window_size - Size of the window in time points (optional)
    %   window_step - Step size for moving the window (optional)
    %   method - Method to compute FCD: 'sliding_window' (default) or 'hilbert'
    %
    % Outputs:
    %   FCD_CDF - Cumulative distribution function of the FCD vector
    %   FCD_matrix - Full FCD correlation matrix
    %   FCD_PDF - Probability density function of the FCD vector
    
    % Set default method if not provided
    if nargin < 5 || isempty(method)
        method = 'sliding_window';
    end
    
    % Check if method is valid
    valid_methods = {'sliding_window', 'hilbert'};
    if ~ismember(method, valid_methods)
        error('Invalid method. Choose either ''sliding_window'' or ''hilbert''');
    end
    
    % Process based on selected method
    if strcmp(method, 'sliding_window')
        % Sliding window method (original implementation)
        if nargin < 3 || isempty(window_size)
            % Calculate default window size for 60 seconds
            window_size = round(60 / TR);
        end
        
        if nargin < 4 || isempty(window_step)
            % Default window step is 1
            window_step = 1;
        end

        [num_timepoints, num_roi] = size(TS);

        % Initialize FCD matrix
        FCD_matrix = zeros(num_roi * (num_roi - 1) / 2, floor((num_timepoints - window_size) / window_step) + 1);

        % Compute FC for each window and extract lower triangular vector
        % Use parfor to parallelize this loop
        parfor idx = 1:floor((num_timepoints - window_size) / window_step) + 1
            j = (idx - 1) * window_step + 1;
            TC_section = TS(j:j + window_size - 1, :);
            FC_section = corr(TC_section);
            FC_vec_section = FC_section(triu(true(size(FC_section, 1)), 1));
            FCD_matrix(:, idx) = FC_vec_section;
        end

        % Compute correlation on the appended matrix
        FCD_matrix = corr(FCD_matrix);

        % Extract the upper triangle of FCD matrix (excluding the main diagonal)
        FCD_vec = FCD_matrix(triu(true(size(FCD_matrix, 1)), 1));

        edges = -1:0.0002:1;
        FCD_PDF = histcounts(FCD_vec, edges, 'Normalization', 'pdf');

        % then turn it into a CDF
        FCD_CDF = cumsum(FCD_PDF);
        
        % Normalize the CDF
        FCD_CDF = FCD_CDF / FCD_CDF(end);
        
    else
        % Hilbert transform method (from computeFCD)
        % Step 1: Compute the Hilbert transform
        hilbert_transformed = hilbert(TS); % Hilbert transform (complex signal)
        
        % Step 2: Extract instantaneous complex argument (theta)
        theta = angle(hilbert_transformed); % Phase of the analytic signal
        
        % Step 3: Compute synchrony (?(i,j,t)) - Vectorized version
        Tmax = size(TS, 1);
        Nregions = size(TS, 2);
        
        % Reshape theta for broadcasting
        theta_i = reshape(theta, [Tmax, Nregions, 1]);
        theta_j = reshape(theta, [Tmax, 1, Nregions]);
        
        % Compute all phase differences at once
        delta_ijt = cos(theta_i - theta_j);
        delta_ijt = permute(delta_ijt, [2 3 1]); % Reshape to [Nregions x Nregions x Tmax]
        
        % Step 4: Compute global synchrony ?(?_u, ?_v) - Vectorized version
        % Precompute the squared norms
        dx = sqrt(sum(sum(delta_ijt.^2, 1), 2)); % [1 x 1 x Tmax]
        dx = squeeze(dx); % Convert to [Tmax x 1]
        
        % Reshape delta_ijt for matrix multiplication
        delta_reshaped = reshape(delta_ijt, [Nregions^2, Tmax]);
        
        % Compute all numerators at once using matrix multiplication
        numerators = delta_reshaped' * delta_reshaped; % [Tmax x Tmax]
        
        % Compute denominators using outer product
        denominators = dx * dx';
        
        % Final FCD computation
        FCD_matrix = numerators ./ denominators;
        
        % Extract the upper triangle of FCD matrix (excluding the main diagonal)
        FCD_vec = FCD_matrix(triu(true(size(FCD_matrix, 1)), 1));
        
        % turn it into a PDF, 10000 bins 
        FCD_PDF = histcounts(sort(FCD_vec), -1:0.0002:1);

        % then turn it into a CDF
        FCD_CDF = cumsum(FCD_PDF);
    end
end
