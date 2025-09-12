function output_mtx = Butterworth_bandpass_matrix(input_mtx, low_f, high_f, TR)
% output_mtx = CBIG_Butterworth_bandpass_matrix(input_mtx, low_f, high_f, TR)
%
% This function applies a Butterworth bandpass filter on the input data matrix.
%
% Inputs:
%   - input_mtx:
%     A T x N matrix, where T is the number of time points, and N is the number of samples.
%     input_mtx(t, n) represents the signal intensity of the n-th sample at time t.
%
%   - low_f:
%     A scalar representing the low cutoff frequency (Hz), e.g., 0.01.
%
%   - high_f:
%     A scalar representing the high cutoff frequency (Hz), e.g., 0.08.
%
%   - TR:
%     A scalar representing the sampling period (seconds), e.g., 2.
%
% Outputs:
%   - output_mtx:
%     A T x N matrix, where T is the number of time points, and N is the number of samples.
%     output_mtx(t, n) represents the filtered signal intensity of the n-th sample at time t.
%
% Example:
% output_mtx = CBIG_Butterworth_bandpass_matrix(input_mtx, 0.01, 0.08, 2);

% Check if input_matrix is valid without nan or inf
if ~any(isnan(input_mtx(:))) && ~any(isinf(input_mtx(:)))
    
    % Check the input frequencies
    if (low_f >= high_f)
        error('ERROR: low_f must be less than high_f.');
    end

    % Calculate the Nyquist frequency
    Nyquist = 1/(2*TR);

    % Check if the high cutoff frequency is above the Nyquist frequency
    if high_f > Nyquist
        error('ERROR: high_f must be less than or equal to the Nyquist frequency.');
    end

    % Design Butterworth filter
    Wn = [low_f/Nyquist high_f/Nyquist]; % Set default bandpass frequencies
    k = 2;
    [bfilt, afilt] = butter(k, Wn, 'bandpass');

    % Apply the filter
    output_mtx = zeros(size(input_mtx)); % Preallocate output matrix
    for i = 1:size(input_mtx, 2) % Loop through each sample
        % Preprocess the signal
        ts = zscore(input_mtx(:, i)); % z-score
        output_mtx(:, i) = filtfilt(bfilt, afilt, ts); % Zero-phase filtering
    end

else
    output_mtx = input_mtx;
end
