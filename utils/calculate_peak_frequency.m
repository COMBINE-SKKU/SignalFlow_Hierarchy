function peak_freqs = calculate_peak_frequency(bold_data, TR)
% CALCULATE_PEAK_FREQUENCY Calculate peak frequency for each brain region
%
% Inputs:
%   bold_data: Time series x Regions matrix of BOLD signals
%   TR: Repetition time (sampling interval) in seconds
%
% Output:
%   peak_freqs: Vector containing peak frequency for each brain region

% Get dimensions
[Tmax, N] = size(bold_data);
Ts = Tmax * TR;
freq = (0:Tmax/2-1) / Ts;
nfreqs = length(freq);

% Initialize power spectrum matrix
PowSpect = zeros(nfreqs, N);

% Calculate power spectrum for each region
for seed = 1:N
    pw = abs(fft(bold_data(:, seed)));
    PowSpect(:, seed) = pw(1:floor(Tmax/2)).^2 / (Tmax / TR);
end

% Find peak frequency for each region
[~, index] = max(PowSpect);
peak_freqs = freq(index);

end 