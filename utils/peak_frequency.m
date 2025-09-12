%% calculate peak frequency
% This script calculates the peak frequency of BOLD signals from macaque data.
% It loads the BOLD data, processes it to compute power spectra, and identifies
% the peak frequency for each area. The results are saved to a .mat file.

clear; close all; clc

% Define directories and file paths
data_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy/data';
bold_file = fullfile(data_dir, 'BOLD(individual)_macaque_individual.mat');
save_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy/data';

% Load BOLD data
BOLDs = importdata(bold_file);
BOLD3d = zeros(1200, 100, 400);

% Reshape BOLD data into 3D matrix
valid_subjects = zeros(1, length(BOLDs));
subject_counter = 0;

for i = 1:length(BOLDs)
    current_data = BOLDs{i}(:,33:end);
    [rows, cols] = size(current_data);
    
    if rows == 1200 && cols == 100
        subject_counter = subject_counter + 1;
        BOLD3d(:, :, subject_counter) = current_data;
        valid_subjects(i) = 1;
    end
end

% Trim BOLD3d to actual size and update subjects count
BOLD3d = BOLD3d(:, :, 1:subject_counter);

% Set parameters
[Tmax, N, subjects] = size(BOLD3d);
TR = 0.72;
Ts = Tmax * TR;
freq = (0:Tmax/2-1) / Ts;
nfreqs = length(freq);
PowSpect = zeros(nfreqs, N, subjects);

% Calculate power spectrum for each subject and seed
for sub = 1:subjects
    data = squeeze(BOLD3d(:, :, sub));
    parfor seed = 1:N
%         bpdata = Butterworth_bandpass_matrix(data(:, seed), 0.04, 0.07, 0.72);
        pw = abs(fft(data(:, seed)));
        PowSpect(:, seed, sub) = pw(1:floor(Tmax/2)).^2 / (Tmax / TR);
    end
end

% Compute mean power across subjects and find peak frequency
Power_Areas = mean(PowSpect, 3);
[~, index] = max(Power_Areas);
f_diff = freq(index);

% Save peak frequency data
save(fullfile(save_dir, 'Schaefer100_peak_freq.mat'), 'f_diff')