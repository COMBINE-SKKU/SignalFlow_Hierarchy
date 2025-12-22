%% Stage 1: Run EC algorithms
% This script runs effective connectivity estimation on macaque fMRI data using
% multiple algorithms on both original and deconvolved BOLD signals.
%
% Main workflow:
% 1. Loads FLN ground truth and individual macaque BOLD signals (19 subjects)
% 2. Performs blind deconvolution using rsHRF to remove hemodynamic response
% 3. Runs 9 EC algorithms: rDCM, VAR, GC, FASK, CCD, BOSS, LiNGAM, GRASP, Patel
% 4. Saves EC matrices for each subject and algorithm for subsequent integration
%
% Author: Younghyun Oh
% Date: 2025-03-04
%
%--------------------------------------------------------------------------%
% Clear workspace
clear; close all; clc

% Load FLN
base_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy';
main_dir = fullfile(base_dir, 'step1_ec_validation/part1_macaque');
FLN_dir = fullfile(main_dir, 'data/FLN40.mat');
FLN = importdata(FLN_dir);
FLN = 1.2*FLN.^0.3; % log-linear transformation
num_regions = size(FLN,1);
num_subjects = 19;

% Load BOLD signals
BOLD_dir = fullfile(main_dir, 'data/BOLD(individual)_macaque.mat');
BOLDs = importdata(BOLD_dir);
TR = 1.6;
Tmax = size(BOLDs{1,1},1);

%--------------------------------------------------------------------------%
%% Run individual EC algorithms
%--------------------------------------------------------------------------%
% Set parameters for blind deconvolution
para.TR = TR;
para.name = 'Canonical HRF (with time and dispersion derivatives)';
temporal_mask = [];
para.T = 3;
para.T0 = 1;
para.order = 3;
para.dt = para.TR/para.T;
para.AR_lag = 1;
para.thr = 1;
para.len = 24;
min_onset_search = 4;
max_onset_search = 8;
para.lag = fix(min_onset_search/para.dt):fix(max_onset_search/para.dt);
T = round(para.len/TR);
data_deconv = zeros(Tmax,num_regions,num_subjects);

% Loop over subjects
for subj = 1:num_subjects
    data = zscore(BOLDs{subj},0,1);
    
    sigma = std(data);
    
     [beta_hrf, bf, ~] = rsHRF_estimation_temporal_basis(data,para,temporal_mask);
     hrf = bf*beta_hrf(1:size(bf,2),:);
     hrf = resample(hrf,1,para.T);
     for roi = 1:num_regions
         data_deconv(:,roi,subj) = zscore(rsHRF_iterative_wiener_deconv(data(:,roi),hrf(:,roi),80));  
     end
end

%% Initialize constants
lambda = 1e+8;
threshold = 0.1;
iter = 50;

% Initialize results storage as arrays of structures
ec_results = struct('rdcm', [], 'var', [], 'gc', [], 'fask', [], ...
    'ccd', [], 'boss', [], 'lingam', [], 'grasp', [], 'patel', []);
ec_results = repmat(ec_results, num_subjects, 1);

% Determine whether to use deconvolved data or not
use_deconv = false;

% Loop over subjects
for subj = 1:num_subjects
    if use_deconv
        BOLD = data_deconv(:,:,subj);
    else
        BOLD = BOLDs{subj};
    end
    
    % run EC algorithms
    ec_results(subj).rdcm = run_rdcm(BOLD, TR);
    ec_results(subj).var = run_var(BOLD', lambda);
    ec_results(subj).gc = run_gc(BOLD', 0, lambda);
    ec_results(subj).fask = run_fask(BOLD, 1e-6, threshold, iter);
    ec_results(subj).ccd = run_ccd(BOLD, threshold, iter);
    ec_results(subj).boss = run_boss(BOLD, threshold, iter);
    ec_results(subj).lingam = run_lingam(BOLD, threshold, iter);
    ec_results(subj).grasp = run_grasp(BOLD, threshold, iter);
    ec_results(subj).patel = run_patel(BOLD, threshold, iter);

end

% Save results
if use_deconv
    save(fullfile(main_dir, 'results/ec_results_deconv.mat'), 'ec_results');
else
    save(fullfile(main_dir, 'results/ec_results.mat'), 'ec_results');
end

%--------------------------------------------------------------------------%
%% Save group average ECs
%--------------------------------------------------------------------------%
group_mean_results = struct();

algorithms = fieldnames(ec_results);
for k = 1:length(algorithms)
    alg = algorithms{k};
    all_subjects_data = arrayfun(@(x) x.(alg), ec_results, 'UniformOutput', false);
    all_subjects_data = cat(3, all_subjects_data{:}); % Concatenate along the third dimension
    all_subjects_data = all_subjects_data/max(all_subjects_data(:));
    
    % Calculate the mean across the third dimension (subjects)
    group_mean_results.(alg) = mean(all_subjects_data, 3);
end

% Create individual connectivity matrix figures
figure('Position', [100 100 400 400]);  % Square figure for square plots

for k = 1:length(algorithms)
    alg = algorithms{k};
    
    % Clear figure for each new plot
    clf;
    
    % Plot connectivity matrix
    imagesc(group_mean_results.(alg));
    density_dir(group_mean_results.(alg))
    colormap(gca, generateColorMap(reshape(group_mean_results.(alg),1,[]),1000));
    
    % Customize plot
    axis square;
    set(gca, 'LineWidth', 1.5, 'Box', 'on', ...
        'XTick', [], 'YTick', [], 'XTickLabel', [], 'YTickLabel', []); % Remove ticks and labels
    
    % Adjust axes position to be centered in figure
    pos = get(gca, 'Position');
    pos(1) = (1-pos(3))/2;  % Center horizontally
    pos(2) = (1-pos(4))/2;  % Center vertically
    set(gca, 'Position', pos);
    
    % Save individual figure
    exportgraphics(gcf, fullfile(main_dir, 'supplementary', ['conn_matrix_', alg, '.png']), ...
        'Resolution', 300, 'BackgroundColor', 'white');
end

% Inputs: FLN_matrix, LiNGAM_matrix, FASK_matrix (all same size)

%% 1. Identify "True" Edges (where FLN is non-zero)
close all
%--------------------------------------------------------------------------%
%% Panel C: Algorithmic Contribution (Group-Level Analysis)
%--------------------------------------------------------------------------%
% Use the normalized group averages calculated above
% Note: We assume beta weights are incorporated or that we are comparing
% raw algorithmic topology. If betas are needed, multiply here.
% For the response letter regarding "Redundancy", comparing raw topological 
% output is often sufficient to show structural differences.

fask_w = abs(group_mean_results.fask);
lingam_w = abs(group_mean_results.lingam);
fln_w = abs(FLN); % Ground Truth

% 1. Setup Data (Focus only on True Biological Edges)
mask = fln_w > 0;
fln_values = fln_w(mask);
fask_vals = fask_w(mask);
lingam_vals = lingam_w(mask);

% 2. Define Linear Bins
num_bins = 20;
% Use linear spacing since FLN is already log-linear transformed
edges = linspace(min(fln_values), max(fln_values), num_bins+1);

% 3. Compute Dominance per Bin
bin_counts = zeros(num_bins, 1);
bin_winner = zeros(num_bins, 1); % 1=FASK, 2=LiNGAM

total_edges = length(fln_values);

for i = 1:num_bins
    % Find biological edges in this strength range
    idx = fln_values >= edges(i) & fln_values < edges(i+1);
    
    % Probability Density (Proportion of total connectome in this bin)
    bin_counts(i) = sum(idx) / total_edges; 
    
    if sum(idx) > 0
        % Calculate Mean Amplitude Contribution in this bin
        m_fask = mean(fask_vals(idx));
        m_lingam = mean(lingam_vals(idx));
        
        % Determine Winner: Who provides more signal for these edges?
        if m_lingam > m_fask
            bin_winner(i) = 2; % LiNGAM Wins (Orange)
        else
            bin_winner(i) = 1; % FASK Wins (Blue)
        end
    else
        bin_winner(i) = NaN; % Empty bin
    end
end

% 4. Plotting Panel C
figure('Position', [100, 100, 500, 300]);
hold on;

% Colors
c_fask = [0, 0.44, 0.74];      % Blue
c_lingam = [0.85, 0.32, 0.09]; % Orange

for i = 1:num_bins
    if ~isnan(bin_winner(i)) && bin_counts(i) > 0
        % Linear width
        w = edges(i+1) - edges(i);
        
        if bin_winner(i) == 1
            col = c_fask;
        else
            col = c_lingam;
        end
        
        rectangle('Position', [edges(i), 0, w, bin_counts(i)], ...
                  'FaceColor', col, ...
                  'EdgeColor', 'k', 'LineWidth', 0.5);
    end
end

% 5. Aesthetics
xlabel('Ground Truth Connection Strength (FLN)');
ylabel('Probability Density');
title('Algorithmic Contribution: FASK vs LiNGAM');

% Custom Legend
h1 = plot(nan, nan, 's', 'MarkerFaceColor', c_fask, 'MarkerEdgeColor', 'k');
h2 = plot(nan, nan, 's', 'MarkerFaceColor', c_lingam, 'MarkerEdgeColor', 'k');
legend([h1, h2], {'FASK Dominant', 'LiNGAM Dominant'}, 'Location', 'NorthEast');

grid on;
box on;
hold off;

% % Save the figure
% exportgraphics(gcf, fullfile(main_dir, 'supplementary', 'PanelC_Contribution.png'), ...
%     'Resolution', 300, 'BackgroundColor', 'white');

%%
%--------------------------------------------------------------------------%
%% Panel C: Algorithmic Contribution (Group Average, Consistent Style)
%--------------------------------------------------------------------------%

% 1. Load Data & Weights
fask_w = abs(group_mean_results.fask);
lingam_w = abs(group_mean_results.lingam);
fln_w = abs(FLN); % Ground Truth

% Setup Data (Focus only on True Biological Edges)
mask = fln_w > 0;
fln_values = fln_w(mask);
fask_vals = fask_w(mask);
lingam_vals = lingam_w(mask);

% 2. Define Linear Bins
num_bins = 100;
% Use linear spacing since FLN is already log-linear transformed
edges = linspace(min(fln_values), max(fln_values), num_bins+1);

% 3. Compute Dominance per Bin
bin_counts = zeros(num_bins, 1);
bin_winner = zeros(num_bins, 1); % 1=FASK, 2=LiNGAM

total_edges = length(fln_values);

for i = 1:num_bins
    % Find biological edges in this strength range
    idx = fln_values >= edges(i) & fln_values < edges(i+1);
    
    % Probability Density
    bin_counts(i) = sum(idx) / total_edges; 
    
    if sum(idx) > 0
        % Calculate Mean Amplitude Contribution
        m_fask = mean(fask_vals(idx));
        m_lingam = mean(lingam_vals(idx));
        
        % Determine Winner
        if m_lingam > m_fask
            bin_winner(i) = 2; % LiNGAM Wins
        else
            bin_winner(i) = 1; % FASK Wins
        end
    else
        bin_winner(i) = NaN;
    end
end

% 4. Create Figure (Matching LOO Figure Dimensions & Style)
% Note: Width increased to 600 to accommodate the distribution
figure('Position', [100 100 600 350]); 
hold on;

% Define Colors (Matched to LOO Box Plot)
% "No LiNGAM" (Blue) implies FASK is the remaining driver -> FASK Dominant = Blue
color_fask = [0.6 0.8 0.9];   % Pastel Blue
% "No FASK" (Orange) implies LiNGAM is the remaining driver -> LiNGAM Dominant = Orange
color_lingam = [0.9 0.8 0.6]; % Pastel Orange

% Plot Bars
for i = 1:num_bins
    if ~isnan(bin_winner(i)) && bin_counts(i) > 0
        w = edges(i+1) - edges(i);
        
        if bin_winner(i) == 1
            col = color_fask;
        else
            col = color_lingam;
        end
        
        % Draw rectangle with box plot styling
        rectangle('Position', [edges(i), 0, w, bin_counts(i)], ...
                  'FaceColor', col, ...
                  'EdgeColor', 'k', ...
                  'LineWidth', 0.5);
    end
end

% 5. Customize Plot Appearance
box on
grid on

xlabel('Ground Truth Connection Strength (FLN)', 'FontSize', 12);
ylabel('Probability Density', 'FontSize', 12);
title('Algorithmic Contribution: FASK vs LiNGAM', 'FontSize', 12, 'FontWeight', 'normal');

% Set limits to match data range tightly
xlim([min(edges) max(edges)]);
ylim([0 max(bin_counts)*1.1]);

% Custom Legend (Invisible markers for legend generation)
h1 = plot(nan, nan, 's', 'MarkerFaceColor', color_fask, 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
h2 = plot(nan, nan, 's', 'MarkerFaceColor', color_lingam, 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
legend([h1, h2], {'FASK Dominant', 'LiNGAM Dominant'}, 'Location', 'NorthEast', 'FontSize', 10);

hold off;

% Save the figure
% exportgraphics(gcf, fullfile(main_dir, 'results/PanelC_Contribution_Final.png'), 'Resolution', 1200);