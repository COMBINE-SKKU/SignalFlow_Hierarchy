%% Brain Network Dynamics Simulation
% This script simulates brain dynamics with different coupling scenarios and parameters.
% It demonstrates how different bifurcation parameters affect node dynamics in
% uncoupled scenarios.
%
% Key parameters:
%   - TR: Repetition time (sampling interval) in seconds
%   - a_values: Bifurcation parameter values controlling node dynamics
%
% Author: [Your Name]
% Date: [Current Date]

%% Initialize parameters
TR = 0.72;               % Repetition time in seconds
Tmax = 2400;             % Maximum simulation time (number of samples)
N = 360;                 % Number of brain regions/nodes
a_values = [-1, -0.01, 1]; % Bifurcation parameter values to test
labels = {'a = -1', 'a = -0.01', 'a = +1'};

% Create output directory for figures if it doesn't exist
figDir = 'figures';
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

%% Simulation: Uncoupled nodes with constant frequency
% f_diff should be defined elsewhere in the workspace
iEC_zero = zeros(N);     % No coupling (G = 0)

bold_all = cell(1,3);

% Simulate with zero EC to illustrate dynamics from local oscillator only
fprintf('Running simulations for uncoupled nodes with different bifurcation parameters...\n');
for i = 1:3
    a_val = a_values(i);
    fprintf('  Simulating uncoupled nodes with a = %.2f\n', a_val);
    bold_sim = run_simulation(iEC_zero, 0, f_diff, TR, Tmax, a_val)';  % Tmax x N
    bold_all{i} = bold_sim;
end

%% Plot single region time series from uncoupled simulations
roi_id = 1;
tvec = (0:Tmax-1) * TR;
fig1 = figure('Position', [100, 100, 800, 600]);

% Set default font to Arial
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

for i = 1:3
    subplot(3,1,i);
    
    % Create publication-quality plot
    plot(tvec, bold_all{i}(:,roi_id), 'LineWidth', 2, 'Color', [0.2 0.5 0.9]);
    
    % Minimalist design - remove all labels and ticks
    % set(gca, 'XTickLabel', [], 'YTickLabel', [], 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on', 'FontName', 'Arial');
    
    % Keep only the plot display range
    xlim([0, 500]);
    grid off;
end

% Make figure clean with white background
set(gcf, 'Color', 'w');
tight_spacing = 0.03;
set(fig1, 'Units', 'normalized');

% Uncomment to save the figure
% exportgraphics(fig1, fullfile(figDir, 'uncoupled_dynamics.png'), 'Resolution', 300);
% exportgraphics(fig1, fullfile(figDir, 'uncoupled_dynamics.pdf'), 'ContentType', 'vector');

fprintf('Analysis complete!\n');
