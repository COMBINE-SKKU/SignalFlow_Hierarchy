%% supplementary_analysis.m
% This script induces directionality into structural connectivity matrices
% and performs topological analyses to compare original and directed networks.
% 
% Main functions:
% 1. Loads SC matrices for Schaefer100 and MMP360 atlases
% 2. Applies randmio_dir_connected to induce directionality (single iteration)
% 3. Compares topological properties (clustering coefficient, degree)
% 4. Creates visualizations of connectivity matrices and scatter plots
% 5. Calculates and reports correlations between original and directed networks
%
% Author: Younghyun Oh
% Date: 2025-09-12
% Version: 1.0

%% Initialize environment
% Clear workspace
clear; close all; clc

% Load structural connectivity matrix and frequency data
main_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy';
result_dir = fullfile(main_dir, 'step1_ec_validation/part2_simulation/supplementary');
atlas_list = {'Schaefer100', 'MMP360'};

%% Induce directionality
% Set up containers
directedSCs = cell(length(atlas_list), 1);
SCs = cell(length(atlas_list), 1);

% Process each atlas
for atlas = 1:length(atlas_list)
    parcellation = atlas_list{atlas};
    SC = importdata(fullfile(main_dir, 'data', [parcellation, '_SC.mat']));
    SCs{atlas} = SC/max(SC(:));
    
    % Apply randmio_dir_connected only once with parameter 1
    directedSCs{atlas} = randmio_dir_connected(SC, 1);
end

%% Topological properties
% Calculate topological properties
for atlas = 1:length(atlas_list)
    parcellation = atlas_list{atlas};
    SC = SCs{atlas};
    directed_SC = directedSCs{atlas};

    % Calculate clustering coefficients
    clustering_sc = clustering_coef_wu(SC);
    clustering_dsc  = clustering_coef_wd(directed_SC);
    corr_clustering_coef = corr(clustering_sc, clustering_dsc);

    % Calculate degree
    SC_out = sum(SC, 1)';
    directed_out = sum(directed_SC, 1)';
    corr_degree = corr(SC_out, directed_out);

    % Create figure with horizontal subplots
    figure('Position', [100 100 700 300])
    
    % Subplot 1: Clustering coefficient
    subplot(1,2,1)
    scatter(clustering_sc, clustering_dsc, 70, 'filled', ...
        'MarkerFaceColor', [0.2 0.6 0.8], ...
        'MarkerFaceAlpha', 0.6, ...
        'MarkerEdgeColor', 'white', ...
        'MarkerEdgeAlpha', 0.5)
    axis square
    box on
    grid on
    set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 4))
    set(gca, 'XTick', linspace(min(get(gca, 'XTick')), max(get(gca, 'XTick')), 4))
    % set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
    
    % Subplot 2: Betweenness centrality
    subplot(1,2,2)
    scatter(SC_out, directed_out, 70, 'filled', ...
        'MarkerFaceColor', [0.8 0.2 0.2], ...
        'MarkerFaceAlpha', 0.6, ...
        'MarkerEdgeColor', 'white', ...
        'MarkerEdgeAlpha', 0.5)
    axis square
    box on
    grid on
    set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 4))
    set(gca, 'XTick', linspace(min(get(gca, 'XTick')), max(get(gca, 'XTick')), 4))
    % set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
    
    % Save figure
    % exportgraphics(gcf, fullfile(result_dir, ['topological_properties_', parcellation, '.png']), 'Resolution', 1200);
end

%% Visualize connectivity matrices
% Create and save visualizations for each atlas
for atlas = 1:length(atlas_list)
    parcellation = atlas_list{atlas};
    SC = SCs{atlas};
    
    % Plot 1: Original SC
    figure('Position', [100 100 300 300])
    imagesc(log(SC))
    colormap(generateColorMap(abs(log(SC(:)+1e-4)),100)) % Using viridis colormap for better visualization
    axis square
    set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
    exportgraphics(gcf, fullfile(result_dir, ['original_sc_', parcellation, '.png']), 'Resolution', 1200);

    % Plot 2: Directed SC
    directed_SC = directedSCs{atlas};
    figure('Position', [100 100 300 300])
    imagesc(log(directed_SC))
    colormap(generateColorMap(abs(log(directed_SC(:)+1e-4)),100)) 
    axis square
    set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
    exportgraphics(gcf, fullfile(result_dir, ['directed_sc_', parcellation, '.png']), 'Resolution', 1200);
end

% Create and save a separate colorbar figure
figure('Position', [100 100 200 300])
c = colorbar('Location', 'west');
c.Ticks = [];  % Remove ticks
c.LineWidth = 1.5;
colormap(generateColorMap(abs(log(directed_SC(:)+1e-4)),100)) 
axis off
exportgraphics(gcf, fullfile(result_dir, 'custom_hot_colorbar.png'), 'Resolution', 1200);

close all

%% Calculate correlations between SC and directed SC degrees
% Create arrays to store correlation values
out_degree_corrs = zeros(length(atlas_list), 1);

% Process each atlas
for atlas = 1:length(atlas_list)
    parcellation = atlas_list{atlas};
    SC = SCs{atlas};
    directed_SC = directedSCs{atlas};
    
    % Calculate degrees
    SC_out = sum(SC, 1)';
    directed_out = sum(directed_SC, 1)';
    
    % Calculate correlation between SC_out and directed_out
    out_degree_corrs(atlas) = corr(SC_out, directed_out);
    
    % Create scatter plot
    figure('Position', [100 100 300 300])
    scatter(SC_out, directed_out, 70, 'filled', ...
        'MarkerFaceColor', [0.8 0.2 0.2], ...
        'MarkerFaceAlpha', 0.6, ...
        'MarkerEdgeColor', 'white', ...
        'MarkerEdgeAlpha', 0.5)
    
    % Add formatting
    hold on
    set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 4))
    set(gca, 'XTick', linspace(min(get(gca, 'XTick')), max(get(gca, 'XTick')), 4))
    set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
    axis square
    box on
    grid on
    hold off

    
    % Save figure
    exportgraphics(gcf, fullfile(result_dir, ['out_degree_comparison_', parcellation, '.png']), 'Resolution', 1200);
end

% Display correlation results
fprintf('Out-degree correlations:\n');
for atlas = 1:length(atlas_list)
    fprintf('%s: %.4f\n', atlas_list{atlas}, out_degree_corrs(atlas));
end
