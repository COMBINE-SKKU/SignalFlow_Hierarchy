%% Save edge strength correlation and degree correlation data for box plot
% This script prepares and saves correlation data for box plot visualization. 
% The data is reshaped, concatenated, and saved as a table. Use BoxPlotPro from 
% MATLAB File Exchange to visualize the data with high-quality box plots.
close all; clear; clc

% Define parameters
atlas = 1; % Atlas index: 1 for Schaefer-100, 2 for MMP-360
num_iter = 100; % Number of iterations used in simulations

% Load precomputed correlation data
EC_correlation_values = importdata(fullfile('fig2_results/EC_correlation_values.mat'));
in_degree_corr_values = importdata(fullfile('fig2_results/EC_in_degree_corr_values.mat'));
out_degree_corr_values = importdata(fullfile('fig2_results/EC_out_degree_corr_values.mat'));

%% Prepare data for edge strength correlation box plot

% Extract correlation values for the selected atlas
correlation_values1 = EC_correlation_values(:, :, atlas); % Data for Schaefer-100
correlation_values2 = EC_correlation_values(:, :, atlas + 1); % Data for MMP-360

% Reshape the data vectors for each measure
data_correlation1 = reshape(correlation_values1, [], 1);
data_correlation2 = reshape(correlation_values2, [], 1);

% Concatenate the reshaped data into one vector
data_all = [data_correlation1; data_correlation2];

% Create group labels for each dataset
group_correlation1 = repmat({'Schaefer'}, size(data_correlation1));
group_correlation2 = repmat({'MMP'}, size(data_correlation2));

% Combine the group labels into one vector
group_all = [group_correlation1; group_correlation2];

% Generate legends for each measure
legend = repmat({'rDCM', 'VAR', 'FASK', 'iEC'}, num_iter, 1);
legend = reshape(legend, [], 1);

% Repeat the legend for each dataset and concatenate with group labels
legend_all = [legend; legend];

% Create a table with the concatenated data, legends, and group labels
result = table(data_all, legend_all, group_all, 'VariableNames', {'Data', 'Legend', 'Group'});

% Save the edge strength correlation data table
save('Figure2b_edge_strength_correlation_table', 'result');

%% Prepare data for degree correlation box plot

% Extract the correlation values for the in-degree and out-degree
correlation_values1 = in_degree_corr_values(:,:,atlas);
correlation_values2 = out_degree_corr_values(:,:,atlas);

% Reshape the data vectors for each measure
data_correlation1 = reshape(correlation_values1, [], 1);
data_correlation2 = reshape(correlation_values2, [], 1);

% Concatenate the reshaped data into one vector
data_all = [data_correlation1; data_correlation2];

% Create group labels for each dataset
group_correlation1 = repmat({'Indegree'}, size(data_correlation1));
group_correlation2 = repmat({'Outdegree'}, size(data_correlation2));

% Combine the group labels into one vector
group_all = [group_correlation1; group_correlation2];

% Generate legends for each measure
legend_all1 = repmat({'rDCM', 'VAR', 'FASK', 'iEC'}, num_iter, 1);
legend_all1 = reshape(legend_all1, [], 1);

% Repeat legend for the second dataset
legend_all2 = repmat({'rDCM', 'VAR', 'FASK', 'iEC'}, num_iter, 1);
legend_all2 = reshape(legend_all2, [], 1);

% Concatenate the legends with the group labels
legend_all = [legend_all1; legend_all2];

% Create a table with the concatenated data, legends, and group labels
result = table(data_all, legend_all, group_all, 'VariableNames', {'Data', 'Legend', 'Group'});

% Save the degree correlation data table
save('Figure2c_degree_correlation_Schaefer_table', 'result');

%% Simple plot for the edge strength correlation
% One should use BoxPlotPro from MATLAB FileExchange to reproduce
% figures from the manuscript with the same quality

% Extract correlation values for the selected atlas
correlation_values1 = EC_correlation_values(:, :, atlas); % Data for Schaefer-100
correlation_values2 = EC_correlation_values(:, :, atlas + 1); % Data for MMP-360

% Reshape the data vectors for each measure
data_correlation1 = reshape(correlation_values1, [], 1);
data_correlation2 = reshape(correlation_values2, [], 1);

% Concatenate the reshaped data into one vector
data_all = [data_correlation1; data_correlation2];

% Create group labels for each dataset
group_correlation1 = repmat({'Schaefer'}, size(data_correlation1));
group_correlation2 = repmat({'MMP'}, size(data_correlation2));

% Combine the group labels into one vector
group_all = [group_correlation1; group_correlation2];

% Generate legends for each measure
num_iter = size(correlation_values1, 1);
legend = repmat({'rDCM', 'VAR', 'FASK', 'iEC'}, num_iter, 1);
legend = reshape(legend, [], 1);

% Repeat the legend for each dataset and concatenate with group labels
legend_all = [legend; legend];

% Combine group labels and legends
combined_group = strcat(group_all, '_', legend_all);

% Create the boxplot and set box colors
figure;
h = boxplot(data_all, combined_group, 'Colors', 'k', 'Symbol', ''); % 'k' sets the box color to black

% Set the colors for each box
col_vals = [...
'fd';'e3';'a2';
'f9';'dc';'e1';
'9b';'c5';'db';
'f6';'69';'45';
'fd';'e3';'a2';
'f9';'dc';'e1';
'9b';'c5';'db';
'f6';'69';'45';];
col_vals = reshape(hex2dec(col_vals), [3 numel(col_vals)/6])' ./ 255;
P = size(col_vals,1);
colors = interp1(1:size(col_vals,1), col_vals, linspace(1,P,8), 'linear');      

% Set the colors for each box
hGroups = findobj(gca, 'Tag', 'Box');
for j = 1:length(hGroups)
    patch(get(hGroups(j), 'XData'), get(hGroups(j), 'YData'), colors(j, :), 'FaceAlpha', .5); % Fill the box with color
end

% Hold the current plot
hold on;

% Add scatter plot points in grey color
unique_groups = {'Schaefer_rDCM';'Schaefer_VAR';'Schaefer_FASK';'Schaefer_iEC';...
    'MMP_rDCM';'MMP_VAR';'MMP_FASK';'MMP_iEC'};
num_groups = length(unique_groups);

% Adjust jitter for better visibility
jitterAmount = 0.1;

for i = 1:num_groups
    % Find the x-position for each group
    x_pos = repmat(i, sum(strcmp(combined_group, unique_groups{i})), 1) + (rand(sum(strcmp(combined_group, unique_groups{i})), 1) - 0.5) * jitterAmount;
    y_pos = data_all(strcmp(combined_group, unique_groups{i}));
    scatter(x_pos, y_pos, 2, 'filled', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5]); % Grey color
end

% Rotate x-ticks and set font size
ax = gca;
ax.XTickLabelRotation = 45;
ax.FontSize = 10; % Adjust the font size as needed
