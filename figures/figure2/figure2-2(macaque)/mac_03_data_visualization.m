%% Simple Visualization of Figure 2e
% This script generates boxplots for correlation values from Figure 2e.
% One should use BoxPlotPro from MATLAB FileExchange to reproduce figures
% from the manuscript with the same quality.

close all; clear; clc;

% Load correlation data
correlation_values = importdata(fullfile('fig2(mac)results', 'EC_correlation_FLN.mat'));

% Reshape the data vectors for each measure
data = reshape(correlation_values, [], 1);

% Generate legends for each measure
num_iter = size(correlation_values, 1);
legend = repmat({'rDCM', 'VAR', 'FASK', 'iEC'}, num_iter, 1);
legend = reshape(legend, [], 1);

% Create the boxplot and set box colors
figure;
h = boxplot(data, legend, 'Colors', 'k', 'Symbol', ''); % 'k' sets the box color to black

% Define colors for each box
col_vals = [...
    'fd'; 'e3'; 'a2';
    'f9'; 'dc'; 'e1';
    '9b'; 'c5'; 'db';
    'f6'; '69'; '45'];
col_vals = reshape(hex2dec(col_vals), [3 numel(col_vals)/6])' ./ 255;
P = size(col_vals, 1);
colors = interp1(1:size(col_vals, 1), col_vals, linspace(1, P, 4), 'linear');

% Set the colors for each box
hGroups = findobj(gca, 'Tag', 'Box');
for j = 1:length(hGroups)
    patch(get(hGroups(j), 'XData'), get(hGroups(j), 'YData'), colors(j, :), 'FaceAlpha', 0.5); % Fill the box with color
end

% Hold the current plot
hold on;

% Add scatter plot points in grey color
unique_groups = {'rDCM', 'VAR', 'FASK', 'iEC'};
num_groups = length(unique_groups);

% Adjust jitter for better visibility
jitterAmount = 0.1;

for i = 1:num_groups
    % Find the x-position for each group
    x_pos = repmat(i, sum(strcmp(legend, unique_groups{i})), 1) + (rand(sum(strcmp(legend, unique_groups{i})), 1) - 0.5) * jitterAmount;
    y_pos = data(strcmp(legend, unique_groups{i}));
    scatter(x_pos, y_pos, 5, 'filled', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', [0.5, 0.5, 0.5]); % Grey color
end

%% Simple Visualization of Figure 2f
% This script generates boxplots for ratio values from Figure 2f.
% One should use BoxPlotPro from MATLAB FileExchange to reproduce figures
% from the manuscript with the same quality.

clear; close all; clc;

% Load ratio data
ratio_array = importdata(fullfile('fig2(mac)results', 'EC_FF_FB.mat'));

% Reshape the data vectors for each measure
data = reshape(ratio_array, [], 1);

% Generate group and legend labels
group = repmat({'FF', 'FF', 'FB', 'FB'}, 10, 1);
group = reshape(group, [], 1);
legend = repmat({'neg', 'pos', 'neg', 'pos'}, 10, 1);
legend = reshape(legend, [], 1);
combined_group = strcat(group, '_', legend);

% Create the boxplot and set box colors
figure;
h = boxplot(data, combined_group, 'Colors', 'k', 'Symbol', ''); % 'k' sets the box color to black

% Define colors for each box
col_vals = [...
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00'];
col_vals = reshape(hex2dec(col_vals), [3 numel(col_vals)/6])' ./ 255;
P = size(col_vals, 1);
colors = interp1(1:size(col_vals, 1), col_vals, linspace(1, P, 4), 'linear');

% Set the colors for each box
hGroups = findobj(gca, 'Tag', 'Box');
for j = 1:length(hGroups)
    patch(get(hGroups(j), 'XData'), get(hGroups(j), 'YData'), colors(j, :), 'FaceAlpha', 0.8); % Fill the box with color
end

% Hold the current plot
hold on;

% Add scatter plot points in grey color
unique_groups = {'FF_neg', 'FF_pos', 'FB_neg', 'FB_pos'};
num_groups = length(unique_groups);

% Adjust jitter for better visibility
jitterAmount = 0.1;

for i = 1:num_groups
    % Find the x-position for each group
    x_pos = repmat(i, sum(strcmp(combined_group, unique_groups{i})), 1) + (rand(sum(strcmp(combined_group, unique_groups{i})), 1) - 0.5) * jitterAmount;
    y_pos = data(strcmp(combined_group, unique_groups{i}));
    scatter(x_pos, y_pos, 5, 'filled', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', [0.5, 0.5, 0.5]); % Grey color
end

