%% Simple Visualization of Figure 3c
% This script generates boxplots for ratio values from Figure 3c.
% One should use BoxPlotPro from MATLAB FileExchange to reproduce figures
% from the manuscript with the same quality.

clear; close all; clc;

% Load ratio data
ratio_array = importdata(fullfile('Figure3c.mat'));

% Predefine the order for Var2 and reshaped_labels
var2_order = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7];
label_order = {'posOut', 'negOut', 'posOut', 'negOut', ...
               'posOut', 'negOut', 'posOut', 'negOut', ...
               'posOut', 'negOut', 'posOut', 'negOut', ...
               'posOut', 'negOut'};

% Mapping of Var2 values to names
var2_names = {'SM', 'Vis', 'DAN', 'VAN', 'Limbic', 'FPN', 'DMN'};

% Initialize empty arrays to store the results
data = [];
labels = {};

% Loop through the desired sequence of Var2 and reshaped_labels
for i = 1:length(var2_order)
    % Extract rows matching the current Var2 and reshaped_labels condition
    matching_rows = ratio_array.Var2 == var2_order(i) & strcmp(ratio_array.reshaped_labels, label_order{i});
    
    % Append the corresponding Var1 data to 'data' array
    data = [data; ratio_array.Var1(matching_rows)];
    
    % Create detailed labels combining Var2 names and reshaped_labels
    % Get the name corresponding to the current Var2 value
    current_var2_name = var2_names{var2_order(i)};
    
    % Create the label for the current Var2 and reshaped_labels
    current_label = strcat(current_var2_name, '_', label_order{i});
    
    % Get the number of matching rows
    num_matching = sum(matching_rows);  
    
    % Repeat the label for all matching rows and append to 'labels' array
    labels = [labels; repmat({current_label}, num_matching, 1)];
end

% Create the boxplot and set box colors
figure;
h = boxplot(data, labels, 'Colors', 'k', 'Symbol', ''); % 'k' sets the box color to black

% Define colors for each box
col_vals = [...
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00'];
col_vals = reshape(hex2dec(col_vals), [3 numel(col_vals)/6])' ./ 255;
P = size(col_vals, 1);
colors = interp1(1:size(col_vals, 1), col_vals, linspace(1, P, 14), 'linear');

% Set the colors for each box
hGroups = findobj(gca, 'Tag', 'Box');
for j = 1:length(hGroups)
    patch(get(hGroups(j), 'XData'), get(hGroups(j), 'YData'), colors(j, :), 'FaceAlpha', 0.8); % Fill the box with color
end

% Hold the current plot
hold on;


% Add scatter plot points in grey color
unique_groups = {...
    'SM_posOut', 'SM_negOut', ...
    'Vis_posOut', 'Vis_negOut', ...
    'DAN_posOut', 'DAN_negOut', ...
    'VAN_posOut', 'VAN_negOut', ...
    'Limbic_posOut', 'Limbic_negOut', ...
    'FPN_posOut', 'FPN_negOut', ...
    'DMN_posOut', 'DMN_negOut'};
num_groups = length(unique_groups);

% Adjust jitter for better visibility
jitterAmount = 0.1;

for i = 1:num_groups
    % Find the x-position for each group
    x_pos = repmat(i, sum(strcmp(labels, unique_groups{i})), 1) + (rand(sum(strcmp(labels, unique_groups{i})), 1) - 0.5) * jitterAmount;
    y_pos = data(strcmp(labels, unique_groups{i}));
    scatter(x_pos, y_pos, 5, 'filled', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', [0.5, 0.5, 0.5]); % Grey color
end

% Tilt the x-ticks by 45 degrees
xtickangle(45);

%% Simple Visualization of Figure 3e
% This script generates boxplots for ratio values from Figure 3e.
% One should use BoxPlotPro from MATLAB FileExchange to reproduce figures
% from the manuscript with the same quality.

clear; close all; clc;

% Load ratio data
ratio_array = importdata(fullfile('Figure3e.mat'));

% Predefine the order for Var2 and Var1 (reshaped_labels equivalent)
var2_order = {'within_unimodal', 'within_unimodal', ...
              'unimodal_to_heteromodal', 'unimodal_to_heteromodal', ...
              'heteromodal_to_unimodal', 'heteromodal_to_unimodal', ...
              'within_heteromodal', 'within_heteromodal'};
var1_order = {'positive', 'negative', ...
              'positive', 'negative', ...
              'positive', 'negative', ...
              'positive', 'negative'};

% Initialize empty arrays to store the results
data = [];
labels = {};

% Loop through the desired sequence of Var2 and reshaped_labels
for i = 1:length(var2_order)
    % Extract rows matching the current Var2 and reshaped_labels condition
    matching_rows = strcmp(ratio_array.Var2, var2_order{i}) & strcmp(ratio_array.Var1, var1_order{i});
    
    % Append the corresponding 'Data' values to 'data' array
    data = [data; ratio_array.Data(matching_rows)];
    
    % Create labels combining Var2 names and Var1 (reshaped_labels equivalent)
    current_label = strcat(var2_order{i}, '_', var1_order{i});
    
    % Get the number of matching rows
    num_matching = sum(matching_rows);
    
    % Append the labels for the matching rows
    labels = [labels; repmat({current_label}, num_matching, 1)];
end

% Create the boxplot and set box colors
figure;
h = boxplot(data, labels, 'Colors', 'k', 'Symbol', ''); % 'k' sets the box color to black

% Define colors for each box (hex values)
col_vals = [...
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00';];  % Color 8

% Convert hex to RGB values and normalize to [0, 1] range
col_vals = reshape(hex2dec(col_vals), [3 numel(col_vals)/6])' ./ 255;
P = size(col_vals, 1);
colors = interp1(1:size(col_vals, 1), col_vals, linspace(1, P, 8), 'linear');

% Set the colors for each box
hGroups = findobj(gca, 'Tag', 'Box');
for j = 1:length(hGroups)
    patch(get(hGroups(j), 'XData'), get(hGroups(j), 'YData'), col_vals(mod(j-1, 8)+1, :), 'FaceAlpha', 0.8); % Fill the box with color
end

% Hold the current plot
hold on;

% Define unique groups in the desired order for scatter plotting
unique_groups_ordered = {...
    'within_unimodal_positive', 'within_unimodal_negative', ...
    'unimodal_to_heteromodal_positive', 'unimodal_to_heteromodal_negative', ...
    'heteromodal_to_unimodal_positive', 'heteromodal_to_unimodal_negative', ...
    'within_heteromodal_positive', 'within_heteromodal_negative'};

% Adjust jitter for better visibility
jitterAmount = 0.1;

% Loop through each group and plot scatter points
for i = 1:length(unique_groups_ordered)
    % Find the x-position for each group based on the ordered labels
    x_pos = repmat(i, sum(strcmp(labels, unique_groups_ordered{i})), 1) + (rand(sum(strcmp(labels, unique_groups_ordered{i})), 1) - 0.5) * jitterAmount;
    y_pos = data(strcmp(labels, unique_groups_ordered{i}));
    
    % Scatter plot for individual data points
    scatter(x_pos, y_pos, 5, 'filled', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', [0.5, 0.5, 0.5]); % Grey color
end

% Tilt the x-ticks by 45 degrees for readability
xtickangle(45);
hold off;