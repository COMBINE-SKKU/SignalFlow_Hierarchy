%% This script is written to replicate the figure 3 of the paper
clear; close all; clc

% Define paths
data_path = '../../data';
data_file_EC = fullfile(data_path, 'MMP_resting_iEC.mat');
data_file_Yeo = fullfile(data_path, 'Yeo7MMP.mat');

% Import Data
EC = importdata(data_file_EC);
modules = importdata(data_file_Yeo);

%====================================%
%% Figure 3a
%====================================%

% Flatten matrix into a vector
X = EC(:);

% Uniform weights (since no specific weighting is given)
weights = ones(size(EC));

% Sort X in descending order
[X_sorted, sort_idx] = sort(X, 'descend');

% Flatten the weight matrix into a vector
W = weights(:);

% Sort weights to match the sorted X values
W_sorted = W(sort_idx);

% Choose the number of largest observations (k)
% Since the matrix is 360x360, let's choose a reasonable k, say 5% of the data
n = length(X_sorted);
k = round(n * 0.05);  % Select the top 5% of the data for tail analysis

% Call the TailWHill function
[Hill, Hillsd, DJV1, DJV2, AM, AMsd, T1, T1sd, T2, T3, D, Dsd] = TailWHill(X_sorted, W_sorted, k);

% Display the Hill estimate (tail index)
disp(['Hill tail index estimate: ', num2str(Hill)]);

% Extract left cortex of iEC matrix
EC_adjusted = EC(1:180, 1:180);
modules_adjusted = modules(1:180);

% Sort out the left cortex according the the 7 networks
[~, sortedIndices] = sort(modules_adjusted);
EC_sorted = EC_adjusted(sortedIndices, sortedIndices);

% Plot sorted EC matrix
figure;
imagesc(EC_sorted);
colormap(generateColorMap(EC(:), 1000));
axis square;
set(gca, 'XTick', [], 'YTick', []);
% savefig('Figure3a.fig');

%====================================%
%% Figure 3b
%====================================%

% Histogram showing edge strength distribution
figure;
histogram(EC_adjusted(:), 100, 'FaceColor', '#808080', 'Normalization', 'probability');

% Set the figure and axes properties
set(gcf, 'Position', [100, 100, 225, 225]); % Updated position
set(gca, 'XTick', [0 0.4 0.8], 'YTick', [0 0.2 0.4]);
set(gca, 'FontName', 'Helvetica', 'FontSize', 10);
xlabel('Edge Strength');
ylabel('Probability');

% Main axis limits
mainAx = gca;
mainAxLimits = axis(mainAx);

% Create an inset axes on the bottom part of the current axes for histogram
% Position inset relative to the new figure size
insetAxHist = axes('Position', [0.6 0.22 0.28 0.1]); 
% Histogram for the inset, focusing on a specific range
histogram(insetAxHist, EC_adjusted(:), 'BinLimits', [0.4, 0.8], 'FaceColor', '#808080', 'Normalization', 'probability');

% Set the axis limits and ticks for the histogram inset
set(insetAxHist, 'XLim', [0.4, 0.8], 'YLim', [0, 0.002]);
set(insetAxHist, 'XTick', [0.4, 0.8], 'YTick', [0, 0.002]);
set(insetAxHist, 'XTicklabel', []);
set(insetAxHist, 'FontSize', 7);
box off

% Create an inset axes on the upper part of the current axes for bar plot
% Position inset relative to the new figure size
insetAxBar = axes('Position', [0.6 0.6 0.28 0.28]); 

% Calculate proportions for the bar plot
totalConnections = numel(EC_adjusted(:));
numNegativeConnections = sum(EC_adjusted(:) < 0);
numPositiveConnections = sum(EC_adjusted(:) > 0);

% Bar plot for the inset
bar(insetAxBar, [1, 2], [numNegativeConnections, numPositiveConnections] / totalConnections, 'FaceColor', 'flat', ...
    'CData', [0 0 1; 1 0 0]); % Blue for negative, red for positive
set(insetAxBar, 'XTickLabel', {'Neg', 'Pos'}, 'XTick', [1 2], 'ylim', [0 0.75], 'YTickLabel', [0 0.6], 'FontSize', 7);
box on

% Save the figure with the insets
% savefig(gcf, 'Figure3b.fig');
% exportgraphics(gcf, 'Figure3b.png', 'Resolution', 2000);

%====================================%
%% Figure 3c 
%====================================%

% Number of nodes in the connectivity matrix 'EC'
num_nodes = size(EC, 2); 

% Initialize arrays to store the ratios for outgoing and incoming connections
ratio_outgoing_pos = zeros(1, num_nodes); % Ratio of positive outgoing connections
ratio_outgoing_neg = zeros(1, num_nodes); % Ratio of negative outgoing connections
ratio_incoming_pos = zeros(1, num_nodes); % Ratio of positive incoming connections
ratio_incoming_neg = zeros(1, num_nodes); % Ratio of negative incoming connections

% Calculate the ratios for each node
for i = 1:num_nodes
    % Outgoing connections for node i (column-wise)
    outgoing_connections = EC(:, i);
    ratio_outgoing_pos(i) = sum(outgoing_connections > 0) / sum(outgoing_connections ~= 0);
    ratio_outgoing_neg(i) = sum(outgoing_connections < 0) / sum(outgoing_connections ~= 0);
    
    % Incoming connections for node i (row-wise)
    incoming_connections = EC(i, :);
    ratio_incoming_pos(i) = sum(incoming_connections > 0) / sum(incoming_connections ~= 0);
    ratio_incoming_neg(i) = sum(incoming_connections < 0) / sum(incoming_connections ~= 0);
end

% Initialize arrays to hold module-specific data
module_values_pos_out = []; % Values for positive outgoing connections by module
module_values_neg_out = []; % Values for negative outgoing connections by module
module_nums = []; % Module numbers corresponding to the values
module_order = [2 1 4 3 5 6 7]; % Custom order of modules

% Collect data for each module in the specified order
for m = 1:length(module_order)
    current_module = module_order(m);
    module_mask = (modules == current_module);
    
    % Collect and store positive outgoing values
    module_values_pos_out = [module_values_pos_out; ratio_outgoing_pos(module_mask)'];
    
    % Collect and store negative outgoing values
    module_values_neg_out = [module_values_neg_out; ratio_outgoing_neg(module_mask)'];
    
    % Store module number for each entry
    module_nums = [module_nums; m * ones(sum(module_mask), 1)];
end

% Combine positive and negative outgoing values
combined_module_vals = [module_values_pos_out, module_values_neg_out];

% Reshape combined values for easier table construction
reshaped_vals = reshape(combined_module_vals, [], 1);

% Generate origin labels ('posOut', 'negOut') for each value
num_values = numel(ratio_outgoing_pos); % Assuming equal numbers of positive and negative values
origin_labels = repmat({'posOut', 'negOut'}, num_values, 1);
reshaped_labels = reshape(origin_labels, [], 1);

% Duplicate module numbers to match reshaped data size
color_labels = [module_nums; module_nums];

% Combine data into a table
result = [array2table([reshaped_vals, color_labels]), cell2table(reshaped_labels)];

% Save the resulting table to a .mat file
% save('Figure3c.mat', 'result');


%====================================%
%% Figure 3d 
%====================================%

% Number of modules
nModules = 22;

% Initialize matrices for positive and negative signal flow
signal_flow_pos = zeros(nModules, nModules);
signal_flow_neg = zeros(nModules, nModules);

% Compute signal flow matrices in parallel
parfor i = 1:nModules
    [signal_flow_neg(:, i), signal_flow_pos(:, i)] = signalflow_modules(EC, i, '22modules');
end

% Plot the signal flow using edge bundling
plot_signal_flow(nModules, signal_flow_pos, signal_flow_neg)

% Save the figure
% savefig(gcf, 'energyflow_22module_rest.fig');  % Save the figure in .fig format
% exportgraphics(gcf, 'energyflow_22module_rest.png', 'Resolution', 1200);  % Export the figure as a high-resolution PNG

%====================================%
%% Figure 3e
%====================================%
%% Cross module analysis
% Define the module information with indices for Unimodal and Heteromodal categories
ModuleInfo = {'Unimodal', 1:11, 'Heteromodal', 12:22};

% Extract indices for unimodal and heteromodal modules
unimodal_idx = ModuleInfo{2};
heteromodal_idx = ModuleInfo{4};

% Reshape the positive signal flow data and filter out zero values
within_unimodal_pos = nonzeros(reshape(signal_flow_pos(unimodal_idx, unimodal_idx), [], 1));
unimodal_to_heteromodal_pos = nonzeros(reshape(signal_flow_pos(heteromodal_idx, unimodal_idx), [], 1));
heteromodal_to_unimodal_pos = nonzeros(reshape(signal_flow_pos(unimodal_idx, heteromodal_idx), [], 1));
within_heteromodal_pos = nonzeros(reshape(signal_flow_pos(heteromodal_idx, heteromodal_idx), [], 1));

% Reshape the negative signal flow data, take the absolute value, and filter out zeros
within_unimodal_neg = nonzeros(abs(reshape(signal_flow_neg(unimodal_idx, unimodal_idx), [], 1)));
unimodal_to_heteromodal_neg = nonzeros(abs(reshape(signal_flow_neg(heteromodal_idx, unimodal_idx), [], 1)));
heteromodal_to_unimodal_neg = nonzeros(abs(reshape(signal_flow_neg(unimodal_idx, heteromodal_idx), [], 1)));
within_heteromodal_neg = nonzeros(abs(reshape(signal_flow_neg(heteromodal_idx, heteromodal_idx), [], 1)));

% Combine positive and negative flow data for bar plot preparation
pos_data = [within_unimodal_pos; unimodal_to_heteromodal_pos; heteromodal_to_unimodal_pos; within_heteromodal_pos];
neg_data = [within_unimodal_neg; unimodal_to_heteromodal_neg; heteromodal_to_unimodal_neg; within_heteromodal_neg];

% Combine all data into one array
Data = [pos_data; neg_data];

% Create group labels for the data (positive and negative)
xGroup1 = repmat({'positive'}, length(pos_data), 1);  % Group for positive data
xGroup2 = repmat({'negative'}, length(neg_data), 1);  % Group for negative data
xGroup = [xGroup1; xGroup2];  % Combine group labels

% Create labels for the different categories within the groups
temp1 = repmat({'within_unimodal'}, length(within_unimodal_pos), 1);
temp2 = repmat({'unimodal_to_heteromodal'}, length(unimodal_to_heteromodal_pos), 1);
temp3 = repmat({'heteromodal_to_unimodal'}, length(heteromodal_to_unimodal_pos), 1);
temp4 = repmat({'within_heteromodal'}, length(within_heteromodal_pos), 1);
temp5 = repmat({'within_unimodal'}, length(within_unimodal_neg), 1);
temp6 = repmat({'unimodal_to_heteromodal'}, length(unimodal_to_heteromodal_neg), 1);
temp7 = repmat({'heteromodal_to_unimodal'}, length(heteromodal_to_unimodal_neg), 1);
temp8 = repmat({'within_heteromodal'}, length(within_heteromodal_neg), 1);
Legend = [temp1; temp2; temp3; temp4; temp5; temp6; temp7; temp8];  % Combine all labels

% Combine the data and labels into a table
result = [array2table(Data), cell2table([xGroup Legend])];

% Save the result table to a .mat file
% save('Figure3e.mat', 'result');


