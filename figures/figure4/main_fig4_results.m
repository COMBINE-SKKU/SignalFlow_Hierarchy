%% This script is written to replicate the figure 4 of the paper
clear; close all; clc

% Define paths
data_path = '../../data';
data_file_EC = fullfile(data_path, 'MMP_resting_iEC.mat');

% Import Data
EC = importdata(data_file_EC);

%====================================%
%% Figure 4b
%====================================%

% compute hierarchical levels
hierarchyLevels = computeHierarchyLevels(EC, 0.15);

% plot the result
surfaceplot(hierarchyLevels,'MMP','both','viridis')
% exportgraphics(gcf,'Figure4b.png','Resolution',2000)

%====================================%
%% Figure 4c
%====================================%

% Import Data from the specified file and create a surface plot
data_file_PC = fullfile(data_path, 'MMP_PC1.mat');
pc1 = importdata(data_file_PC);
surfaceplot(pc1, 'MMP', 'both', 'pc1');  
% exportgraphics(gcf, 'Figure5c_pc1.png', 'Resolution', 2000);  

% Normalize the data for the scatter plot
range = max(pc1(:)) - min(pc1(:));  
pc1 = ((pc1 - min(pc1(:))) / range) * 10;  

% Prepare data for scatter plot
x1 = pc1;  
y1 = hierarchyLevels;  

% Calculate the correlation between the x and y data
xyCorr = corr(x1, y1);  

% Create the scatter plot
figure;
scatter(x1, y1, 'o', 'SizeData', 50, 'MarkerFaceColor', [0.15, 0.15, 0.15], ...
    'MarkerEdgeColor', [1, 1, 1], 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.8); 
hold on;

% Calculate mean and standard deviation for x and y data
xMean = mean(x1);  
xStd = std(x1); 
xLine = round(xMean + xStd, 1);  

yMean = mean(y1); 
yStd = std(y1);  
yLine = round(yMean + yStd, 1);  

% Display vertical and horizontal lines on the scatter plot
xline(xLine, '--');  
yline(yLine, '--');  

% Identify points of interest based on their position relative to the lines
HigherThanPC = find(y1 > yLine & x1 < xLine);  
scatter(x1(HigherThanPC), y1(HigherThanPC), 'r');  

LowerThanPC = find(y1 < yLine & x1 > xLine);  
scatter(x1(LowerThanPC), y1(LowerThanPC), 'b');  

% Optionally, fit a linear model and plot the fitted line
mdl = fitlm(x1, y1);  
plot(x1, mdl.Fitted, 'k-');  

% Save the figure
% savefig('Figure4c.fig');

% Create a surface plot highlighting points above yLine and to the left of xLine
z = zeros(360, 1);  
z(HigherThanPC) = 1;  
surfaceplot(z, 'MMP', 'both', 'viridis2');  
% exportgraphics(gcf, 'Figure4c_pc1_hierarchy_difference1.png', 'Resolution', 2000);  

% Create a surface plot highlighting points below yLine and to the right of xLine
z = zeros(360, 1); 
z(LowerThanPC) = 1;
surfaceplot(z, 'MMP', 'both', 'viridis2'); 
% exportgraphics(gcf, 'Figure4c_pc1_hierarchy_difference2.png', 'Resolution', 2000);  

%====================================%
%% Figure 4e
%====================================%
% Define the file path for the cortical types data
data_file_corTypes = fullfile(data_path, 'cortical_types.mat');

% Import cortical types data from the specified file
corticalTypes = importdata(data_file_corTypes);

% Create a surface plot of the cortical types data
surfaceplot(corticalTypes, 'MMP', 'both', 'tpl');

% Export the surface plot as a high-resolution image
% exportgraphics(gcf, 'Figure4e.png', 'Resolution', 2000);

% Set the number of cortical regions and their corresponding labels
num_regions = 6;
cyto_label = {'Konicortex', 'Eulaminate-III', 'Eulaminate-II', 'Eulaminate-I', 'Dysgranular', 'Agranular'};

% Generate and plot boxplots for different cytoarchitectonic groups based on hierarchy levels
plotCytoGroups(hierarchyLevels, corticalTypes, cyto_label, 'tpl', num_regions);

% Save the boxplot figure
% savefig('Figure4eBoxplot.fig');

cyto_data = [];
group_count = 0;
group_labels = [];
% Compute median values and prepare boxplot data for each region
for i = 1:num_regions
    val_cyto = hierarchyLevels(corticalTypes == i);
    group_count = group_count+1;
    cyto_data = [cyto_data ;val_cyto];
    group_labels = [group_labels; repmat(group_count, length(val_cyto), 1)];
end

glm_model = fitglm(group_labels, cyto_data, 'linear');

disp(glm_model);

%====================================%
%% Figure 4f
%====================================%
% Number of modules
data_file_modules = fullfile(data_path, 'module27(ordered).mat');

% Import 
MMPmodule27 = importdata(data_file_modules);
surfaceplot(MMPmodule27, 'MMP', 'both', 'module27');
exportgraphics(gcf,'Figure4f(brainmap).png','Resolution',2000)

nModules = 27;

% Initialize matrices for positive and negative signal flow
signal_flow_pos = zeros(nModules, nModules);
signal_flow_neg = zeros(nModules, nModules);

% Compute signal flow matrices in parallel
parfor i = 1:nModules
    [signal_flow_neg(:, i), signal_flow_pos(:, i)] = signalflow_modules(EC, i, '27modules');
end

% Plot the signal flow using edge bundling
plot_signal_flow(nModules, signal_flow_pos, signal_flow_neg)

% Save the figure
% savefig(gcf, 'energyflow_27module_rest.fig');  % Save the figure in .fig format
% exportgraphics(gcf, 'energyflow_27module_rest.png', 'Resolution', 1200);  % Export the figure as a high-resolution PNG

%% joy plot 
nModules = 360;

% Initialize matrices for positive and negative signal flow
signal_flow_pos = zeros(nModules, nModules);
signal_flow_neg = zeros(nModules, nModules);

% Compute signal flow matrices in parallel
parfor i = 1:nModules
    [signal_flow_neg(:, i), signal_flow_pos(:, i)] = signalflow_modules(EC, i, 'nomodule');
end

nModules = 27;
outgoing_neg = abs(sum(signal_flow_neg,1));
outgoing_pos = sum(signal_flow_pos,1);

avg_out_neg = zeros(1,nModules);
avg_out_pos = zeros(1,nModules);

for i = 1:nModules
    temp_pos = outgoing_pos(MMPmodule27==i);
    temp_neg = outgoing_neg(MMPmodule27==i);
    
    avg_out_pos(i) = mean(temp_pos);
    avg_out_neg(i) = mean(temp_neg);
end

%%
figure;

vector1 = avg_out_neg; 
vector2 = avg_out_pos; 

% Normalize the data
min_neg = min(vector1);
max_neg = max(vector1);
normalized_neg = (vector1 - min_neg) / (max_neg - min_neg);

min_pos = min(vector2);
max_pos = max(vector2);
normalized_pos = (vector2 - min_pos) / (max_pos - min_pos);

% Interpolate the normalized vectors to have the same number of points
num_points = max(length(normalized_neg), length(normalized_pos));
x_interp = linspace(1, num_points, num_points);

% Add spacing between the curves (joy plot effect)
spacing = 0.7; % Adjust this value to change the space between lines
normalized_pos = normalized_pos + spacing;

% Plot the lines of the vectors
hold on;
plot(x_interp, normalized_neg, '-o', 'Color', [0.0824, 0.3765, 0.5098], 'LineWidth', 1.5);
plot(x_interp, normalized_pos, '-o', 'Color', [0.7373, 0.0078, 0.0078], 'LineWidth', 1.5);

% Shade the area under the lines
fill([x_interp fliplr(x_interp)], [zeros(1, num_points) fliplr(normalized_neg)], [0.0824, 0.3765, 0.5098], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([x_interp fliplr(x_interp)], [spacing * ones(1, num_points) fliplr(normalized_pos)], [0.7373, 0.0078, 0.0078], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Customize the plot to give an impression like a joy plot
set(gca, 'XColor', 'none'); % Hide X axis
set(gca, 'YColor', 'none'); % Hide Y axis
box off;
hold off;
set(gca, 'Visible', 'off');  % Hide the axes 

% savefig(gcf, 'Figure4f_joy.fig');  

