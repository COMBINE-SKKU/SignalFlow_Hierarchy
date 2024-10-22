%% This script is written to replicate the figure 5 of the paper
clear; close all; clc

% Define paths
data_path = '/combinelab/03_user/younghyun/01_project/01_HierarchyMapping/data';
EC_rest = importdata(fullfile(data_path, 'MMP_resting_iEC.mat'));

% Import Data
EC_movie = importdata(fullfile(data_path,'MMP_movie_iEC.mat'));

% Import Data
EC_pain = importdata(fullfile(data_path,'MMP_pain_iEC.mat'));

%====================================%
%% Figure 5a
%====================================%

% compute hierarchical levels
hierarchyRest = computeHierarchyLevels(EC_rest, 0.15);
hierarchyMovie = computeHierarchyLevels(EC_movie, 0.15);
hierarchyPain = computeHierarchyLevels(EC_pain, 0.15);

% plot the result
surfaceplot(hierarchyRest,'MMP','both','viridis')
exportgraphics(gcf,'Figure5aMovie.png','Resolution',2000)

% plot the result
surfaceplot(hierarchyPain,'MMP','both','viridis')
exportgraphics(gcf,'Figure5aPain.png','Resolution',2000)

%====================================%
%% Figure 5b
%====================================%

% Import Cortical Zones Atlas
CortZones = importdata(fullfile(data_path,'/cortical_zones.mat'));
surfaceplot(CortZones,'MMP','both','zone4')
exportgraphics(gcf,'Figure5b(brainmap).png','Resolution',2000)

% Assuming rest_hierarchyLevels, pain_hierarchyLevels, and movie_hierarchyLevels are defined
rest_ratio = zeros(4,1);
pain_ratio = zeros(4,1);
movie_ratio = zeros(4,1);

rest_sem = zeros(4,1);
pain_sem = zeros(4,1);
movie_sem = zeros(4,1);

% Calculate median hierarchy levels and SEM for each category
for cyto_type = 1:4
    indx = cyto_downsample == cyto_type;
    rest_data = hierarchyRest(indx);
    pain_data = hierarchyPain(indx);
    movie_data = hierarchyMovie(indx);
    
    rest_ratio(cyto_type) = median(rest_data);
    pain_ratio(cyto_type) = median(pain_data);
    movie_ratio(cyto_type) = median(movie_data);
    
    % Calculate SEM
    rest_sem(cyto_type) = std(rest_data) / sqrt(numel(rest_data));
    pain_sem(cyto_type) = std(pain_data) / sqrt(numel(pain_data));
    movie_sem(cyto_type) = std(movie_data) / sqrt(numel(movie_data));
end

% Plot the results with error bars
figure;
errorbar(1:4, rest_ratio, rest_sem, '-o'); hold on;
errorbar(1:4, movie_ratio, movie_sem, '-o'); hold on;
errorbar(1:4, pain_ratio, pain_sem, '-o');
legend('rest', 'movie', 'pain')
xticks([1, 2, 3, 4]);
xticklabels({'primary', 'unimodal', 'heteromodal', 'limbic'});
title('Median Hierarchy Levels');
xlabel('Brain Region Categories');
ylabel('Hierarchy Level');
savefig('Figure5b.fig')

%====================================%
%% Figure 5c
%====================================%
% ----------------------------- Rest Signal Flow (for reference) ----------------------------- %
% Number of modules
nModules = 27;

% Initialize matrices for positive and negative signal flow
signal_flow_pos_rest = zeros(nModules, nModules);
signal_flow_neg_rest = zeros(nModules, nModules);

% Compute signal flow matrices in parallel
parfor i = 1:nModules
    [signal_flow_neg_rest(:, i), signal_flow_pos_rest(:, i)] = signalflow_modules(EC_rest, i, '27modules');
end

% ----------------------------- Pain Signal Flow ----------------------------- %

% Initialize matrices for positive and negative signal flow
signal_flow_pos_pain = zeros(nModules, nModules);
signal_flow_neg_pain = zeros(nModules, nModules);

% Compute signal flow matrices in parallel
parfor i = 1:nModules
    [signal_flow_neg_pain(:, i), signal_flow_pos_pain(:, i)] = signalflow_modules(EC_pain, i, '27modules');
end

% Plot the signal flow using edge bundling
plot_signal_flow(nModules, signal_flow_pos_pain, signal_flow_neg_pain)

% Save the figure
savefig(gcf, 'energyflow_27module_pain.fig');  % Save the figure in .fig format

% ----------------------------- Movie Signal Flow ----------------------------- %

% Initialize matrices for positive and negative signal flow
signal_flow_pos_movie = zeros(nModules, nModules);
signal_flow_neg_movie = zeros(nModules, nModules);

% Compute signal flow matrices in parallel
parfor i = 1:nModules
    [signal_flow_neg_movie(:, i), signal_flow_pos_movie(:, i)] = signalflow_modules(EC_movie, i, '27modules');
end

% Plot the signal flow using edge bundling
plot_signal_flow(nModules, signal_flow_pos_movie, signal_flow_neg_movie)

% Save the figure
savefig(gcf, 'energyflow_27module_movie.fig');  % Save the figure in .fig format

%====================================%
%% Delta lines for the signal flow
%====================================%

% Number of modules
nModules = 360;

% Initialize matrices for positive and negative signal flow
signal_flow_pos_rest = zeros(nModules, nModules);
signal_flow_neg_rest = zeros(nModules, nModules);

% Compute signal flow matrices in parallel
parfor i = 1:nModules
    [signal_flow_neg_rest(:, i), signal_flow_pos_rest(:, i)] = signalflow_modules(EC_rest, i, 'nomodule');
end

% ----------------------------- Pain Signal Flow ----------------------------- %

% Initialize matrices for positive and negative signal flow
signal_flow_pos_pain = zeros(nModules, nModules);
signal_flow_neg_pain = zeros(nModules, nModules);

% Compute signal flow matrices in parallel
parfor i = 1:nModules
    [signal_flow_neg_pain(:, i), signal_flow_pos_pain(:, i)] = signalflow_modules(EC_pain, i, 'nomodule');
end


% ----------------------------- Movie Signal Flow ----------------------------- %

% Initialize matrices for positive and negative signal flow
signal_flow_pos_movie = zeros(nModules, nModules);
signal_flow_neg_movie = zeros(nModules, nModules);

% Compute signal flow matrices in parallel
parfor i = 1:nModules
    [signal_flow_neg_movie(:, i), signal_flow_pos_movie(:, i)] = signalflow_modules(EC_movie, i, 'nomodule');
end
%%
% ----------------------------- Module Averaging with Raw Values ----------------------------- %

% Load the 26 modules atlas
modules = importdata(fullfile(data_path, 'module27(ordered).mat'));
unique_modules = unique(modules);

% Initialize variables to store average and raw signal flow values for each module
avg_signal_flow_pos_rest = zeros(length(unique_modules), 1);
avg_signal_flow_neg_rest = zeros(length(unique_modules), 1);
avg_signal_flow_pos_pain = zeros(length(unique_modules), 1);
avg_signal_flow_neg_pain = zeros(length(unique_modules), 1);
avg_signal_flow_pos_movie = zeros(length(unique_modules), 1);
avg_signal_flow_neg_movie = zeros(length(unique_modules), 1);

raw_signal_flow_pos_rest = cell(length(unique_modules), 1);
raw_signal_flow_neg_rest = cell(length(unique_modules), 1);
raw_signal_flow_pos_pain = cell(length(unique_modules), 1);
raw_signal_flow_neg_pain = cell(length(unique_modules), 1);
raw_signal_flow_pos_movie = cell(length(unique_modules), 1);
raw_signal_flow_neg_movie = cell(length(unique_modules), 1);


% Calculate the average and store raw signal flow values for each module
for module = 1:length(unique_modules)
    
    % Extract indices of each module
    moduleIndx = (modules == unique_modules(module));
    
    % Extract raw signal flow values for each state
    raw_signal_flow_pos_rest{module} = signal_flow_pos_rest(:, moduleIndx);
    raw_signal_flow_neg_rest{module} = signal_flow_neg_rest(:, moduleIndx);
    
    raw_signal_flow_pos_pain{module} = signal_flow_pos_pain(:, moduleIndx);
    raw_signal_flow_neg_pain{module} = signal_flow_neg_pain(:, moduleIndx);
    
    raw_signal_flow_pos_movie{module} = signal_flow_pos_movie(:, moduleIndx);
    raw_signal_flow_neg_movie{module} = signal_flow_neg_movie(:, moduleIndx);
    
    % Compute average signal flow values for each state
    avg_signal_flow_pos_rest(module) = mean(raw_signal_flow_pos_rest{module}(:));
    avg_signal_flow_neg_rest(module) = mean(raw_signal_flow_neg_rest{module}(:));
    
    avg_signal_flow_pos_pain(module) = mean(raw_signal_flow_pos_pain{module}(:));
    avg_signal_flow_neg_pain(module) = mean(raw_signal_flow_neg_pain{module}(:));
    
    avg_signal_flow_pos_movie(module) = mean(raw_signal_flow_pos_movie{module}(:));
    avg_signal_flow_neg_movie(module) = mean(raw_signal_flow_neg_movie{module}(:));
    
end

% Calculate the difference between rest vs. movie and movie vs. rest
diff_avg_signal_flow_pos_rest_movie =  avg_signal_flow_pos_movie - avg_signal_flow_pos_rest;
diff_avg_signal_flow_neg_rest_movie = abs(avg_signal_flow_neg_movie) - abs(avg_signal_flow_neg_rest);

% Calculate the difference between rest vs. pain and pain vs. rest
diff_avg_signal_flow_pos_rest_pain = avg_signal_flow_pos_pain - avg_signal_flow_pos_rest;
diff_avg_signal_flow_neg_rest_pain = abs(avg_signal_flow_neg_pain) - abs(avg_signal_flow_neg_rest);

%%
% ----------------------------- Delta line plot ----------------------------- %
figure;

% Example data (Replace with your vectors)
vector1 = diff_avg_signal_flow_neg_rest_pain; 
vector2 = diff_avg_signal_flow_pos_rest_pain; 

% Normalize the data to retain negative values
max_abs_neg = max(abs(vector1));
normalized_neg = vector1 / max_abs_neg;

% Define thresholds for 1.645 and 1.96
threshold1 = 1.645;
threshold2 = 1.96;

% Compute mean and standard deviation for normalized_neg
mu = mean(normalized_neg);
sigma = std(normalized_neg);
threshold_low_neg1 = mu - threshold1 * sigma;  % Lower threshold for 1.645
threshold_high_neg1 = mu + threshold1 * sigma; % Upper threshold for 1.645
threshold_low_neg2 = mu - threshold2 * sigma;  % Lower threshold for 1.96
threshold_high_neg2 = mu + threshold2 * sigma; % Upper threshold for 1.96

max_abs_pos = max(abs(vector2));
normalized_pos = vector2 / max_abs_pos;

mu = mean(normalized_pos);
sigma = std(normalized_pos);
threshold_low_pos1 = mu - threshold1 * sigma;  % Lower threshold for 1.645
threshold_high_pos1 = mu + threshold1 * sigma; % Upper threshold for 1.645
threshold_low_pos2 = mu - threshold2 * sigma;  % Lower threshold for 1.96
threshold_high_pos2 = mu + threshold2 * sigma; % Upper threshold for 1.96

% --- Compute p-values based on normalized_neg and normalized_pos --- %
% Compute z-scores for each data point
z_neg = (vector1 - 0) / std(vector1);
z_pos = (vector2 - 0) / std(vector2);

% Convert z-scores to p-values (two-tailed test)
p_neg = 2 * (1 - normcdf(abs(z_neg)));
p_pos = 2 * (1 - normcdf(abs(z_pos)));

% False Discovery Rate (FDR) Correction using Benjamini-Hochberg procedure
[~, ~, ~, fdr_p_neg]=fdr_bh(p_neg);
[~, ~, ~, fdr_p_pos]=fdr_bh(p_pos);

% Define FDR significance threshold (typically alpha = 0.05)
alpha = 0.05;

% Add spacing between the curves (joy plot effect)
spacing = 1.5; % Adjust this value to change the space between lines
normalized_pos = normalized_pos + spacing;

xData = 1:27;
num_points = 27;

% Interpolation
x_interp = linspace(min(xData), max(xData), 500); % Increase the number of points for smoother curves
normalized_neg_interp = interp1(xData, normalized_neg, x_interp, 'linear');
normalized_pos_interp = interp1(xData, normalized_pos, x_interp, 'linear');

% Convert hex colors to RGB
% color_above = [184, 184, 255] / 255; % "ED6a5a" pain
color_above = [155, 193, 189] / 255; % "9bc1bd" movie
color_below =  [237, 106, 90] / 255;% "b8b8ff"


% Plot the lines of the vectors
hold on;
plot(xData, normalized_neg, '-o', 'Color', [0.0824, 0.3765, 0.5098], 'LineWidth', 1.5); %
plot(xData, normalized_pos, '-o', 'Color', [0.7373, 0.0078, 0.0078],  'LineWidth', 1.5); %

% Add dashed lines at 0 value for the first line and spacing for the second line
plot(x_interp, zeros(size(x_interp)), '--k', 'LineWidth', 1); % 'Color', [0.0824, 0.3765, 0.5098],
plot(x_interp, spacing * ones(size(x_interp)), '--k', 'LineWidth', 1); % 'Color', [0.7373, 0.0078, 0.0078], 

hold off;

% Add symbols based on thresholds for significance
for i = 1:length(xData)
    % Threshold 1.645 (*)
    if (normalized_neg(i) < threshold_low_neg1 || normalized_neg(i) > threshold_high_neg1) && ...
            ~(normalized_neg(i) < threshold_low_neg2 || normalized_neg(i) > threshold_high_neg2) % Exclude those in 1.96 range
        text(xData(i), normalized_neg(i) + 0.1, '*', 'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 14);
    end
    if (normalized_pos(i) < threshold_low_pos1 + spacing || normalized_pos(i) > threshold_high_pos1 + spacing) && ...
            ~(normalized_pos(i) < threshold_low_pos2 + spacing || normalized_pos(i) > threshold_high_pos2 + spacing)
        text(xData(i), normalized_pos(i) + 0.1, '*', 'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 14);
    end

    % Threshold 1.96 (**)
    if (normalized_neg(i) < threshold_low_neg2 || normalized_neg(i) > threshold_high_neg2) && ...
            fdr_p_neg(i) >= alpha
        text(xData(i), normalized_neg(i) + 0.1, '**', 'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 14);
    end
    if normalized_pos(i) < threshold_low_pos2 + spacing || normalized_pos(i) > threshold_high_pos2 + spacing && ...
            fdr_p_pos(i) >= alpha
        text(xData(i), normalized_pos(i) + 0.1, '**', 'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 14);
    end
    
    % Threshold for FDR 0.05(***)
    if fdr_p_neg(i) < alpha
        text(xData(i), normalized_neg(i) + 0.1, '***', 'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 14);
    end
    if fdr_p_pos(i) < alpha
        text(xData(i), normalized_pos(i) + 0.1, '***', 'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 14);
    end
end

% Customize the plot to give an impression like a joy plot
set(gca, 'XColor', 'none'); % Hide X axis
set(gca, 'YColor', 'none'); % Hide Y axis
box off;
hold off;
set(gca, 'Visible', 'off');  % Hide the axes
% savefig('diff_rest_movie.fig')
