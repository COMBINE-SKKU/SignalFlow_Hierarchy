clear; close all; clc

% Define paths
data_path = '../../data';
data_file_EC = fullfile(data_path, 'MMP_resting_iEC.mat');

% Import Data
EC = importdata(data_file_EC);

% Define the file path for the cortical types data
data_file_corTypes = fullfile(data_path, 'cortical_types.mat');

% Import cortical types data from the specified file
corticalTypes = importdata(data_file_corTypes);

%% Supple impact of negative connections

% compute hierarchical levels (intact)
hierarchyIntact = computeHierarchyLevels(EC, 0.15);

% compute hierarchical levels (positive only)
hierarchyPositive = computeHierarchyLevels(max(EC,0), 0.15);

num_regions = 6;
cyto_label = {'Konicortex', 'Eulaminate-III', 'Eulaminate-II', 'Eulaminate-I', 'Dysgranular', 'Agranular'};

% Initialize variables for table data
data = [];
cyto_labels_for_table = [];
source_labels = [];

% Compute mean values for cyto groups
Intact_data = cell(1, num_regions);
Positive_data = cell(1, num_regions);

for i = 1:num_regions
    % Data for Intact
    val_cyto_intact = hierarchyIntact(corticalTypes == i);
    Intact_data{i} = val_cyto_intact;

    % Data for Positive only
    val_cyto_positive = hierarchyPositive(corticalTypes == i);
    Positive_data{i} = val_cyto_positive;

    % Append intact data to the table arrays
    data = [data; val_cyto_intact];  % First column: data
    cyto_labels_for_table = [cyto_labels_for_table; repmat(cyto_label(i), length(val_cyto_intact), 1)];  % Second column: cyto_label
    source_labels = [source_labels; repmat({'Intact'}, length(val_cyto_intact), 1)];  % Third column: 'Intact'
    
    % Append positive-only data to the table arrays
    data = [data; val_cyto_positive];  % First column: data
    cyto_labels_for_table = [cyto_labels_for_table; repmat(cyto_label(i), length(val_cyto_positive), 1)];  % Second column: cyto_label
    source_labels = [source_labels; repmat({'Positive Only'}, length(val_cyto_positive), 1)];  % Third column: 'Positive Only'
end

% Create the table
result_table = table(data, cyto_labels_for_table, source_labels, ...
    'VariableNames', {'Data', 'CytoLabel', 'Source'});
% save('Supple9.mat', 'result_table');

%% Paired t-tests and Bonferroni correction for multiple comparisons

% Initialize variables to store p-values for t-tests
p_values = zeros(1, num_regions);
median_differences = zeros(1, num_regions);  % Store the median differences

% Perform paired t-test for each cyto_label
for i = 1:num_regions
    val_cyto_intact = Intact_data{i};
    val_cyto_positive = Positive_data{i};
    
    % Ensure both vectors have the same length for paired t-test
    len = min(length(val_cyto_intact), length(val_cyto_positive));
    val_cyto_intact = val_cyto_intact(1:len);
    val_cyto_positive = val_cyto_positive(1:len);
    
    % Perform paired t-test
    [~, p] = ttest(val_cyto_intact, val_cyto_positive);
    p_values(i) = p;
    
    % Calculate the median for Intact and Positive Only data
    median_intact = median(Intact_data{i});
    median_positive = median(Positive_data{i});
    
    % Calculate the difference
    median_differences(i) = median_intact - median_positive;
    end

% Multiple comparison correction using Bonferroni correction
alpha = 0.05;
corrected_alpha = alpha / num_regions;  % Bonferroni correction factor

% Apply correction
significant_indices = p_values < corrected_alpha;
