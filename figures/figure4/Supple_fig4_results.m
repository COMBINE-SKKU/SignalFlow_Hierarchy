clear; close all; clc

% Define paths
data_path = '/combinelab/03_user/younghyun/01_project/01_HierarchyMapping/data';
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
save('Supple9.mat', 'result_table');

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

%% Supple sensory hierarchy
% Compute hierarchical levels
hierarchyLevels = computeHierarchyLevels(EC, 0.15);

% Define visual areas and their names
Vis = [181, 186, 184, 185, 193, 199, 332, 197, 196, 183, ...
       187, 198, 202, 343, 333, 340, 334, 337, 339, 200 ,201, ...
       318, 338, 336, 182, 203];
VisNames = {'V1', 'V4', 'V2', 'V3', 'V3A', 'V3B', 'V6A', 'IPS1', 'V7', 'V6', ...
            'V8', 'FFC', 'PIT', 'VVC', 'VMV1', 'VMV2', 'VMV3', 'FST', 'LO3', 'LO1', ...
            'LO2', 'PH', 'V3CD', 'V4t', 'MST', 'MT'};

% Initialize hierarchy levels for visualization
hierarchyVis = zeros(360, 1);
hierarchyVis(Vis) = hierarchyLevels(Vis);

% Extract the hierarchical levels for the visual areas and sort them
[hierarchyLevelsSorted, sortIdx] = sort(hierarchyLevels(Vis));

% Sort VisNames according to the sorted hierarchical levels
VisNamesSorted = VisNames(sortIdx);

% Create bar plot of hierarchical levels in ascending order
figure;
bar(hierarchyLevelsSorted);
set(gca, 'XTickLabel', VisNamesSorted, 'XTick', 1:length(VisNamesSorted));
xlabel('Visual Areas');
ylabel('Hierarchical Levels');
title('Hierarchical Levels of Visual Areas (Ascending Order)');
xtickangle(45); % Rotate x-axis labels for better readability
grid on;

surfaceplot(hierarchyVis,'MMP','both','viridis2')

% Define auditory areas and their names
AuditoryAreas = [283, 285, 354, 353, 304, 204, 284, 303, 308, 309, 356, 310, 287, 355, 305];
AuditoryNames = {'52', 'PFcm', 'LBelt', 'MBelt', 'PBelt', 'A1', 'RI', 'STGa', 'STSda', ...
                 'STSdp', 'STSva', 'STSvp', 'TA2', 'A4', 'A5'};

% Initialize hierarchy levels for visualization
hierarchyAuditory = zeros(360, 1);
hierarchyAuditory(AuditoryAreas) = hierarchyLevels(AuditoryAreas);

% Extract the hierarchical levels for the auditory areas and sort them
[hierarchyLevelsSorted, sortIdx] = sort(hierarchyLevels(AuditoryAreas));

% Sort AuditoryNames according to the sorted hierarchical levels
AuditoryNamesSorted = AuditoryNames(sortIdx);

% Create bar plot of hierarchical levels in ascending order
figure;
bar(hierarchyLevelsSorted);
set(gca, 'XTickLabel', AuditoryNamesSorted, 'XTick', 1:length(AuditoryNamesSorted));
xlabel('Auditory Areas');
ylabel('Hierarchical Levels');
title('Hierarchical Levels of Auditory Areas (Ascending Order)');
xtickangle(45); % Rotate x-axis labels for better readability
grid on;

surfaceplot(hierarchyAuditory,'MMP','both','viridis2')