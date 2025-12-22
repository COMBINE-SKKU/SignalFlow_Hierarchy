% Clear workspace
clear; close all; clc;

%% Load MMP module information
% Set path
data_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy/data';
module_csv_file = fullfile(data_dir, 'Glasser_22_modules.xlsx');

% Load module data
module_data = readtable(module_csv_file);

% Extract module names
module_names = module_data.Names;

% Create left and right hemisphere module names
left_modules = strcat('L_', module_names);
right_modules = strcat('R_', module_names);

% Concatenate left and right hemisphere names
module_names = [right_modules;left_modules];

%% Load f-tract data
region_name_path = '/combinelab/03_user/younghyun/03_data/MNI-HCP-MMP1/export/amplitude';
region_name_file = fullfile(region_name_path, 'amplitude.csv');

% Load the CSV file with region names as strings
opts = detectImportOptions(region_name_file);
opts = setvartype(opts, 'char');  % Read all columns as char initially
data = readtable(region_name_file, opts);

% Extract region names from the first column
region_names = data{:,1};

ftract = str2double(data{:,2:end});  % Convert to numeric matrix

% Handle any NaN values that were in the CSV
ftract(isnan(ftract)) = 0;

% Sort matrix based on module_names order
[~, sort_idx] = ismember(module_names, region_names);

% Check if all module names were found
if any(sort_idx == 0)
    warning('Some module names were not found in region names');
    % Display missing modules for debugging
    missing_modules = module_names(sort_idx == 0);
    disp('Missing modules:');
    disp(missing_modules);
end

% Reorder matrix rows and columns according to the sorting index
ftract = ftract(sort_idx, sort_idx)';
region_names = region_names(sort_idx);

%% Compare with iEC
% Load iEC matrix (replace with your actual path to iEC data)
data_dir2 = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy/step1_ec_validation/part3_empirical_human/data';
iEC = importdata(fullfile(data_dir, 'MMP360_resting_iEC.mat'));
FC = importdata(fullfile(data_dir2, 'MMP360_fc_fcd_td_test.mat')).fc_emp_test;

% Calculate column-wise correlations and p-values
column_correlations = zeros(size(ftract, 2), 1);
column_pvalues = zeros(size(ftract, 2), 1);
for i = 1:size(ftract, 2)
    [r, p] = corr(ftract(:,i), iEC(:,i), 'Rows', 'complete');
    column_correlations(i) = r;
    column_pvalues(i) = p;
end

% Determine region with median correlation
medianCorrelation = nanmedian(column_correlations);
[~, medianIdx] = min(abs(column_correlations - medianCorrelation));
disp(['Region with median correlation: ', region_names{medianIdx}]);
disp(['Correlation: ', num2str(column_correlations(medianIdx))]);
disp(['p-value: ', num2str(column_pvalues(medianIdx))]);

%% Permutation test for each parcel
n_parcels = size(ftract, 2);
n_perms = 1000;

% Calculate correlations for each parcel with different measures
iEC_corrs = zeros(n_parcels, 1);
FC_corrs = zeros(n_parcels, 1);

for i = 1:n_parcels
    iEC_corrs(i) = corr(ftract(:,i), iEC(:,i), 'Rows', 'complete');
    FC_corrs(i) = corr(ftract(:,i), FC(:,i), 'Rows', 'complete');
end

% Perform statistical tests to compare correlations
[p_iEC_vs_FC] = signrank(iEC_corrs, FC_corrs, 'tail', 'right');

% Display mean correlations and p-values
fprintf('Mean iEC correlation: r = %.2f ± %.2f\n', nanmean(iEC_corrs), nanstd(iEC_corrs));
fprintf('Mean FC correlation: r = %.2f ± %.2f\n', nanmean(FC_corrs), nanstd(FC_corrs));
fprintf('iEC vs. FC paired t-test: p = %.4f\n', p_iEC_vs_FC);

% Combine all correlations for permutation testing
actual_corrs = iEC_corrs;
p_values = zeros(n_parcels, 1);

% Permutation test for each parcel
parfor i = 1:n_parcels
    null_corrs = zeros(n_perms, 1);
    parcel_data = iEC(:,i);
    target = ftract(:,i);
    
    for j = 1:n_perms
        % Randomly permute the order of one variable
        shuffled_data = parcel_data(randperm(length(parcel_data)));
        null_corrs(j) = corr(shuffled_data, target, 'Rows', 'complete');
    end
    
    % Calculate p-value as fraction of permutations with |correlation| >= |observed|
    p_values(i) = sum(abs(null_corrs) >= abs(iEC_corrs(i))) / n_perms;
end

% Apply FDR correction
[h, ~, ~] = fdr_bh(p_values, 0.05, 'pdep', 'yes');

% Calculate percentage of robust correlations
percent_significant = sum(h) / n_parcels * 100;

fprintf('Permutation test significance: p < %.4f\n', mean(p_values));
fprintf('Percentage of parcels with significant correlations: %.1f%%\n', percent_significant);

%% Figure 1: Boxplot with scatter comparing iEC and FC
figure('Position', [100 100 400 400]);

% Define colors for each box
iEC_color = [0.8 0.6 0.9];
FC_color = [0.5 0.7 0.9];  % Light blue color for FC

% Prepare data for boxplot
all_corrs = [iEC_corrs, FC_corrs];
h = boxplot(all_corrs, 'Colors', 'k', 'Symbol', '', 'Labels', {'iEC', 'FC'});
hold on;

set(h, 'LineWidth', 1.5);
set(findobj(gca, 'type', 'line', 'Tag', 'Median'), 'Color', 'k', 'LineWidth', 2);

% Apply custom box colors using patch approach
for i = 1:2
    % Define color for current box
    if i == 1
        box_color = iEC_color;
    else
        box_color = FC_color;
    end
    
    % Get box data and create patch
    box_data = get(h(5,i), 'XData');
    box_data_y = get(h(5,i), 'YData');
    if ~isempty(box_data) && ~isempty(box_data_y)
        patch(box_data, box_data_y, box_color, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    end
end

% Add jittered data points for each group
for i = 1:2
    % Add jittered data points
    x = i + (randn(n_parcels,1)*0.05);
    scatter(x, all_corrs(:,i), 10, [0.3 0.3 0.3], 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
end

% Add significance annotations
if p_iEC_vs_FC < 0.05
    if p_iEC_vs_FC < 0.001
        sig_str = '***';
    elseif p_iEC_vs_FC < 0.01
        sig_str = '**';
    else
        sig_str = '*';
    end
    line([1, 2], [max(all_corrs(:))*1.05, max(all_corrs(:))*1.05], 'Color', 'k');
    text(1.5, max(all_corrs(:))*1.07, sig_str, 'HorizontalAlignment', 'center');
end

set(gca, 'XTickLabel', {}, 'YTickLabel', {});
set(gca, 'YTick', linspace(min(all_corrs(:)), max(all_corrs(:)), 4));
set(gca, 'YLim', [min(all_corrs(:))*0.95, max(all_corrs(:))*1.15]);
set(gca, 'LineWidth', 1.5);
box on;
grid on
exportgraphics(gcf, 'boxplot_iEC_FC.png', 'Resolution', 1200);

%% Figure 2: iEC matrix visualization - Left cortex only
figure('Position', [100 100 400 400]);
left_idx = 1:180;  % Left cortex indices
left_iEC = iEC(left_idx, left_idx);  % Extract left cortex submatrix
imagesc(left_iEC);  % Show only left cortex
colormap(generateColorMap(left_iEC(:), 1000));
axis square;
set(gca, 'XTickLabel', {}, 'YTickLabel', {}, 'XTick', [], 'YTick', []);
set(gcf, 'Color', 'w');
% exportgraphics(gcf, 'iEC_left_cortex.png', 'Resolution', 1200);

% Figure 3: ftract matrix visualization - Left cortex only
figure('Position', [100 100 400 400]);
left_ftract = ftract(left_idx, left_idx);  % Extract left cortex submatrix
imagesc(left_ftract);  % Show only left cortex
colormap(generateColorMap(left_ftract(:), 1000));
axis square;
set(gca, 'XTickLabel', {}, 'YTickLabel', {}, 'XTick', [], 'YTick', []);
set(gcf, 'Color', 'w');
% exportgraphics(gcf, 'ftract_left_cortex.png', 'Resolution', 1200);

%% Compare iEC original vs. iEC symmetric using simulation
f_diff = importdata(fullfile(data_dir, 'MMP360_peak_freq.mat'));
fc_emp_test = importdata(fullfile(data_dir2, 'MMP360_fc_fcd_td_test.mat')).fc_emp_test;
fcd_emp_test = importdata(fullfile(data_dir2, 'MMP360_fc_fcd_td_test.mat')).fcd_emp_test;
iEC_sym = (iEC + iEC')/2;
TR = 0.72;
Tmax = 2400;

num_iter = 50;
gofs = zeros(num_iter, 2);
triu_ind_fc = find(triu(ones(size(fc_emp_test)),1));
fc_emp_triu = fc_emp_test(triu_ind_fc);


parfor i = 1:num_iter
    bold_sim = run_simulation(iEC, 1, f_diff, TR, Tmax, -0.01)';
    fc_sim = corr(bold_sim);
    fcd_sim = calculate_fcd(bold_sim, TR, ceil(60/TR), 1);
    fc_corr = corr(fc_sim(triu_ind_fc), fc_emp_triu);
    fcd_ks = max(abs(fcd_emp_test - fcd_sim));
    gofs(i,1) = fc_corr - fcd_ks;
end

parfor i = 1:num_iter
    bold_sim = run_simulation(iEC_sym, 1, f_diff, TR, Tmax, -0.01)';
    fc_sim = corr(bold_sim);
    fcd_sim = calculate_fcd(bold_sim, TR, ceil(60/TR), 1);
    fc_corr = corr(fc_sim(triu_ind_fc), fc_emp_triu);
    fcd_ks = max(abs(fcd_emp_test - fcd_sim));
    gofs(i,2) = fc_corr - fcd_ks;
end

save('gofs_iEC_symmetric.mat', 'gofs');

%% Boxplot of GOFS
figure('Position', [100 100 400 400]);

% Define colors for each box
original_color = [0.8 0.6 0.9];
symmetric_color = [0.5 0.7 0.9];

% Prepare data for boxplot
h = boxplot(gofs, 'Colors', 'k', 'Symbol', '', 'Labels', {'Original', 'Symmetric'});
hold on;

% Apply custom box colors using patch approach
for i = 1:2
    % Define color for current box
    if i == 1
        box_color = original_color;
    else
        box_color = symmetric_color;
    end
    
    % Get box data and create patch
    box_data = get(h(5,i), 'XData');
    box_data_y = get(h(5,i), 'YData');
    if ~isempty(box_data) && ~isempty(box_data_y)
        patch(box_data, box_data_y, box_color, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    end
end

% Add jittered data points for each group
for i = 1:2
    % Add jittered data points
    x = i + (randn(size(gofs,1),1)*0.05);
    scatter(x, gofs(:,i), 10, [0.3 0.3 0.3], 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
end

% Perform statistical test
[p_val] = signrank(gofs(:,1), gofs(:,2), 'tail', 'both');

% Add significance annotations
if p_val < 0.05
    if p_val < 0.001
        sig_str = '***';
    elseif p_val < 0.01
        sig_str = '**';
    else
        sig_str = '*';
    end
    line([1, 2], [max(gofs(:))*1.05, max(gofs(:))*1.05], 'Color', 'k');
    text(1.5, max(gofs(:))*1.07, sig_str, 'HorizontalAlignment', 'center');
end

% set(gca, 'XTickLabel', {}, 'YTickLabel', {});
set(gca, 'YTick', linspace(min(gofs(:)), max(gofs(:)), 4));
set(gca, 'YLim', [min(gofs(:))*0.95, max(gofs(:))*1.15]);
set(gca, 'LineWidth', 1.5);
box on;
grid on
% exportgraphics(gcf, 'boxplot_gofs_iEC_symmetric.png', 'Resolution', 1200);

