function plotCytoGroups(hierarchy_data, corticalTypes, cyto_label, type, num_regions)
% plotCytoGroups - Plots a boxplot of cytoarchitectonic groups with scatter points.
% 
% Inputs:
%   hierarchy_data - Numeric data related to cortical hierarchy.
%   corticalTypes  - Vector indicating the cortical type for each region.
%   cyto_label     - Cell array of labels for cytoarchitectonic groups.
%   type           - String indicating the colormap type ('campbell', 'economo', 'tpl').
%   num_regions    - Number of regions to be plotted.
%
% This function creates a boxplot of cytoarchitectonic groups, with the data 
% categorized based on cortical types. It also overlays scatter points colored 
% according to the selected colormap.

% Initialize arrays for storing median values and box plot data
median_values = zeros(1, num_regions);
box_data_cell = cell(1, num_regions);

% Compute median values and prepare boxplot data for each region
for i = 1:num_regions
    val_cyto = hierarchy_data(corticalTypes == i);
    median_values(i) = median(val_cyto);
    box_data_cell{i} = val_cyto;
end

% Prepare sorted box plot data
box_data = [];
group_labels = [];
sorted_cyto_label = cell(1, num_regions);
group_count = 0; % Counter for the actual groups with data

for i = 1:num_regions
    val_cyto = box_data_cell{i};
    if ~isempty(val_cyto)
        group_count = group_count + 1;
        box_data = [box_data; val_cyto];
        group_labels = [group_labels; repmat(group_count, length(val_cyto), 1)];
        sorted_cyto_label{group_count} = strrep(cyto_label{i}, 'L_', '');
    end
end
sorted_cyto_label = sorted_cyto_label(1:group_count);

% Create the boxplot
figure;
boxplot(box_data, group_labels, 'Labels', sorted_cyto_label, 'Color', 'k', 'symbol', '');
xlabel('Cytoarchitectonic Groups');
ylabel('Values');
title('Box Plot of Cytoarchitectonic Groups');
xtickangle(45);

% Select and load the appropriate colormap
switch type
    case 'campbell'
        clrData = '/combinelab/03_user/younghyun/04_software/colormaps/campbellColormap.mat';
        values = importdata(clrData);
        cyto_label = {'Prefrontal', 'Frontal', 'InterPrecentral', 'Precentral', ...
                      'Postcentral', 'InterPostcentral', 'Parietal', 'Visuopsychic', ...
                      'Visuosensory', 'Temporal', 'Audiopsychic', 'Audiosensory', ...
                      'Olfactory', 'LimbicA', 'LimbicB', 'LimbicC', 'Insula'};
        hex_colors = arrayfun(@(i) sprintf('#%02X%02X%02X', floor(values(i,1:3) * 255)), 2:length(values), 'UniformOutput', false);
        color_map = containers.Map(cyto_label, hex_colors);
    case 'economo'
        color_map = containers.Map({'agranular', 'frontal', 'parietal', 'polar', 'granular'}, ...
                                   {'#7F297F', '#35689B', '#A8D38E', '#FDCD0B', '#F3EC1A'});
    case 'tpl'
        color_map = containers.Map({'Konicortex', 'Eulaminate-III', 'Eulaminate-II', 'Eulaminate-I', ...
                                    'Dysgranular', 'Agranular'}, ...
                                   {'#5880aa', '#8ba7b0', '#abbaa2', '#d9c598', '#fdd3c8', '#fdccfd'});
end

% Add scatter points with jitter to boxplot
hold on;
jitter_amount = 0.3; % Define jitter amount
for i = 1:group_count
    label = sorted_cyto_label{i};
    if isKey(color_map, label)
        rgb_color = hex2rgb(color_map(label)); % Convert hex to RGB
    end
    y = box_data(group_labels == i);
    x = i + (rand(size(y)) - 0.5) * jitter_amount; % Jitter x positions
    scatter(x, y, 20, 'filled', 'MarkerFaceColor', rgb_color, 'MarkerFaceAlpha', 0.8);
end
hold off;

end

% Helper function to convert hex to RGB
function rgb = hex2rgb(hex_str)
    hex_str = char(hex_str(2:end)); % Remove the hash (#)
    rgb = reshape(sscanf(hex_str, '%2x') / 255, 1, 3);
end
