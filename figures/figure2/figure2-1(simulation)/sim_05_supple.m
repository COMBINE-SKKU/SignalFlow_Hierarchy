%% Supplementary Fig. 3 Panel A
% This script generates a line plot showing the evolution of beta values
% of each algorithm across iterations.
clc; close all; clear;

% Load beta histories data
betaHistories = importdata(fullfile('fig2_results/EC_beta_histories.mat'));

% Select the atlas
atlas = 1;
betaHistoriesAtlas = betaHistories{atlas};

% Extract beta values for each algorithm
rDCM_beta = squeeze(betaHistoriesAtlas(:, 1, :))';
VAR_beta = squeeze(betaHistoriesAtlas(:, 2, :))';
FASK_beta = squeeze(betaHistoriesAtlas(:, 3, :))';

% Define the scale factor for error bars
scaleFactor = 0.1;

% Define optimization steps
optim_step = 1:50;

% Create figure
figure;

% Plot rDCM beta values
rDCMColor = [1, 0.89, 0.64];
shadedErrorBar(optim_step, rDCM_beta, {@mean, @(x) std(x) * scaleFactor}, ...
               'lineprops', {'-k', 'Color', rDCMColor}, ...
               'transparent', true, ...
               'patchSaturation', 0.5);

% Plot VAR beta values
VARColor = [0.98, 0.86, 0.88];
shadedErrorBar(optim_step, VAR_beta, {@mean, @(x) std(x) * scaleFactor}, ...
               'lineprops', {'-k', 'Color', VARColor}, ...
               'transparent', true, ...
               'patchSaturation', 0.5);

% Plot FASK beta values
FASKColor = [0.61, 0.77, 0.86];
shadedErrorBar(optim_step, FASK_beta, {@mean, @(x) std(x) * scaleFactor}, ...
               'lineprops', {'-k', 'Color', FASKColor}, ...
               'transparent', true, ...
               'patchSaturation', 0.5);

% Set y-axis limits
ylim([0 50]);

% Add legend
legend({'rDCM', 'VAR', 'FASK'}, 'Location', 'Best');

%% Supplementary Fig. 3 Panel B
% This script generates images of matrices and histograms of edge strength
% distribution for structural connectivity matrices of each atlas.
clc; close all; clear;

% List of atlases
atlases = {'Schaefer', 'MMP'};

% Define data path
data_path = 'data';

% set atlas variable
atlas = 1;
atlas_name = atlases{atlas};

% Load structural connectivity matrix
SC = importdata(fullfile(data_path, sprintf('%s_SC.mat', atlas_name)));

% Normalize the structural connectivity matrix
SC_normalized = SC / max(SC(:));
kurto = kurtosis(nonzeros(SC_normalized(:)));

% Create a figure and histogram for the edge strength distribution
figure;
histogram(nonzeros(SC_normalized(:)), 100, 'FaceColor', '#808080', 'Normalization', 'probability');

% Set the figure and axes properties
set(gcf, 'Position', [100, 100, 225, 225]);
set(gca, 'XTick', [0 0.4 0.8], 'YTick', [0 0.2 0.4 0.6]);
set(gca, 'FontName', 'Helvetica', 'FontSize', 10);
xlabel('Edge Strength');
ylabel('Probability');
savefig(sprintf('SC_%s.fig', atlas_name));

% Apply log transformation to the structural connectivity matrix
SC_transformed = log(SC_normalized + 0.1) + 100;

% Create an image of the transformed structural connectivity matrix
figure;
imagesc(SC_transformed);
colormap(generateColorMap(SC_transformed(:), 1000));
axis square;
set(gca, 'XTick', [], 'YTick', []);
savefig(sprintf('SC_%s_mat.fig', atlas_name));
