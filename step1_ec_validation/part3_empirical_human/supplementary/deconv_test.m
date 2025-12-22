clear; close all; clc;

main_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy/step1_ec_validation/part3_empirical_human';

ec_orig = importdata(fullfile(main_dir, 'data', 'Schaefer100_ec_results.mat')).ec_results;

ec_deconv = importdata(fullfile(main_dir, 'run_env_nike', 'ec_results_deconv_Schaefer100.mat')).ec_results;

% Calculate group averages for ec_deconv
deconv_avg_var = mean(cat(3, ec_deconv.var), 3);
deconv_avg_fask = mean(cat(3, ec_deconv.fask), 3);
deconv_avg_lingam = mean(cat(3, ec_deconv.lingam), 3);

% Calculate group averages for ec_orig (using indices 1:20 and taking 33:end,33:end)
var_cell = cell(20,1);
fask_cell = cell(20,1);
lingam_cell = cell(20,1);

for i = 1:20
    var_cell{i} = ec_orig(i).var(33:end,33:end);
    fask_cell{i} = ec_orig(i).fask(33:end,33:end);
    lingam_cell{i} = ec_orig(i).lingam(33:end,33:end);
end

orig_avg_var = mean(cat(3, var_cell{:}), 3);
orig_avg_fask = mean(cat(3, fask_cell{:}), 3);
orig_avg_lingam = mean(cat(3, lingam_cell{:}), 3);

% Calculate correlations
corr_var = corr(deconv_avg_var(:), orig_avg_var(:));
corr_fask = corr(deconv_avg_fask(:), orig_avg_fask(:));
corr_lingam = corr(deconv_avg_lingam(:), orig_avg_lingam(:));

% Display results
fprintf('Correlation between deconv and orig:\n');
fprintf('VAR: %.4f\n', corr_var);
fprintf('FASK: %.4f\n', corr_fask);
fprintf('LINGAM: %.4f\n', corr_lingam);


% Beta values
betas = [7, 12, 0.5];
orig_iec = betas(1) * orig_avg_var + betas(2) * orig_avg_fask + betas(3) * orig_avg_lingam;
deconv_iec = betas(1) * deconv_avg_var + betas(2) * deconv_avg_fask + betas(3) * deconv_avg_lingam;
orig_iec = orig_iec(:)/max(abs(orig_iec(:)));
deconv_iec = deconv_iec(:)/max(abs(deconv_iec(:)));
corr_iec = corr(orig_iec(:), deconv_iec(:));

fprintf('Correlation between orig and deconv:\n');
fprintf('IEC: %.4f\n', corr_iec);

%% Scatter plot between orig and deconv iEC
figure('Position', [100 100 400 300]);
scatter(orig_iec(:), deconv_iec(:), 40, [0.4 0.4 0.4], 'filled', ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);

% Add regression line
hold on;
p = polyfit(orig_iec(:), deconv_iec(:), 1);
x_fit = linspace(min(orig_iec(:)), max(orig_iec(:)), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, '--', 'Color', 'r', 'LineWidth', 1.5);

% Customize plot
box on;
grid on;

% Set axis limits with some padding separately for x and y
x_min = min(orig_iec(:)) * 1.1;
x_max = max(orig_iec(:)) * 1.1;
y_min = min(deconv_iec(:)) * 1.1;
y_max = max(deconv_iec(:)) * 1.1;
xlim([x_min x_max]);
ylim([y_min y_max]);

% Customize appearance
set(gca, 'LineWidth', 1.5, ...
    'XTick', linspace(x_min, x_max, 5), ...
    'YTick', linspace(y_min, y_max, 5), ...
    'TickLength', [0.02 0.02]);

% % Save figure
% exportgraphics(gcf, fullfile(main_dir, 'run_env_nike', 'deconv_test.png'), ...
%     'Resolution', 1200, 'BackgroundColor', 'white');
