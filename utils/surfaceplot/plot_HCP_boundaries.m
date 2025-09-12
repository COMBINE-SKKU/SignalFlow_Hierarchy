% plot boundaries of selected areas on top of fct harmonics (on the anatomy)
% border files were downloaded from BALSA and converted to metric files
% using wb_command as follows ("\\" means nothing):
% wb_command -border-to-vertices \\
% /home/localadmin/manifolds/surface_plots/S900.R.inflated_MSMAll.32k_fs_LR.surf.gii \\
% /home/localadmin/manifolds/surface_plots/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil.32k_fs_LR.border\\
% borders_right.func.gii
% link to documentation for this function: https://www.humanconnectome.org/software/workbench-command/-border-to-vertices
% files for somatotopic sub-areas:
% Q1-Q6_RelatedParcellation210.?.SubAreas.32k_fs_LR.border (replace '?'
% with 'L' and 'R' for left and right hemisphere
% more areas/subsets of areas can easily be implemented (see below)
% code by Katharina Glomb
% katharina.glomb@gmail.com

function plot_data_with_borders = plot_HCP_boundaries(plot_data_without_borders,clr,atlas)

addpath(genpath('/local_raid1/03_user/younghyun/02_data/parcellations'));

switch atlas
    case "MMP360"
        load('MMP_32k_very_inflated.mat','borders_left','borders_right');
    case "Schaefer100"
        load('borders_Schaefer100.mat','borders_left','borders_right');
    case "Markov"
        load('borders_Markov.mat','borders_left','borders_right');
end



plot_borders_left = any(borders_left,2);
plot_borders_right = any(borders_right,2);
plot_all_borders = [plot_borders_left;plot_borders_right];

plot_data_with_borders = plot_data_without_borders;

if size(plot_data_without_borders,2) == 1
    plot_data_with_borders(plot_all_borders) = max(plot_data_without_borders)+(max(plot_data_without_borders)-min(plot_data_without_borders))/65;
% explicitly assign colors 
elseif size(plot_data_without_borders,2) == 3
    if isstr(clr)
        clr = rgb(clr);
    end
    plot_data_with_borders(plot_all_borders,:) = repmat(clr,[sum(plot_all_borders),1]);
end