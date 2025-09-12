function surfaceplot(source,atlas,hemi,colormaptype)

% set path
dir_atlas = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy/utils/surfaceplot/data';
addpath(genpath('/local_raid1/01_software/toolboxes/toolboxes/cifti-matlab'))
switch atlas
    case "MMP360"
    % MMP 360
    aname = fullfile(dir_atlas, '/Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
    rois = ft_read_cifti(aname).indexmax;
    case "Schaefer100"
    %Schaefer 100
    aname = fullfile(dir_atlas,'/Schaefer2018_100Parcels_7Networks_order.dlabel.nii');
    rois = ft_read_cifti(aname).parcels;  
end

% Load Atlas
atlasroi = [ gifti(fullfile(dir_atlas, '100206.R.atlasroi.32k_fs_LR.shape.gii')).cdata; ...
        gifti(fullfile(dir_atlas, '100206.R.atlasroi.32k_fs_LR.shape.gii')).cdata];

surf_type = 'very_inflated_MSMAll';
file_name.left  = ['100206.L.', surf_type, '.32k_fs_LR.surf.gii'];
file_name.right = ['100206.R.', surf_type, '.32k_fs_LR.surf.gii'];

temp = gifti(fullfile(dir_atlas, file_name.left));
vertices.left   = temp.vertices;
faces.left      = temp.faces;
clear temp;

temp = gifti(fullfile(dir_atlas, file_name.right));
vertices.right   = temp.vertices;
faces.right      = temp.faces;

vertices.all    = [vertices.left; vertices.right];
faces.all       = [faces.left; faces.right+size(vertices.left,1)];

Cmap = zeros(64984,1);

switch atlas
    case "MMP360"
        for i =1:360
            ind = rois == i;
            Cmap(ind,1) = source(i);
        end
    case "Schaefer100"
        for i = 1:100
            ind = rois == i;
            Cmap(ind,1) = source(i);
        end
end

Cmap(atlasroi == 0) = 0;

n = 100;
switch colormaptype
    case 'viridis'
        clrs = viridis(n);
        n = size(clrs,1);
    case 'viridis2'
        clrs = viridis2(n);    
    case 'BlRd'
        clrs = generateColorMap(Cmap(:),n);
%         clrs = flipud(redblue);
    case 'Yeo7'
        clrs = importdata('../../data/7NetworksColors.mat')/255;
        n = 8;
    case 'inferno'
        clrs = slanCM('inferno');
        n = size(clrs,1);
    case 'parula'
        clrs = colormap('parula');
        n = size(clrs,1);
    case 'economo'
        clrs = EKcolormap(6);
        n = size(clrs,1);
    case 'tpl'
        clrs = tplcolormap(7);
        n = size(clrs,1);
    case 'campbell'
        clrs = Campbell(18);
        n = size(clrs,1);
    case 'mesulam'
        clrs = mesulamclr(5);
        n = size(clrs,1);        
    case 'pc1'
        clrs = pc1color;
        n = size(clrs,1);
    case 'module27'
        clrs_temp = importdata('../../data/module27(ordered)clrs.mat');
        clrs = zeros(28,3);
        clrs(1,:) = [0.0039, 0.0039, 0.0039];
        clrs(2:end,:) = clrs_temp;
        n=28;
    case 'module22'
        clrs = MMPcolor;
        n = 23;
    case 'zone4'
        clrs = importdata('../../data/cortical_zones_clrs.mat');
        n=5;
end

clr_IDs = round((Cmap-min(Cmap))/max(Cmap-min(Cmap))*(n-1))+1;
plot_colors = clrs(clr_IDs,:);
medial_wall = [128, 128, 128]./255;

% Find the indices where atlasroi is 0
medial_indices = find(atlasroi == 0);

% Fill those indices in plot_colors with the color values given by medial_wall
plot_colors(medial_indices, :) = repmat(medial_wall, length(medial_indices), 1);


plot_data_mf_w_borders = plot_HCP_boundaries(plot_colors,'black',atlas);
switch hemi
    case 'both'

        % Set figure
        f = figure;
        f.Units='centimeters';
        f.Position = [22, 22, 6, 7.54];

        % Options for subtightplot
        gapleft = 0.0000005;
        gapright = 0.0000005;
        margin_height = [0.05, 0.05]; % lower, upper
        margin_width = [0.0000005, 0.0000005];  % left, right
        
        nRows = 3;
        nCols = 2;

        % left lateral
        subtightplot(nRows,nCols, 1,[gapleft gapright], margin_height, margin_width);
        options.face_vertex_color = plot_data_mf_w_borders(1:length(plot_data_mf_w_borders)/2,:);
        plot_mesh(vertices.left, faces.left, options);
        view([270 0])
        material('dull')
        camlight(190,-180)
        
        % right lateral
        subtightplot(nRows,nCols, 2,[gapleft gapright], margin_height, margin_width);
        options.face_vertex_color = plot_data_mf_w_borders(length(plot_data_mf_w_borders)/2+1:end,:);
        plot_mesh(vertices.right, faces.right, options);
        view([90 0])
        material('dull')
        camlight(145,215)

        % left medial
        subtightplot(nRows,nCols, 3,[gapleft gapright], margin_height, margin_width);
        options.face_vertex_color = plot_data_mf_w_borders(1:length(plot_data_mf_w_borders)/2,:);
        plot_mesh(vertices.left, faces.left, options);
        view([90 0])
        material('dull')
        camlight(-190,240)

        % right medial
        subtightplot(nRows,nCols, 4,[gapleft gapright], margin_height, margin_width);
        options.face_vertex_color = plot_data_mf_w_borders(length(plot_data_mf_w_borders)/2+1:end,:);
        plot_mesh(vertices.right, faces.right, options);
        view([270 0])
        material('dull')
        camlight(-190,-180)
        
         if min(Cmap)<0 && max(Cmap)>0
            subtightplot(nRows,nCols,5,[0.04 0.04], [0.025 0], [0.05 0.01]);
            colormap(clrs)
            axis off
            h=colorbar('Location', 'SouthOutside');
            set(h,'XTick',[0,abs(min(Cmap))/(abs(min(Cmap))+abs(max(Cmap))),1])
            % set labels
            xlabs = [min(Cmap),0,max(Cmap)];
            set(h,'XTickLabel',xlabs)
            set(h, 'Position', [0.2495    0.2853    0.4815    0.0468]);
            set(h, 'AxisLocation' ,'in');
         elseif min(Cmap)<0 && max(Cmap)<=0
            subtightplot(nRows,nCols,5,[0.04 0.04], [0.025 0], [0.05 0.01]);
            colormap(clrs)
            axis off
            h=colorbar('Location', 'SouthOutside');
            set(h,'XTick',[0,abs(min(Cmap))/(abs(min(Cmap))+abs(max(Cmap))),1])
            % set labels
            xlabs = [min(Cmap),0,max(Cmap)];
            set(h,'XTickLabel',xlabs)
            set(h, 'Position', [0.2495    0.2853    0.4815    0.0468]);
            set(h, 'AxisLocation' ,'in');
         else
            subtightplot(nRows,nCols,5,[0.04 0.04], [0.025 0], [0.05 0.01]);
            colormap(clrs)
            axis off
            h=colorbar('Location', 'SouthOutside');
            set(h,'XTick',[0,1])
            xlabs = [min(Cmap),max(Cmap)];
            set(h,'XTickLabel',xlabs)
            set(h, 'Position', [0.2495    0.2853    0.4815    0.0468]);
            set(h, 'AxisLocation' ,'in');
         end
    
             
    case 'single'

        % Set figure
        f = figure;
        f.Units='centimeters';
        f.Position = [22, 22, 12.5, 3.5];

        % Options for subtightplot
        gapleft = 0.0000005;
        gapright = 0.0000005;
        margin_height = [0.05, 0.05]; % lower, upper
        margin_width = [0.0000005, 0.0000005];  % left, right
        
        nRows = 1;
        nCols = 3;

        % left lateral
        subtightplot(nRows,nCols, 1,[gapleft gapright], margin_height, margin_width);
        options.face_vertex_color = plot_data_mf_w_borders(1:length(plot_data_mf_w_borders)/2,:);
        plot_mesh(vertices.left, faces.left, options);
        view([270 0])
        material('dull')
        camlight(190,-180)

        % left medial
        subtightplot(nRows,nCols, 2,[gapleft gapright], margin_height, margin_width);
        options.face_vertex_color = plot_data_mf_w_borders(1:length(plot_data_mf_w_borders)/2,:);
        plot_mesh(vertices.left, faces.left, options);
        view([90 0])
        material('dull')
        camlight(-190,240)
        
         if min(Cmap)<0 && max(Cmap)>0
            subtightplot(nRows,nCols,3,[0.04 0.04], [0.025 0], [0.05 0.01]);
            colormap(clrs)
            axis off
            h=colorbar('Location', 'East');
            set(h,'XTick',[0,abs(min(Cmap))/(abs(min(Cmap))+abs(max(Cmap))),1])
            % set labels
            xlabs = [min(Cmap),0,max(Cmap)];
            set(h,'XTickLabel',xlabs)
            set(h, 'Position', [0.6822    0.0842    0.0287    0.8717]);
            set(h, 'AxisLocation' ,'in');
        elseif min(Cmap)<0 && max(Cmap)<=0
            subtightplot(nRows,nCols,3,[0.04 0.04], [0.025 0], [0.05 0.01]);
            colormap(clrs)
            axis off
            h=colorbar('Location', 'East');
            set(h,'XTick',[0 1])
            % set labels
            xlabs = [min(Cmap), 0];
            set(h,'XTickLabel',xlabs)
            set(h, 'Position', [0.6822    0.0842    0.0287    0.8717]);
            set(h, 'AxisLocation' ,'in');
         else
            subtightplot(nRows,nCols,3,[0.04 0.04], [0.025 0], [0.05 0.01]);
            colormap(clrs)
            axis off
            h=colorbar('Location', 'East');
            set(h,'XTick',[0,1])
            xlabs = [0,max(Cmap)];
            set(h,'XTickLabel',xlabs)
            set(h, 'Position', [0.6822    0.0842    0.0287    0.8717]);
            set(h, 'AxisLocation' ,'in');
         end
   
end
