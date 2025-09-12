function [plot_data_with_borders,parcel_IDs,parcel_names] = plot_HCP_boundaries_flat(plot_data_without_borders,which_areas,clr,atlas)

addpath(genpath('/local_raid1/03_user/younghyun/02_data/parcellations'));
if strcmp(which_areas,'SSM_yeo_som')
    if ~exist('borders_HCP_w_som.mat','file')
        error('Create file with borders of HCP somatotopic sub-areas using get_border_to_mat.m')
    else
        load('borders_HCP_w_som','borders_left','borders_right'); % borders of regions 1,2,3a,3b,5 replaced by hand, face, eye, foot, trunk
    end
    which_areas = 'SSM_yeo';
else
%     if ~exist('borders_HCP.mat','file')
%         error('Create file with borders of HCP parcellation using get_border_to_mat.m')
%     else
        switch atlas
            case "MMP"
                load('borders_HCP_flat.mat','borders_left','borders_right');
            case "Schaefer"
                load('borders_Schaefer100.mat','borders_left','borders_right');
            case "Markov"
                load('borders_Markov_32k_flat.mat','borders_left','borders_right');
        end
%     end
end

if strcmp(which_areas,'all')
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
    
    parcel_IDs = [];
    parcel_names = [];
else
    
    if strcmp(which_areas,'somatotopy')        
        if ~exist('borders_somatotopy.mat','file')
            error('Create file with borders of HCP omatotopic sub-areas using get_border_to_mat.m')
        else
            som = load('borders_somatotopy.mat','borders_left','borders_right');
        end
    end
    
    % cortical areas, grouped according to HCP parcellation paper (names and
    % HCP atlas labels)
    % neuroanatomical SI
    % visual
    early_vis = {'V1','V2','V3','V4'};
    early_vis_IDs = [1,4,5,6];
    dorsal_stream = {'V3A','V7','V3B','V6','V6A','IPS1'}; % "where" 
    dorsal_stream_IDs = [13,16,19,3,152,17];
    ventral_stream = {'V8','VVC','VMV1','VMV2','VMV3','PIT','FFC'}; % "what"
    ventral_stream_IDs = [7,163,153,160,154,22,18];
    MT_plus = {'V3CD','LO1','LO2','LO3','MT','MST','V4t','FST','PH'};
    MT_plus_IDs = [158,20,21,23,159,2,156,157,138];
    % put them together
    all_vis_IDs = [early_vis_IDs,dorsal_stream_IDs,ventral_stream_IDs,MT_plus_IDs];
    
    % somatosensory and motor
    % abbreviations:
    % PLMCC = paracentral lobular and mid cingulate cortex
    % POC = posterior opercular cortex
    SSM_core = {'4','3a','3b','1','2'};
    SSM_core_IDs = [8,53,9,51,52];
    PLMCC = {'5L','5m','5mv','24dd','24dv','6mp','6ma','SCEF'};
    PLMCC_IDs = [39,36,37,40,41,55,44,43];
    premot = {'6a','6d','FEF','PEF','55b','6v','6r'};
    premot_IDs = [96,54,10,11,12,56,78];
    POC = {'43','FOP1','OP4','OP2-3','OP1','PFcm'};
    POC_IDs = [99,113,100,102,101,105];
    % put them together
    all_SSM_IDs = [SSM_core_IDs,PLMCC_IDs,premot_IDs,POC_IDs];
    % higher hand area 24dd
    area_24dd = {'24dd'};
    area_24dd_ID = 40;
    
    % POS2 - not known what function, but has strong FC with RSC (both in PCC)
    area_POS2 = {'POS2'};
    area_POS2_ID = 15;   
    
    % somatotopic areas
    % these are separate from the parcellation - label IDs don't correspond
    somatotopy = {'foot','trunk','hand','eye','face'};
    somatotopy_IDs = [1,2,3,4,5];
    
    % auditory
    early_aud = {'A1','MBelt','LBelt','PBelt','RI'};
    early_aud_IDs = [24,173,174,124,104];
    ass_aud = {'A4','A5','STSdp','STSda','STSvp','STSva','TA2','STGa'};
    ass_aud_IDs = [175,125,129,128,130,176,107,123];
    % put them together
    all_AUD_IDs = [early_aud_IDs,ass_aud_IDs];
    
    % MULTIMODAL AREAS
    % audiovisual integration
    % tpo=temporo-parieto-occipital
    tpo_junction = {'TPOJ1','TPOJ2','TPOJ3','STV','PSL'};
    tpo_junction_IDs = [139,140,141,28,25];
    
    % sensory-motor pathway (visual to motor/somatosensory)
    % SPC = superior parietal cortex
    SPC = {'LIPv','LIPd','VIP','AIP','MIP','7PC','7AL','7Am','7PL','7Pm'};
    SPC_IDs = [48,95,49,117,50,47,42,45,46,29];
    
    % add borders by setting corresponding vertices to a value slightly above
    % the max value
    if strcmp(which_areas,'VIS')
        parcel_IDs = all_vis_IDs;
        parcel_names = [early_vis,dorsal_stream,ventral_stream];
    elseif strcmp(which_areas,'earlyVIS')
        parcel_IDs = early_vis_IDs;
        parcel_names = early_vis;
    elseif strcmp(which_areas,'dorsal_stream')
        parcel_IDs = dorsal_stream_IDs;
        parcel_names = dorsal_stream;
    elseif strcmp(which_areas,'ventral_stream')
        parcel_IDs = ventral_stream_IDs;
        parcel_names = ventral_stream;
    elseif strcmp(which_areas,'MT_plus')
        parcel_IDs = MT_plus_IDs;
        parcel_names = MT_plus;
    elseif strcmp(which_areas,'SSM')
        parcel_IDs = all_SSM_IDs;
        parcel_names = [SSM_core,PLMCC,premot,POC];
    elseif strcmp(which_areas,'premot')
        parcel_IDs = premot_IDs;
        parcel_names = premot;
    elseif strcmp(which_areas,'somatotopy')
        parcel_IDs = somatotopy_IDs;
        parcel_names = somatotopy;
    elseif strcmp(which_areas,'AUD')
        parcel_IDs = all_AUD_IDs;
        parcel_names = [early_aud,ass_aud];
    elseif strcmp(which_areas,'earlyAUD')
        parcel_IDs = early_aud_IDs;
        parcel_names = early_aud;
    elseif strcmp(which_areas,'TPO')
        parcel_IDs = tpo_junction_IDs;
        parcel_names = tpo_junction;
    elseif strcmp(which_areas,'SPC')
        parcel_IDs = SPC_IDs;
        parcel_names = SPC;
    elseif strcmp(which_areas,'24dd')
        parcel_IDs = area_24dd_ID;
        parcel_names = area_24dd;
    elseif strcmp(which_areas,'POS2')
        parcel_IDs = area_POS2_ID;
        parcel_names = area_POS2;
    elseif strcmp(which_areas,'LIPv_VIP_complex')
        parcel_IDs = [48,49];
        parcel_names = {'LIPv','VIP'};
    elseif strcmp(which_areas,'55b')
        parcel_IDs = 12;
        parcel_names = '55b';
    elseif strcmp(which_areas,'SCEF')
        parcel_IDs = 43;
        parcel_names = 'SCEF';
    elseif strcmp(which_areas,'PCC')
        parcel_names = {'DVT','ProS','POS1','POS2','RSC','v23ab','d23ab','31pv','31pd','31a','23d','23c','PCV','7m'}; % POS: transition early vis - PCC
        parcel_IDs = [142,121,31,15,14,33,34,35,161,162,32,38,27,30];
    elseif strcmp(which_areas,'dlPFC')
        parcel_names = {'8C','8Av','i6-8','s6-8','SFL','8BL','9p','9a','8Ad','p9-46v','a9-46v','46','9-46d'};
        parcel_IDs = [73,67,97,98,26,70,71,87,68,83,85,84,86];
    elseif strcmp(which_areas,'orbital_polar_FC')
        parcel_names = {'47s','47m','a47r','11l','13l','a10p','p10p','10pp','10d','OFC','pOFC'};
        parcel_IDs = [94,66,77,91,92,89,170,90,72,93,166];
    elseif strcmp(which_areas,'insula_frontal_opercular')
        parcel_names = {'52','PI','Ig','PoI1','PoI2','FOP2','FOP3','MI','AVI','AAIC','Pir','FOP4','FOP5'};
        parcel_IDs = [103,178,168,167,106,115,114,109,111,112,110,108,169];
        
        % yeo RSNs
    elseif ~isempty(which_areas) && strcmp(which_areas(end-2:end),'yeo')
        load('parcels_for_plotting_som','parcel_IDs','parcel_names')
        if strcmp(which_areas,'VIS_yeo')
            yp = 1;
        elseif strcmp(which_areas,'SSM_yeo')
            yp = 2;
        elseif strcmp(which_areas,'DAtt_yeo')
            yp = 3;
        elseif strcmp(which_areas,'VAtt_yeo')
            yp = 4;
        elseif strcmp(which_areas,'limbic_yeo')
            yp = 5;
        elseif strcmp(which_areas,'FPN_yeo')
            yp = 6;
        elseif strcmp(which_areas,'DMN_yeo')
            yp = 7;
        end
        parcel_IDs = parcel_IDs{1}{yp}(parcel_IDs{1}{yp}<181);
        parcel_names = parcel_names{1}{yp};
        
    % nothing 
    elseif strcmp(which_areas,'')
        parcel_IDs = [];
        parcel_names = {};
    else
        error('Areas not recognized.');
    end
    
    plot_borders_left = any(borders_left(:,parcel_IDs),2);
    plot_borders_right = any(borders_right(:,parcel_IDs),2);
    plot_all_borders = [plot_borders_left;plot_borders_right];
    
    plot_data_with_borders = plot_data_without_borders;
    % automatically map colors
    if size(plot_data_without_borders,2) == 1
        plot_data_with_borders(plot_all_borders) = max(plot_data_without_borders)+(max(plot_data_without_borders)-min(plot_data_without_borders))/65;
    % explicitly assign colors 
    elseif size(plot_data_without_borders,2) == 3
        if isstr(clr)
            clr = rgb(clr);
        end
        plot_data_with_borders(plot_all_borders,:) = repmat(clr,[sum(plot_all_borders),1]);
    end
end
end