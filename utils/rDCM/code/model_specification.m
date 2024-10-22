%model specification
function [DCM] = model_specification(dms,a,b,c,d)

DCM.Y.y = dms;
DCM.Y.dt = 1;
for Nr_region = 1:size(DCM.Y.y,2)
    DCM.Y.name{Nr_region} = ['region_' num2str(Nr_region)];
end
DCM.a = a;
DCM.b = b;
DCM.c = c;
DCM.d = d;
DCM.U.u = zeros(size(DCM.Y.y,1)*16,size(DCM.Y.y,2));
DCM.U.dt = 16;
for Nr_input = 1:size(DCM.U.u,2)
    DCM.U.name{Nr_input} = ['input_' num2str(Nr_input)];
end
% specify number of datapoints (per regions)
DCM.v = size(DCM.Y.y,1);

% specify number of regions
DCM.n = size(DCM.Y.y,2);
