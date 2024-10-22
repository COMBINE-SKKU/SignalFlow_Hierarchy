function map = pc1color

% % Define the number of points you want to sample
% first = round(256 * 0.9);
% second = 256 - first;
% 
% % Create a linear space to sample from the colormaps
% % MATLAB's indexing starts at 1, hence the adjustment in linspace
% x_viridis = linspace(1, first, first);
% x_YlOrBr = linspace(1, 256, second);
% 
% % Interpolate to sample colors from the colormaps
% % 'interp1' is used for 1-D linear interpolation
% colors2 = interp1(1:256, slanCM('viridis'), x_viridis);
% colors3 = interp1(1:256, slanCM('YlOrBr'), x_YlOrBr);


startFraction_viridis = 0;
endFraction_viridis = 1;
startFraction_YlOrBr = 0.25;
endFraction_YlOrBr = 1.0;

% Calculate the start and end indices for truncation
startIndex_viridis = round(256 * startFraction_viridis) + 1;
endIndex_viridis = round(256 * endFraction_viridis);
startIndex_YlOrBr = round(256 * startFraction_YlOrBr) + 1;
endIndex_YlOrBr = round(256 * endFraction_YlOrBr);

clr1 = slanCM('viridis');
clr2 = slanCM('YlOrBr');

% Truncate the colormaps
colors2 = clr1(startIndex_viridis:endIndex_viridis, :);
colors3 = clr2(startIndex_YlOrBr:endIndex_YlOrBr, :);
% Combine the sampled colors to create a new colormap
map = [colors2; colors3];
