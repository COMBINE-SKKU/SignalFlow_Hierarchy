function map = Campbell(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

clrData = '/combinelab/03_user/younghyun/04_software/colormaps/campbellColormap.mat';
values = importdata(clrData);


P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');