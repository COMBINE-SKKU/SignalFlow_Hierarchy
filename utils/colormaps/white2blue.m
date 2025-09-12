function map = white2blue(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'ff';'ff';'ff';
'ff';'ff';'ff';
'ff';'ff';'ff';
% '80';'80';'80';
'e8';'e8';'fa';
'd2';'d2';'f5';
'bc';'bc';'f0';
'a6';'a6';'eb';
'90';'90';'e7';
'79';'79';'e2';
'63';'63';'dd';
'4d';'4d';'d8';
'37';'37';'d3';
'21';'21';'cf';


]; % black

values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');