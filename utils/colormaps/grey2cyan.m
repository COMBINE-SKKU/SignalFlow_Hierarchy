function map = grey2cyan(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'80';'80';'80';
'80';'80';'80';
'80';'80';'80';
'00';'63';'79';
'00';'78';'92';
'00';'8d';'ac';
'00';'a2';'c5';
'00';'a2';'c5';
'66';'c7';'dc';
'66';'c7';'dc';
'00';'ff';'f9';
'00';'ff';'ff';
]; 

values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');