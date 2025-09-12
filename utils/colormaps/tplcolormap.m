function map = tplcolormap(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'80';'80';'80';
'58';'70';'aa';
'8b';'a7';'b0';
'ab';'ba';'a2';
'd9';'c5';'98';
'fd';'93';'77';
'fd';'cc';'fd';


];  % red


values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');