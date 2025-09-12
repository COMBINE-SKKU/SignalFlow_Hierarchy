function map = grey2white(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'80';'80';'80';
'80';'80';'80';
'80';'80';'80';
'e9';'3e';'3a';
'e9';'3e';'3a';
'ed';'68';'3c';
'ed';'68';'3c';
'f3';'90';'3f';
'f3';'90';'3f';
'fd'; 'c7'; '0c';
'fd'; 'c7'; '0c';
'ff';'ff';'00']; 

values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');