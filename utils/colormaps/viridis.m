function map = viridis(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
% '80';'80';'80';
% '80';'80';'80';
'24';'01';'54';
'2c';'01';'54';
'44';'01';'54';
'44';'39';'83';
'44';'39';'83';
'31';'68';'8e';
'21';'91';'8c';
'21';'91';'8c';
'35';'b7';'79';
'90';'d7';'43';
'fd';'e7';'25'
'ff';'ff';'00';
'ff';'ff';'00';
];  % red


values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');