function map = pastel1(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'80';'80';'80';
'2c';'5d';'37';
'e3';'c5';'15';
'ee';'51';'b1';
'a5';'9c';'d3';
'4b';'2d';'9f';
% '83';'3e';'5b';
% 'f7';'c9';'94';
];  % red


values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');