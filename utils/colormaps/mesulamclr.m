function map = mesulamclr(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'80';'80';'80';
'96';'AD';'DC';
'F8';'F5';'4D';
'F3';'B5';'A3';
'AA';'C6';'56';
];  % red


values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');