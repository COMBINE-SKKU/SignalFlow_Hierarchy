function map = blue_grey(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'80'; '80'; '80';  % black
'21'; '21'; 'cf';
];  % yellow

values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');
