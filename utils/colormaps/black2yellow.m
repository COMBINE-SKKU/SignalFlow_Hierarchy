function map = black2yellow(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'0a'; '00'; '01';  % black
'f5'; '0e'; '02';  
'f6'; '55'; '03';  % 
'f9'; '9e'; '03';  %
'fc'; 'd4'; '01';  %
'ff'; 'ff'; '01';  %
% 'ff'; 'ff'; '4e';  % orange
];  % yellow

values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');

