function map = blue2red(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'4c'; '4c'; '4c';  % grey
'1c'; '1d'; 'ff'; % Blue
'59'; '59'; 'ff';  % light-blue
'b8'; 'b8'; 'ff';  
'eb'; 'eb'; 'ff';  % white-ish
'fb'; 'cf'; 'cf';  % light-red
'f8'; '96'; '98';  %
'f6'; '6a'; '6a';  %
'f5'; '4e'; '4e';  % orange
'f5'; '0e'; '02'];  % red


values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');