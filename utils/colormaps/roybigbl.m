function map = roybigbl(n)

 

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'00'; 'ff'; 'ff';  % cyan
'15'; 'b2'; '00';  % green
'72'; '00'; '98';  % purple
'00'; '00'; '97';  % dark-blue
'0a'; '00'; '00';  % black
'6b'; '00'; '00';  % brown
'ca'; '00'; '00';  % dark-red
'ff'; '77'; '00';  % orange
'f8'; 'f6'; '00'];  % yellow


values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');