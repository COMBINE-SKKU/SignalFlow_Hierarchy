function map = vidden(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'00'; '00'; '00';  % black
'5c'; '0a'; '38';  % dark-purple
'3f'; '3f'; '64';  % navy
'74'; '74'; 'bb';  % light-purple
'01'; 'fc'; '01';  % green
'4d'; 'c4'; '0c';  % dark-green
'ff'; 'ce'; '00';  % gold
'ff'; '76'; '00';  % orange
'fb'; '03'; '00'];  % red


values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');