function map = EKcolormap(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'80';'80';'80';
'7F';'29';'7F';
'35';'68';'9B';
'A8';'D3';'8E';
'FD';'CD';'0B';
'F3';'EC';'1A';


];  % red


values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');