function map = red2blue(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'd3';'42';'54';
'd7';'53';'61';
'e6';'81';'89';
'f1';'ab';'af';
'fc';'df';'e1';
'ff';'ff';'ff';
'e3';'f1';'fb';
'b8';'db';'ec';
'96';'c3';'e0';
'4e';'8e';'c0';
'46';'7b';'a3'];  % red


values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');