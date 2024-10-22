function map = cyan2black(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'00'; 'ff'; 'ff';
'00';'ff';'f9';
'00';'a2';'c5';
'03';'6c';'b0';
'00';'5b';'96';
'0a'; '00'; '01';  % black
]; 

values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');


