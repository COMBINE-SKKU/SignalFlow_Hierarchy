function map = cyan2white(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...

'00';'ff';'ff';
'00';'ff';'f9';
'00';'a2';'c5';
'01';'62';'9c';
'03';'39';'6c';
'01';'1f';'4b';
% '80';'80';'80';
% '80';'80';'80';
'55'; '01'; '02'; 
'bf'; '08'; '00';  % light-blue
'f5'; '0e'; '02';  
'f6'; '55'; '03';  % white-ish
'f9'; '9e'; '03';  % light-red
'fc'; 'd4'; '01';  %
'ff'; 'ff'; '01';  %
'ff'; 'ff'; '4e';  % orange
];  % red


values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');