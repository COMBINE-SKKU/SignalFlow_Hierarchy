function map = black2cyan(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...

'0a'; '00'; '01'; 
% '01';'62';'9c';
'00';'a2';'c5';
'00';'a2';'c5';
% '00'; 'cc'; 'ff'; % light blue
'00';'ff';'f9';
'00'; 'ff'; 'ff';
% '00'; 'ff'; 'ff';
 
 %version 2
% '00';'00';'00';
% '00';'19';'19';
% '00';'33';'33';
% '00';'4c';'4c';
% '00';'66';'66';
% '00';'7f';'7f';
% '00';'99';'99';
% '00';'b2';'b2';
% '00';'cc';'cc';
% '00';'e5';'e5';
% '00';'ff';'ff';



]; % black

values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');