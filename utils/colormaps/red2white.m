function map = red2white(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'bc';'02';'02';
'c2';'1b';'1b';
% 'c9';'34';'34';
'd0';'4d';'4d';
% 'd6';'67';'67';
'dd';'80';'80';
% 'e4';'99';'99';
'ea';'b3';'b3';
% '80';'80';'80';
% 'f1';'cc';'cc';
'f8';'e5';'e5';
'ff';'ff';'ff';

]; % black

values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');