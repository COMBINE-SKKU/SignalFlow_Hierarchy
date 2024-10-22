function map = white2red(n)


if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

values = [...
'ff'; 'fb'; 'ee';  
'ff'; 'f4'; 'd7';  
'ff'; 'f4'; 'd7';  
'ff'; 'ec'; 'c4';  
'ff'; 'ec'; 'c4';  
'fe'; 'd7'; '9a'; 
'fe'; 'd7'; '9a'; 
'fe'; 'c9'; '82';  
'fe'; 'c9'; '82';  
'ff'; 'b7'; '69';  
'ff'; 'b7'; '69';  
'ff';'9f';'45';
'ff';'9f';'45';
'ff'; '8f'; '2b'; 
'ff'; '8f'; '2b'; 
'ff'; '77'; '0f';  
'ff'; '7a'; '0d'; 

]; 


values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');
