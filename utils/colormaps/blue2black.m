function map = blue2black(n)

if nargin < 1
   n = 8; % Defining the default number of levels to be 8
end

% Defining the shades ranging from black to blue to cyan
values = [...
'00'; '00'; '00';  % Black
'00'; '00'; '40';  % Dark blue
'00'; '00'; '80';  % Blue
'00'; '20'; 'bf';  % Blue-cyan
'00'; '40'; 'df';  % Brighter blue-cyan
'00'; '60'; 'ef';  % Cyan-blue
'10'; '80'; 'f7';  % Light cyan-blue
'30'; 'a0'; 'ff';  % Cyanish-blue
'00'; 'c0'; 'ff';  % Cyan
'00'; 'ff'; 'ff';  % Bright cyan
];

values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');
