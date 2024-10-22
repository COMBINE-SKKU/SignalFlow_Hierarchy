function map = red2black(n)

% Defining the shades ranging from black to red to orange to yellow
values = [...
'00'; '00'; '00';  % Black
'40'; '00'; '00';  % Dark red
'80'; '00'; '00';  % Red
'bf'; '20'; '00';  % Red-orange
'df'; '40'; '00';  % Brighter red-orange
'ef'; '60'; '00';  % Orange
'ff'; '80'; '00';  % Light orange
'ff'; 'a0'; '20';  % Yellowish-orange
'ff'; 'c0'; '00'; % Yellow
'ff'; 'ff'; '00'; % Bright Yellow
];


values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');