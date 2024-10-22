function map = magma(n)

if nargin < 1
   n = size(get(gcf, 'Colormap'), 1);
end

% Magma colormap values
values = [...
'00';'00';'04'; % Dark indigo
'1c';'03';'29'; % Deep purple
'34';'02';'3f'; % Dark purple
'4e';'02';'4f'; % Purple
'6e';'0b';'4c'; % Medium purple
'8c';'20';'4b'; % Red-purple
'aa';'3b';'3f'; % Dark red
'c2';'55';'30'; % Reddish brown
'd8';'6f';'24'; % Orange-red
'e9';'89';'1d'; % Orange
'f5';'a7';'2a'; % Light orange
'fc';'cc';'5c'; % Yellowish
'ff';'f7';'b7'  % Light yellow
];  % red

values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');
