function map = MMP22moduleColors


% Additional hex values from switch case
values = [...
'80'; '80'; '80';
'00'; '00'; '80';  % Corresponds to '#000080'  
'01'; '00'; 'FF';  % Corresponds to '#0100FF'
'80'; '80'; 'FF';  % Corresponds to '#8080FF'
'60'; '00'; 'C0';  % Corresponds to '#6000C0'
'26'; '80'; '80';  % Corresponds to '#268080'
'59'; 'FF'; '00';  % Corresponds to '#59FF00'
'C0'; 'FF'; '80';  % Corresponds to '#C0FF80'
'80'; 'FF'; '7F';  % Corresponds to '#80FF7F'
'80'; 'FF'; '03';  % Corresponds to '#80FF03'
'F5'; '0E'; '02';  % Corresponds to '#F50E02'
'80'; '03'; '00';  % Corresponds to '#800300'
'F8'; 'A0'; 'C0';  % Corresponds to '#F8A0C0'
'40'; '30'; '40';  % Corresponds to '#403040'
'20'; '20'; '20';  % Corresponds to '#202020'
'C0'; '60'; '20';  % Corresponds to '#C06020'
'57'; 'FF'; 'C0';  % Corresponds to '#57FFC0'
'80'; '80'; '80';  % Corresponds to '#808080'
'80'; '20'; '60';  % Corresponds to '#802060'
'60'; '40'; '20';  % Corresponds to '#604020'
'C0'; 'C0'; 'C0';  % Corresponds to '#C0C0C0'
'30'; '20'; '10';  % Corresponds to '#302010'
'E0'; 'C0'; 'C0'   % Corresponds to '#E0C0C0'
];



values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

P = size(values,1);

map = interp1(1:size(values,1), values, linspace(1,P,23), 'linear');