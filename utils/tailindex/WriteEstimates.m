%% WriteEstimates
%
% First version: Richard Tol, 8 November 2011
% This version: Richard Tol, 8 November 2011

display('Write estimates');

vFileToOpen = 'C:\Users\rtol\Google Drive\TailIndex\tailindices.xlsx';
%vEst = [HillEst' DekkersEst' AlvesEst' CG15' CG20' CG25' LPN15' LPN20' LPN25'];
vEst = [HillEst'];
xlswrite(vFileToOpen, vEst, 'TailIndex', 'B2');

clear v*