function CCD = run_ccd(BOLD, threshold, iter)
% RUN_CCD Performs subsampling of the BOLD signal using CCD
% and computes the average of the subsampled signals.
%
% BOLD   : BOLD signal time series
% threshold : Threshold for the CCD algorithm
% iter : Number of iterations for the CCD algorithm

% Initialize a cell array to store the CCD results
N = size(BOLD,2);
CCDs = zeros(N,N);

% Perform subsampling using CCD in parallel   
parfor j = 1:iter
    temp= subsampleCCD(BOLD, j);
    CCDs = CCDs+temp;
end

% Compute the mean of the concatenated results
CCD = CCDs ./ iter;

% threshold
CCD(CCD<threshold) = 0;

% Delete temporary files
delete BOLD* CCD* causal-cmd.log*
end

function temp = subsampleCCD(BOLD, j)
% SUBSAMPLECCD Subsamples the BOLD signal and applies the CCD algorithm.
%
% BOLD   : BOLD signal time series
% j      : Index for differentiating file names

% Select subsamples
nobs = size(BOLD, 1);
BOLDsignal = datasample(BOLD, ceil(nobs/2), 1, 'Replace', false);

% Write the BOLD signal to a text file
writematrix(BOLDsignal, sprintf('BOLDsignal%d', j))

% Run CCD algorithm using the command line
command = sprintf('java  -jar causal-cmd-1.12.0-jar-with-dependencies.jar --algorithm ccd --test sem-bic-test --data-type continuous --dataset BOLDsignal%d.txt --delimiter comma --no-header --prefix CCD%d', j, j);
system(command)

% Convert the output to a connectivity matrix
temp = Tetrad2Matrix(sprintf('CCD%d_out.txt', j), 'pattern','ccd',0)';

end