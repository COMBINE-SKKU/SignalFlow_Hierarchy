function FASK = subsample_wrapper(BOLD, oAlpha, alpha, threshold, iter)
% SUBSAMPLE_WRAPPER Performs subsampling of the BOLD signal using FASK and computes the average of the subsampled signals.
%
% BOLD   : BOLD signal time series
% i      : Index for differentiating file names
% oAlpha : Orientation alpha value for the FASK algorithm
% alpha  : Significance level for the FASK algorithm

% Initialize a cell array to store the FASK results
N = size(BOLD,2);
FASKs = zeros(N,N);

% Perform subsampling using FASK in parallel   
parfor j = 1:iter
    temp= subsampleFASK(BOLD, j, oAlpha, alpha);
    FASKs = FASKs+temp;
end

% Compute the mean of the concatenated results
FASK = FASKs ./ iter;

% threshold
FASK(FASK<threshold) = 0;

% Delete temporary files
delete BOLD* FASK* FAS*
end

function temp = subsampleFASK(BOLD, j, oAlpha, alpha)
% SUBSAMPLEFASK Subsamples the BOLD signal and applies the FASK algorithm.
%
% BOLD   : BOLD signal time series
% j      : Index for differentiating file names
% i      : Index for differentiating file names
% oAlpha : Orientation alpha value for the FASK algorithm
% alpha  : Significance level for the FASK algorithm

% Select subsamples
nobs = size(BOLD, 1);
BOLDsignal = datasample(BOLD, ceil(nobs/2), 1, 'Replace', false);

% Write the BOLD signal to a text file
writematrix(BOLDsignal, sprintf('BOLDsignal%d', j))

% Run FASK algorithm using the command line
command = sprintf('java  -jar causal-cmd-1.4.1-jar-with-dependencies.jar --algorithm fask --test fisher-z-test --orientationAlpha %d --faskAdjacencyMethod 1 --alpha %d --data-type continuous --dataset BOLDsignal%d.txt --delimiter comma --no-header --prefix FASK%d --score sem-bic-score', oAlpha, alpha, j, j);
system(command)

% Convert the output to a connectivity matrix
temp = Tetrad2Matrix(sprintf('FASK%d.txt', j), 'directed')';

end