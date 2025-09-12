function LINGAM = run_lingam(BOLD, threshold, iter)
% RUN_LINGAM Performs subsampling of the BOLD signal using LiNGAM
% and computes the average of the subsampled signals.
%
% BOLD   : BOLD signal time series
% threshold : Threshold for the LiNGAM algorithm
% iter : Number of iterations for the LiNGAM algorithm

% Initialize a cell array to store the LiNGAM results
N = size(BOLD,2);
LINGAMs = zeros(N,N);

% Perform subsampling using LiNGAM in parallel   
parfor j = 1:iter
    temp= subsampleLiNGAM(BOLD, j);
    LINGAMs = LINGAMs+temp;
end

% Compute the mean of the concatenated results
LINGAM = LINGAMs ./ iter;

% threshold
LINGAM(LINGAM<threshold) = 0;

% Delete temporary files
delete BOLD* LINGAM* causal-cmd.log*
end

function temp = subsampleLiNGAM(BOLD, j)
% SUBSAMPLELINGAM Subsamples the BOLD signal and applies the LiNGAM algorithm.
%
% BOLD   : BOLD signal time series
% j      : Index for differentiating file names

% Select subsamples
nobs = size(BOLD, 1);
BOLDsignal = datasample(BOLD, ceil(nobs/2), 1, 'Replace', false);

% Write the BOLD signal to a text file
writematrix(BOLDsignal, sprintf('BOLDsignal%d', j))

% Run LiNGAM algorithm using the command line
command = sprintf('java  -jar causal-cmd-1.12.0-jar-with-dependencies.jar --algorithm direct-lingam --score sem-bic-score --data-type continuous --dataset BOLDsignal%d.txt --delimiter comma --no-header --prefix LINGAM%d', j, j);
system(command)

% Convert the output to a connectivity matrix
temp = Tetrad2Matrix(sprintf('LINGAM%d_out.txt', j), 'directed','lingam',0)';

end