function BOSS = run_boss(BOLD, threshold, iter, varargin)
% SUBSAMPLE_BOSS Performs subsampling of the BOLD signal using BOSS and computes the average of the subsampled signals.
%
% BOLD      : BOLD signal time series
% threshold : Threshold for the BOSS algorithm
% iter      : Number of iterations for subsampling
% Optional parameters:
%   num_start : Number of random restarts (integer, default=1)
%   bes       : Use BES score (0 or 1, default=0)

% Parse optional arguments
num_start = 1;
bes = 0;
if ~isempty(varargin)
    num_start = varargin{1};
    if length(varargin) > 1
        bes = varargin{2};
    end
end

num_start = num_start;
bes = bes;

% Initialize a cell array to store the BOSS results
N = size(BOLD,2);
BOSSs = zeros(N,N);

% Perform subsampling using BOSS in parallel   
parfor j = 1:iter
    temp = subsampleBOSS(BOLD, j, num_start, bes);
    BOSSs = BOSSs+temp;
end

% Compute the mean of the concatenated results
BOSS = BOSSs ./ iter;

% threshold
BOSS(BOSS<threshold) = 0;

% Delete temporary files
delete BOLD* BOSS* causal-cmd.log*
end

function temp = subsampleBOSS(BOLD, j, num_start, bes)
% SUBSAMPLEBOSS Subsamples the BOLD signal and applies the BOSS algorithm.
%
% BOLD      : BOLD signal time series
% j         : Index for differentiating file names
% num_start : Number of random restarts
% bes       : Use BES score flag

% Select subsamples
nobs = size(BOLD, 1);
BOLDsignal = datasample(BOLD, ceil(nobs/2), 1, 'Replace', false);

% Write the BOLD signal to a text file
writematrix(BOLDsignal, sprintf('BOLDsignal%d', j))

% Run BOSS algorithm using the command line
command = sprintf('java -jar causal-cmd-1.12.0-jar-with-dependencies.jar --algorithm boss --numThreads 30 --numStarts %d --useBes %d --data-type continuous --dataset BOLDsignal%d.txt --delimiter comma --no-header --prefix BOSS%d --score sem-bic-score', num_start, bes, j, j);
system(command)

% Convert the output to a connectivity matrix
temp = Tetrad2Matrix(sprintf('BOSS%d_out.txt', j), 'directed','boss',0)';
end