function GRASP = run_grasp(BOLD, threshold, iter, varargin)
% RUN_GRASP Performs subsampling of the BOLD signal using GRASP
% and computes the average of the subsampled signals.
%
% BOLD      : BOLD signal time series
% threshold : Threshold for the GRASP algorithm
% iter      : Number of iterations for the GRASP algorithm
% Optional parameters:
%   num_start : Number of random restarts (integer, default=1)

% Parse optional arguments
num_start = 1;
if ~isempty(varargin)
    num_start = varargin{1};
end

% Initialize a cell array to store the GRASP results
N = size(BOLD,2);
GRASPs = zeros(N,N);

% Perform subsampling using GRASP in parallel   
parfor j = 1:iter
    temp = subsampleGRASP(BOLD, j, num_start);
    GRASPs = GRASPs+temp;
end

% Compute the mean of the concatenated results
GRASP = GRASPs ./ iter;

% threshold
GRASP(GRASP<threshold) = 0;

% Delete temporary files
delete BOLD* GRASP* causal-cmd.log*
end

function temp = subsampleGRASP(BOLD, j, num_start)
% SUBSAMPLEGRASP Subsamples the BOLD signal and applies the GRASP algorithm.
%
% BOLD      : BOLD signal time series
% j         : Index for differentiating file names
% num_start : Number of random restarts

% Select subsamples
nobs = size(BOLD, 1);
BOLDsignal = datasample(BOLD, ceil(nobs/2), 1, 'Replace', false);

% Write the BOLD signal to a text file
writematrix(BOLDsignal, sprintf('BOLDsignal%d', j))

% Run GRASP algorithm using the command line
command = sprintf('java  -jar causal-cmd-1.12.0-jar-with-dependencies.jar --algorithm grasp --test sem-bic-test --score sem-bic-score --data-type continuous --dataset BOLDsignal%d.txt --delimiter comma --no-header --prefix GRASP%d --numStart %d', j, j, num_start);
system(command)

% Convert the output to a connectivity matrix
temp = Tetrad2Matrix(sprintf('GRASP%d_out.txt', j), 'directed','grasp',0)'; 

end