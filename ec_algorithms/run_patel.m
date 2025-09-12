function Patel = run_patel(BOLD, threshold, iter)
% RUN_PATEL Performs subsampling of the BOLD signal using Patel
% and computes the average of the subsampled signals.
%
% BOLD   : BOLD signal time series
% threshold : Threshold for the Patel algorithm
% iter : Number of iterations for the Patel algorithm

% Initialize a cell array to store the Patel results
N = size(BOLD,2);
Patels = zeros(N,N);

% Perform subsampling using Patel in parallel   
parfor j = 1:iter
    temp= subsamplePatel(BOLD, j);
    Patels = Patels+temp;
end

% Compute the mean of the concatenated results
Patel = Patels ./ iter;

% threshold
Patel(Patel<threshold) = 0;

% Delete temporary files
delete BOLDsignal* FAS* causal-cmd.log*
end

function temp = subsamplePatel(BOLD, j)
% SUBSAMPLEPATEL Subsamples the BOLD signal and applies the Patel algorithm.
%
% BOLD   : BOLD signal time series
% j      : Index for differentiating file names

% Select subsamples
nobs = size(BOLD, 1);
BOLDsignal = datasample(BOLD, ceil(nobs/2), 1, 'Replace', false);

% Write the BOLD signal to a text file
writematrix(BOLDsignal, sprintf('BOLDsignal%d', j))

% Run FAS-stable algorithm
command = sprintf('java -jar causal-cmd-1.12.0-jar-with-dependencies.jar --algorithm fas --stableFAS --test sem-bic-test --data-type continuous --dataset BOLDsignal%d.txt --delimiter comma --no-header --prefix FAS%d', j, j);
system(command)

% Convert the output to a connectivity matrix
fas = Tetrad2Matrix(sprintf('FAS%d_out.txt', j), 'undirected','fas',0);

% Run Patel algorithm
temp = patel_tau(BOLDsignal,fas);

end

function [mat] = patel_tau(data,FAS)
%--------------------------------------------------------------------------
% this function implements patel's tau method described in Ramsey et al.,
% 2014 paper [Non-Gaussian methods and high-pass filters in the 
% estimation of effective connections]
%It takes for input timeseries data and adjacency matrix
%--------------------------------------------------------------------------

BOLDsignal = data;
n = size(BOLDsignal,2);
mat = zeros(n,n);
temp = tril(FAS);

[mi,mj] = find(temp>0); 

% calculate percentile and map data into interval [0,1]
for i = 1:size(mi,1)
    
    %convert first node in adjacency matrixcl
    temp = BOLDsignal(:,mi(i));
    Xvals = prctile(temp,[10 90]);
    a1 = (temp-Xvals(1,1))/(Xvals(1,2)-Xvals(1,1));
    a2 = min(a1,1);
    X2 = max(a2,0);
    
    %conver second node in adjacency matrix
    temp = BOLDsignal(:,mj(i));
    Yvals = prctile(temp,[10 90]);
    b1 = (temp-Yvals(1,1))/(Yvals(1,2)-Yvals(1,1));
    b2 = min(b1,1);
    Y2 = max(b2,0);
    
    
    theta1 = dot(X2,Y2);
    theta2 = dot(X2,(1-Y2));
    theta3 = dot((1-X2),Y2);
    if theta2 > theta3
        tau = (1-(theta1+theta3))/(theta1+theta2);
    else
        tau = (theta1+theta2)/((theta1+theta3)-1);
    end
    
    if tau < 0
        mat(mi(i),mj(i)) = 1;
    else
        mat(mj(i),mi(i)) = 1;
    end
                
end
end
    