function signal = linear_propagation(EC, z)
% linear_propagation - Simulates the energy flow in a dynamic system over time.
%
% This function computes the energy flow of a system using an exponential
% matrix model (exponential of the matrix EC) to simulate how the system
% evolves over time when starting from an initial state 'z'.
%
% INPUTS:
%   EC - (NxN matrix) Effective connectivity matrix representing the system.
%   z  - (Nx1 vector) Initial state vector of the system.
%
% OUTPUT:
%   signal - (Tmax x N matrix) Time-series of the system's energy flow,
%            where each row represents the state of the system at a specific
%            time point, and each column represents the state of each element.

    % Define parameters
    N = size(EC, 1);       % Number of elements (size of the system)
    dt = 0.1;              % Time step
    ts = 0.1:dt:10;        % Time vector from 0.1 to 10 with a step of 3.1
    Tmax = length(ts);     % Total number of time points
    
    % Initialize the output signal matrix
    signal = zeros(Tmax, N);  % Pre-allocate signal matrix (Tmax x N)
    
    % Loop over each time point to calculate the energy flow
    for t = 1:Tmax
        % Scale the effective connectivity (EC) matrix by the current time step
        EC_s = EC * ts(t);
        
        % Compute the matrix exponential of the scaled EC matrix
        EC_e = expm(EC_s);
        
        % Compute the new state vector by applying EC_e to the initial state z
        z2 = EC_e * z;
        
        % Store the new state in the signal matrix
        signal(t, :) = z2';
    end

end