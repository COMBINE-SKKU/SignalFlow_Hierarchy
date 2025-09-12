function xs = solve_hopf_sde(omega, a, wC, dt, Tmax, TR, sig)
% SOLVE_HOPF_SDE Simulates time series based on the Hopf model using a
% stochastic differential equation (SDE) solver.
%
% omega : Oscillatory frequency (in radians) for each node (calculated beforehand)
% a     : Hopf bifurcation parameter that transitions between oscillatory and noisy signals
% wC    : Weighted coupling matrix
% dt    : Time step
% Tmax  : Max number of volumes in the fMRI time course
% TR    : Repetition Time for each volume
% sig   : Strength of the noise terms

% Determine the number of nodes
N = size(wC, 1);

% Set up the degree matrix
sumC = repmat(sum(wC, 2), 1, 2); % for sum Cij*xj

% Initialize the time series
xs = zeros(Tmax, N);

% Initialize the z variable, which captures both x and y variables
z = 0.1 * ones(N, 2);

% Calculate the noise term step for the Euler-Maruyama method
dsig = sqrt(dt*sig);

% Perform initial simulation to discard the first 3000 time steps
for t = 0:dt:1000
    suma = wC * z - sumC .* z;
    z = euler_maruyama_update(z, dt, suma, a, omega, dsig, N);
end

% Initialize variables for actual modeling
timeVector = 0:dt:((Tmax - 1) * TR);
downsampleRate = floor(TR / dt);

% Perform the actual modeling (x = BOLD signal, y = some other oscillation)
nn = 1; % Initialize the index for storing downsampled data
for tind = 1:length(timeVector)
    suma = wC * z - sumC .* z;
    z = euler_maruyama_update(z, dt, suma, a, omega, dsig, N);
    
    % Downsample to the same TR
    if mod(tind-1, downsampleRate) == 0 && nn <= Tmax
        xs(nn, :) = z(:, 1)';
        nn = nn + 1;
    end
end
end

function z = euler_maruyama_update(z, dt, suma, a, omega, dsig, N)
    zz = z(:,end:-1:1);

    % Calculate the drift term
    drift = a .* z + zz .* omega - z .* (z .* z + zz .* zz) + suma;
    
    % Calculate the diffusion term
    diffusion = dsig * randn(N, 2);
    
    % Euler-Maruyama method update
    z = z + dt * drift + diffusion;
end