function BOLD = run_simulation(J, varargin)
% RUN_SIMULATION Simulates BOLD signals using the specified model
% BOLD = run_simulation(J, varargin)
% J: connectivity matrix
% varargin: 
% 1. G: global coupling parameter
% 2. f_diff: peak frequency vector
% 3. TR: repetition time
% 4. Tmax: total simulation time
% 5. scale_a: intrinsic coupling strength (default: -0.02)
% 6. sig: noise level (default: 0.01)


% Normalize the structural connectivity matrix (SC)
J = J-diag(diag(J)); % Remove self-connections (diagonal elements)

% Check input arguments
if isempty(varargin)
    assert(false, 'G, f_diff, TR, Tmax are required')
else
    G = varargin{1};
    f_diff = varargin{2};
    TR = varargin{3};
    Tmax = varargin{4};
    % Get optional parameters or use defaults
    if length(varargin) >= 5
        scale_a = varargin{5};
    else
        scale_a = -0.01;
    end
    if length(varargin) >= 6
        sig = varargin{6};
    else
        sig = 0.01;
    end
end

% Set simulation parameters
dt = 0.1*TR/2; % Time step

% Set parameters for Hopf model
N = size(J, 1);
wC = G * J; % Scale the connectivity matrix

% Define node frequencies and convert them to radians
omega = repmat(2 * pi * f_diff', 1, 2);
omega(:,1) = -omega(:,1);

% Set intrinsic coupling strength
a = scale_a * ones(N, 2);

% Solve the Hopfield model using a stochastic differential equation (SDE)
xs = solve_hopf_sde(omega, a, wC, dt, Tmax, TR, sig);

% Format and normalize the output BOLD signals
BOLD = xs';
BOLD = zscore(BOLD, 0, 2);
end

