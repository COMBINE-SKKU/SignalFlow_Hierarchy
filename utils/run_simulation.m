function BOLD = run_simulation(J, varargin)
% RUN_SIMULATION Simulates BOLD signals using the specified model
% given a connectivity matrix J and type of model: 'Hopf', 'Linear', or 'Non-Linear'.
% If 'Hopf' model is chosen, peak frequency vector f_diff is required.
% Global coupling parameter G is optional; if not provided, default is 2.

% Normalize the structural connectivity matrix (SC)
J = J-diag(diag(J)); % Remove self-connections (diagonal elements)

% Set simulation parameters
sig = 0.04; % Noise level
dt = 0.1; % Time step
TR = 1; % Repetition time (sample rate of BOLD signal)
Tmax = 2400; % Total simulation time

% Check if G is provided, if not, set default to 2
if isempty(varargin)
    G = 2;
else
    G = varargin{1};
end

% Check if f_diff is provided
if length(varargin) < 2
    error('The peak frequency vector f_diff is required for Hopf model.')
end

f_diff = varargin{2};

% Set parameters for Hopf model
N = size(J, 1);
wC = G * J; % Scale the connectivity matrix

% Define node frequencies and convert them to radians
omega = repmat(2 * pi * f_diff, 2, 1);

% Set intrinsic coupling strength
a = -0.02 * ones(N, 2);

% Solve the Hopfield model using a stochastic differential equation (SDE)
xs = solve_hopf_sde(omega', a, wC, dt, Tmax, TR, sig);

% Format and normalize the output BOLD signals
BOLD = xs';
BOLD = zscore(BOLD, 0, 2);
end