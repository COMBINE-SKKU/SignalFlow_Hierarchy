function CC = FC_approx_Hopf(J, G, FC_target, atlas, state)
% FC_APPROX_HOPF Computes the functional connectivity approximation using the Hopf model.
%
% J          : Input connectivity matrix
% G          : Scaling factor for the connectivity matrix
% TargetBOLD : Target BOLD series
% atlas      : String representing the atlas to use ('Schaefer', 'MMP', or 'Markov')
% type       : String indicating the type of model ('Hopf', 'Linear', or 'Non-Linear')

N = size(J,1);

% Load the peak frequency data based on the chosen atlas
switch atlas
    case 'Schaefer'
        f_diff = importdata(fullfile('../../data/Schaefer_peak_freq.mat'));
    case 'MMP'
        switch state
            case 'resting'
                f_diff = importdata(fullfile('../../data/MMP_peak_freq.mat'));
            case 'movie'
                f_diff = importdata(fullfile('../../data/MMP_movie_peak_freq.mat'));
            case 'pain'
                f_diff = importdata(fullfile('../../data/MMP_pain_peak_freq.mat'));
        end
        
end

f_diff = f_diff(1:N);

% Run the simulation based on the chosen model
SimulatedBOLD = run_simulation(J, G, f_diff);

% Process the BOLD signal
SimulatedBOLD = zscore(SimulatedBOLD, 0, 2)';

% Calculate static FC
FC_simulated = corr(SimulatedBOLD);
FC_simulated = FC_simulated - diag(diag(FC_simulated));
CC = calculate_correlation(FC_target, FC_simulated,1,0,1,'Pearson');
end