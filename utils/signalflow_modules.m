function [signal_flow_neg, signal_flow_pos] = signalflow_modules(EC, roi, module_type)
% SIGNALFLOW_MODULES Computes the negative and positive signal flow for a specified ROI (region of interest) 
% and module type.
%
%   [signal_flow_neg, signal_flow_pos] = signalflow_modules(EC, roi, module_type)
%
%   Inputs:
%       EC          - A square matrix representing effective connectivity (EC) data.
%       roi         - The region of interest for which the signal flow is computed.
%       module_type - A string specifying the module type ('22modules' or '26modules').
%
%   Outputs:
%       signal_flow_neg - A vector containing the mean negative signal flow for each module.
%       signal_flow_pos - A vector containing the mean positive signal flow for each module.

    % Load the appropriate module information based on the module_type input
    switch module_type
        case '22modules'
            modules = importdata('/combinelab/03_user/younghyun/01_project/01_HierarchyMapping/data/module.txt');
            num_modules = 22;
        case '27modules'
            modules = importdata('/combinelab/03_user/younghyun/01_project/01_HierarchyMapping/data/module27(ordered).mat');
            num_modules = 27;
        case 'nomodule'
            modules = 1:360;
            num_modules = 360;
        otherwise
            error('Invalid module type. Choose either "22modules" or "26modules".');
    end

    % Identify the indices corresponding to the given ROI in the modules
    seed = find(modules == roi);

    % Create a state vector with a 1 for the seed ROI and 0 elsewhere
    z = zeros(360, 1);
    z(seed) = 1;

    % Normalize the EC matrix and remove self-connections
    EC = normalize_connectivity(EC, 0);
    EC = EC - eye(size(EC, 1));

    % Compute the energy flow for the system
    signal = energy_flow(EC, z);

    % Initialize vectors to store the mean positive and negative signal flow
    signal_flow_neg = zeros(num_modules, 1);
    signal_flow_pos = zeros(num_modules, 1);

    % Iterate through each ROI to compute signal flow
    for current_roi = 1:num_modules
        % Find the indices for the current ROI in the modules
        seed = find(modules == current_roi);

        % Compute the mean positive signal flow for the current ROI
        signal_flow_pos(current_roi) = mean(max(signal(100, seed), 0));

        % Compute the mean negative signal flow for the current ROI
        signal_flow_neg(current_roi) = mean(min(signal(100, seed), 0));
    end

end
