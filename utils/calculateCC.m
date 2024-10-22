function meanCC = calculateCC(EC, target, G, atlas, state)
    % CALCULATECC computes mean correlation coefficient (CC) values
    %   MEANCC = CALCULATECC(EC, TARGET, G, ATLAS, STATE) computes the mean
    %   correlation coefficient values using FC_approx_Hopf function. EC is
    %   the input connectivity matrix, TARGET is the target matrix for
    %   comparison, G is a parameter (not used in this example), ATLAS is
    %   the resolution of the brain atlas, and STATE is the state of fMRI
    %   data.
    %
    %   Example:
    %   meanCC = calculateCC(EC, targetMatrix, 1, 'Schaefer', 'resting');
    %
    %   Inputs:
    %   - EC: Input connectivity matrix (typically a 2D matrix)
    %   - TARGET: Target matrix for comparison (typically a 2D matrix)
    %   - G: Parameter (optional, not used in this example)
    %   - ATLAS: Resolution of the brain atlas (optional, can be 'MMP' or 'Schaefer')
    %   - STATE: State of fMRI data (optional, can be 'resting' or 'pain')
    %
    %   Output:
    %   - MEANCC: Mean correlation coefficient values (typically a 1D vector)

    % Initialize an array to store CC values
    CC_values = zeros(10, 1);

    % Calculate CC values in a parallel loop
    parfor i = 1:10
        CC_values(i) = FC_approx_Hopf(EC, G, target, atlas, state);
    end

    % Calculate the mean CC value
    meanCC = mean(CC_values);
end
