function rdcm = run_rdcm(BOLD, TR)
% RUN_RDCM Estimates the resting-state effective connectivity using DCM for fMRI.
%
% BOLD: BOLD signal time series
% TR  : Repetition time (sampling interval) of the fMRI data

% Prepare the input data structure
Y.y = BOLD;
Y.dt = TR;

% Model specification for the resting-state DCM
rdcm = tapas_rdcm_model_specification(Y, [], []);

% Estimate the effective connectivity (A) using the DCM for fMRI
[output, ~] = tapas_rdcm_estimate(rdcm, 'r', [], 1);

% Extract the estimated connectivity matrix and remove the diagonal elements
rdcm = output.Ep.A;
rdcm = rdcm - diag(diag(rdcm));
end