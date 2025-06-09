%% A class for saving simulation results
% This class ensures that saved variables have consistent names across runs.
% It can be easily expanded to include additional parameters as needed.

classdef SaveClass
    properties
        signal          % Final time-domain signal
        spectrum        % Final frequency-domain spectrum
        Save_signal     % Time evolution of the signal
        Save_detuning   % Detuning values at each save point
        Save_thermal    % (Optional) Thermal drift or related parameter
        Dint            % Integrated dispersion profile
        kappa           % 2 Pi x linewidth, can be vector or scalar
        Q               % Quality factor
        X               % Normalised input power
    end
end