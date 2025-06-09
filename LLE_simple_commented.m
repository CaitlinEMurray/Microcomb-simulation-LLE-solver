%% Lugiato–Lefever Equation (LLE) Simulation
% Caitlin E. Murray, Chawaphon Prayoonyong, and Bill Corcoran, 2025
% Written in MATLAB version 2022b
%
% This script performs a simple simulation of soliton formation using the 
% LLE. It is solved using a symmetrical split-step Fourier method. It is 
% currently set to use a constant-power CW input and generates solitons by 
% varying the frequency offset between the pump laser and the microring 
% resonator (i.e., detuning).
%
% Features included:
% - Third-order dispersion (D3) is implemented.
% - Easily able to be extended to include higher-order dispersion terms or a measured 
%   dispersion profile.
%
% Features not included:
% - Raman scattering, shock noise, and thermal drift are not included 
%   in this version, but users can add these effects if desired.
%
% Parallel Execution:
% - Compatible with MATLAB's Parallel Computing Toolbox, to impliment c
%   hange the main for loop to a parfor
%
% Naming Conventions:
% - Variable names follow the conventions used in:
%   Pasquazi et al., *"Micro-combs: A novel generation of optical sources"* (2018).
%
% Other useful references:
% - Kovach et al., *"Emerging material systems for integrated optical Kerr frequency combs"* (2020)
% - Kippenberg et al., *"Dissipative Kerr solitons in optical microresonators"* (2018)

%%
clear all;
close all;
clc;

%% Fundamental constants
hbar = 1.0545718e-34;   % Reduced Planck constant (Js)
c = 299792458;          % Speed of light (m/s)
i = sqrt(-1);           % Imaginary unit

%% Basic ring parameters
FSR = 25.172e9;                 % Free Spectral Range (Hz)
f0 = c / (1550e-9);             % Central frequency (Hz)
Q = 1.8e6;                      % Quality factor

linewidth = f0 / Q;             % Linewidth (Hz)
kappa_avg = linewidth * 2 * pi; % Average linewidth (rad/s)
tr = 1 / FSR;                   % Round trip time (s)
alpha = tr * kappa_avg / 2;     % Loss coefficient (rad)
norm_t = 1 / (alpha / tr);      % Normalization factor for time (s/rad)

%% General simulation parameters
nF = 1024;              % Number of frequency modes — ensure sufficient for good resolution, 
                        % but not too many, or the FFT will start to repeat
mu = (-nF/2:1:nF/2-1);  % Mode number vector

save_step = 5000;       % Loop steps between saving
calc_step = 1;          % FSR step in equation

dt = calc_step * tr;    % Time step (s)
dt_norm = dt / norm_t;  % Normalized time step 
h = dt_norm / 2;        % Half of time step 

%% Dispersion constants, including optional AMX
kappa_all = ones(1, nF) * kappa_avg;        % Vector of kappa values, can replace with measured values (rad/s)
kappa_all = kappa_all(:).';                 % Ensure this is a row vector

Dint = zeros(1, nF);                        % Pre-allocate dispersion vector 
D2 = 42e3;                                  % Second order dispersion coefficient (Hz/mode^2), for conversion from propergation coefficent, see Material and device parameters below
D3 = 0;                                     % Third order coefficient (Hz/mode^3)
Dint = 1/2 * D2 * mu.^2 + 1/6 * D3 * mu.^3; % Integrated dispersion profile (Hz)

AMX_strength = 0 * 1e6;                     % AMX strength (Hz)
AMX_loc = 142.5;                            % AMX location (mode number)
AMX = -AMX_strength ./ (mu - AMX_loc)/4;    % AMX dispersion term (Hz)
Dint = Dint + AMX;                          % Include AMX in integrated dispersion, this can be replaced with a vector of measured values (Hz)

Dint_norm = 4 * Dint * pi / kappa_avg;      % Normalized dispersion profile
Dint_norm = Dint_norm(:).';                 % Ensure Dint_norm is a row vector
%% Plot integrated dispersion profile
plot(mu, (Dint) * 1e-6);
xlabel('Mode number (\mu)');
ylabel('Dint/(2\pi) (MHz)');
xlim([-100 100]);

%% Material and device parameters

% Some reorganisation of the equations might be needed, depending on which
% parameters are known for a given device.

% These values are used for the magnitude of noise added and to convert
% normalised power to real power.

% An alternative way to convert normalised power to real power for a device 
% with known dispersion would be from the initial primary comb spacing for 
% a given input power and convert back from there.


theta = alpha;                                              % Coupling coefficient
neff = 1.6;                                                 % Effective refractive index
n2 = 1.15e-19;                                              % Nonlinear coefficient (m^2/W)
gamma_v = 0.233;                                            % Gamma parameter for high index doped silica (w^-1m^-1)

L = c / (FSR * neff);                                       % Length of ring (m)
beta2 = D2 / (-c / neff * (FSR^2 * 2 * pi));                % Second-order propagation coefficient
                          
beta3 = neff * (3 * D2^2 - 2 * D3 * FSR * pi) / ...
        (8 * c * FSR^4 * pi^3);                             % Third-order propagation coefficient
                                                            % See Kovach et al., "Emerging material systems for integrated optical Kerr frequency combs" (2020)
                                                            % for conversion equations

Aeff = n2 * (2 * pi * f0) / (c * gamma_v);                  % Effective mode area
mode_volume = Aeff * L;                                     % Mode volume
g = hbar * (2 * pi * f0)^2 * c * n2 / neff^2 / mode_volume; % Kerr frequency shift per photon
norm_E_1 = sqrt(kappa_avg / (2 * g));                       % Normalized electric field, signal_normalised = norm_E_1 * signal_original
norm_E_2 = sqrt(gamma_v *L / alpha);                        % Normlaised electric field, signal_normalised = norm_E_2 * signal_original, where signal is in sqrt(W)

%% Noise
% Noise amplitude to be added. Currently set to add white noise to both the
% initial signal and the pump.

noise_amp = sqrt(1/2/dt_norm) / norm_E_1;    % Half a photon of noise
pump_noise_amp = noise_amp;                  % Pump Noise
sig_noise_amp =  noise_amp;                  % Signal noise


%% Power values, iterated over

P_lis = 3; % List of normalized input power

for k=1:length(P_lis) % Use 'for' to run in serries, 'parfor' for parralel however requires parallel toolbox

    %% Input power  
    % Current input power is assumed to be entirely continuous wave (CW).
    % This can be changed in the pump conditions section below.
    
    X = P_lis(k);      % Normalised input power
    S = sqrt(X);       % Electric field magnitude is the square root of power
    
    %% Simulation parameters
    % The approximate end point for the soliton regime
    end_point = pi^2 * S^2 / 8;
    
    % Stop point is set to be 20% higher than the expected end point to
    % encapsulate the full dynamics
    stop_point = end_point * 1.2;
    
    %% Detuning
    % detuning = 2 * pi / (FSR * alpha) * (f - f0) 
    % f0 is the resonance frequency, f is the laser frequency

    detuning_start = -4;         % Starting detuning value. As a general rule, this should be less than 1 - sqrt(X - 1)
    detuning_end   = stop_point; % Ending detuning is set to the stop point

    % Total simulation time, rounded to ensure an integer number of save points.
    total_time_steps = round(1000000 / save_step) * save_step;

    time = total_time_steps * calc_step * tr;
    N = total_time_steps / calc_step;                      % Loop length
    
    delta_detuning = (detuning_end - detuning_start) / N;  % Detuning increment per step
    detuning = detuning_start;                             % Initial detuning value

    sweep_speed = 2* pi * alpha * delta_detuning * FSR / (tr * calc_step); % Sweep speed based on detuning (Hz/s)
    
    %% Pump and signal initial conditions
    % CW input is currently used. Spectrum is nF times larger than signal 
    % due to MATLAB FFT scaling. If using an alternative pumping scheme,
    % the input may need to be renormalised.
    
    tE_in = repmat(S, length(mu), 1);       % CW input in time domain
    E_in = fft(tE_in, [], 1);               % Converted to frequency domain (scaled by nF)
    fE_in_original = fftshift(E_in, 1);     % Re-ordered FFT to match mode numbering (mu)
    
    %% Noise
    pump_noise = random('norm', 0, pump_noise_amp, [length(mu), 1]) + ...
                 1i * random('norm', 0, pump_noise_amp, [length(mu), 1]);
    
    fE_in = fE_in_original + pump_noise;
    
    signal = random('norm', 0, sig_noise_amp, [length(mu), 1]) + ...
             1i * random('norm', 0, sig_noise_amp, [length(mu), 1]);
    
    %% Pre-allocation of arrays for more efficient saving
    SaveSignal = zeros(round(N / save_step), length(mu));   % Array for saving signal
    SaveDetuning = zeros(round(N / save_step), 1);          % Array for saving detuning values
    
    %% Loop initialisation
    spectrum = fftshift(fft(signal, [], 1), 1);
    add      = fE_in * h;    % Amount of pump added each loop, scaled by time step

     for j=1:round(N)
        
         if rem(j,save_step)==0 
            %Store the values mid simulation

            signal = ifft(fftshift(spectrum, 1), [], 1);
            

            SaveSignal(floor(j / save_step), :) = signal;
            SaveDetuning(floor(j / save_step)) = detuning;

            % Recalculate pump_noise
            pump_noise = random('norm', 0, pump_noise_amp, [length(mu), 1]) + ...
                         1i * random('norm', 0, pump_noise_amp, [length(mu), 1]);

            fE_in = fE_in_original + pump_noise;
            add = fE_in * h;

         end

        % Increase laser wavelength via changing the detuning
        detuning = detuning + delta_detuning;
        
        % Linear part
        lin_part = (-(kappa_all / kappa_avg).' - 1i * detuning - 1i * Dint_norm.');  % kappa_all/kappa_avg is the -1 in normalized version
        exp_prop = exp(lin_part * h);
        spectrum = exp_prop .* (spectrum + add);
        % Note: spectrum energy is nF times the signal energy due to MATLAB FFT defaults
        
        % Non-linear part
        signal = ifft(fftshift(spectrum));
        signal = exp(1i * (abs(signal).^2) * dt_norm) .* signal;
        
        % Next linear part
        spectrum = fftshift(fft(signal), 1);
        spectrum = exp_prop .* (spectrum + add);


     end

    % Compute final signal
    signal = ifft(fftshift(spectrum, 1), [], 1); 
    
    % Saving the important bits
    file_name = ['Ring_Z_', num2str(X, 3)];
    file_name = strrep(file_name, '.', '_');
    file_name = [file_name, '.mat'];
    
    Save_struct = SaveClass;
    Save_struct.signal = signal(:);
    Save_struct.spectrum = spectrum(:);
    Save_struct.Save_signal = SaveSignal(:,:);
    Save_struct.Save_detuning = SaveDetuning;
    Save_struct.Dint = Dint;
    Save_struct.Q = Q;
    Save_struct.X = X;

    parsave(file_name, Save_struct);

   
end
  


