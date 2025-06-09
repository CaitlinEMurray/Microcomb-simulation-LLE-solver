%% A quick file to load a saved file, plot the solitons formed, and see the time/frequency domain evolution
clear all; close all; clc;
%% Load file
load('Ring_Z_3.mat')

%% Extract saved data
SaveDetuning = save_obj.Save_detuning;
SaveSignal = save_obj.Save_signal;
P = save_obj.X;

%% Mode number setup
nF = length(SaveSignal(1,:)); 
mu = (-nF/2:1:nF/2-1);     % Mode numbers

%% Find soliton point
soliton_detuning = P;  % Take spectrum when detuning = normalised power, this (typically) is in the soliton regeme
[val, soliton_point_loc] = min(abs(SaveDetuning - soliton_detuning));

spectrum = fftshift(fft(SaveSignal(soliton_point_loc,:))) / nF;
signal   = abs(SaveSignal(soliton_point_loc,:));

%% Plot: Time domain snapshot
subplot(3,1,1)
hold on
plot(linspace(-pi, pi, nF), signal)
xlabel('Position in ring (radians)')
ylabel('Power (a.u.)')
xticks([-pi, 0, pi])
xticklabels({'-\pi', '0', '\pi'})
xlim([-pi, pi])
yticks([])

%% Plot: Frequency domain snapshot
subplot(3,1,2)
hold on
stem(mu, db(spectrum) - max(db(spectrum)), 'Marker', 'none', 'BaseValue', -100);
ylim([-50, 0])
xlabel('Mode number (\mu)')
ylabel('Power (dB)')
xlim([-100 100])

%% Plot: Intracavity power vs detuning
subplot(3,1,3)
hold on
intracavity_power = sum(abs(SaveSignal).^2, 2) / nF;
plot(SaveDetuning, intracavity_power)
plot(SaveDetuning(soliton_point_loc), intracavity_power(soliton_point_loc), '^')
legend(["Intracavity power", "Point spectrum taken at"],Location="northwest")
yticks([]) 
ylabel('Intracavity power (a.u.)')
xlabel('Detuning')

%% Surface plot: Time domain evolution
[X, Y] = meshgrid(mu, SaveDetuning);

figure;
surf(X/nF * 2 * pi, Y, abs(SaveSignal), 'EdgeColor', 'none');
view(2)
shading interp
colormap jet
cbar = colorbar;
ylabel(cbar, 'Power (a.u.)')  

xlabel('Position in ring (radians)')
xticks([-pi, 0, pi])
xticklabels({'-\pi', '0', '\pi'})
xlim([-pi, pi])
ylim([SaveDetuning(1), SaveDetuning(end)])
ylabel('Detuning')
title('Time domain')


%% Surface plot: Frequency domain evolution
figure;
surf(X, Y, db(abs(fftshift(fft(SaveSignal, [], 2) / nF, 2))), 'EdgeColor', 'none');
view(2)
shading interp
colormap jet
cbar = colorbar;
ylabel(cbar, 'Power (dB)') 

xlabel('Mode number (\mu)')
xlim([-nF/2 nF/2])
ylim([SaveDetuning(1), SaveDetuning(end)])
ylabel('Detuning')
title('Frequency domain')
