% Short description
%   Creates the constants needed for the code
%
% Calling sequence
%   run("constantDefinition.m")
%
% Parameters
%   none
%
% Description
%   constantDefinition creates the constants that will be used in the code

% Authors
%   Martin Le Guennec - Univ. Montpellier - France
%
% Versions
%   1.0.0 -- M. Le Guennec -- 2023-05-03
%     First version

%% Inform the user 
disp("Creating constants")
disp(" ")

%% EMG
SAMP_FREQ_EMG = 2148;      % Sampling frequency (Hz) of the EMG
SAMP_PERIOD_EMG = 0.4655;  % Sampling period (s) of the LPT 

%% LPT
SAMP_FREQ_LPT = 518;       % Sampling frequency (Hz) of the LPT
SAMP_PERIOD_LPT = 0.0019;  % Sampling period (s) of the LPT

