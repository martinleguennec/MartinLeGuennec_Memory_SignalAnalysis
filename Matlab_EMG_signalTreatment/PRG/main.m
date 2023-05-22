%% Short description
%   Initializes then executes the code
%
% Calling sequence
%   run("main.m")
%
% Parameters
%   none
%
% Description
%   Change the current folder to the path for main.m before running it
%   main initializes the workspace by executing InitTRT, then create useful
%   variables
%   It will then execute the dedicated scripts to create constants, load
%   all the data in a Data structure and perform the analysis of the EMG 
%   signals by calculating their RMS enveloppe and their IMNF
%   The wanted variables will be extracted in each script and stored in the
%   Data structure
%   Finally, the usefull variables will be stored in a table that will be
%   exported as a csv file

% Authors
%   Martin Le Guennec - Univ. Montpellier - France
%
% Versions
%   1.0.0 -- M. Le Guennec -- 2023-05-03
%     First version


%% Initialize

PRG_PATH = pwd;  % Create variable with path for programms

% Run InitTRT, it will initialize the workspace
fullPathInitTRT = fullfile(PRG_PATH, "InitTRT.m");
run(fullPathInitTRT)

% Run InitTRT, it will initialize the workspace
fullPathConstantDefinition = fullfile(PRG_PATH, "constantDefinition.m");
run(fullPathConstantDefinition)

%% Load files

Data = loadFiles(DAT_PATH);
Data = Data.Data;

% Detect the repetitions in the position and EMG files and extract some
% variables about the position
fullPathIdentifyRepetitions = fullfile(PRG_PATH, "identifyRepetitions.m");
run(fullPathIdentifyRepetitions)

%% EMG Signal treatment

% First, we want to calculate the RMS enveloppe of the signal 
fullPathRmsCalculation = fullfile(PRG_PATH, "rmsCalculation.m");
run(fullPathRmsCalculation)

% Then, we want the instantaneous mean frequencies values
fullPathWaveletTransform = fullfile(PRG_PATH, "waveletTransform.m");
run(fullPathWaveletTransform)

%% Export the data as csv file 
fullPathCreateTableForR = fullfile(PRG_PATH, "createTableForR.m");
run(fullPathCreateTableForR)


