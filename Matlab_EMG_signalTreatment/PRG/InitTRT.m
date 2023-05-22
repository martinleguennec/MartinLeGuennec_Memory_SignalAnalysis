% Short description
%   Initialises the working environment
%
% Calling sequence
%   run("InitTRT.m")
%
% Parameters
%   PRG_PATH  : character,  the full path where the scripts and functions
%                           are stored
%
% Description
%   InitTRT must be called at the very begining of the main script, it  :
%     1. Clears the worspace
%     2. Creates variables with paths for the repositories with data, 
%        scripts/functions and results
%     3. Loads the scripts and functions contained in the repository PRG

% Authors
%   Denis Mottet - Univ. Montpellier - France
%   Martin Le Guennec - Univ. Montpellier - France
%
% Versions
%   Version 2.0.0 -- D. Mottet -- 2011-06-17
%     First version for Scilab
%   Version 2.1.0 -- D. Mottet -- 2019-10-10
%     recursive getd on directories within PRG_PATH 
%   Version 2.1.1 -- D. Mottet -- 2021-03-01
%     Update to scilab 6 
%   Version 2.1.2 -- M. Le Guennec -- 2023-05-03
%     Adaptation for Matlab
%
% Original source
%   https://github.com/DenisMot/ScilabDataAnalysisTemplate/blob/master/PRG/InitTRT.sce


%% Clear workspace

% Close all docked and undocked array editors
desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
Titles  = desktop.getClientTitles;
for k = 1:numel(Titles)
   Client = desktop.getClient(Titles(k));     
   if ~isempty(Client) & ...
      strcmp(char(Client.getClass.getName), 'com.mathworks.mde.array.ArrayEditor')
      Client.close();
   end
end

clear all;  % Clear the workspace
close all;  % Close all figures
clc;        % Clear the console

% Inform the user
disp('==========================================================================================')
fprintf('<strong>Initializing</strong>'); disp(" ")

%% Organize the workspace on the hard disk

PRG_PATH = pwd;                           % This is the PRG file containing the current file
cd(fullfile(PRG_PATH));
WRK_PATH = fullfile(PRG_PATH, "..");      % We go up from one folder
RES_PATH = fullfile(PRG_PATH, "../RES");  % Results repository, that is within WRK
DAT_PATH = fullfile(PRG_PATH, "../DAT");  % Data repository, that is within WRK

%% Load the function folders

addpath(genpath(PRG_PATH));

%% Inform the user

disp("  Working directory :")
disp("    " + pwd())
disp(" ")
disp("  Sub-folders loaded")
disp('==========================================================================================')