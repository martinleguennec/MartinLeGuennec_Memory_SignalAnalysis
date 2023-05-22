%% Short description
%   Calculates the RMS enveloppe of the EMG signals
%
% Calling sequence
%   run("rmsCalculation.m")
%
% Parameters
%   Data             : structure,  contains all the EMG values, created 
%                                  with loadFiles.m function
%   rmsMovingWindow  : number,     determined at the beginning of the
%                                  script ; determines the length (in ms)
%                                  of the window of the RMS enveloppe
%
% Description
%   rmsCalculation loops through the EMG files of sets of all subjects and
%   calculates, for each repetition, the RMS enveloppe of the EMG signal
%   with a moving window whom length is determined at the beginning of the
%   code
%   The RMS enveloppe of each subject is normalized with respect to the 
%   maximum value of the RMS enveloppe of the PRE set that he performed
%   (consisting of 3 repetitions)

% Authors
%   Martin Le Guennec - Univ. Montpellier - France
%
% Versions
%   1.0.0 -- M. Le Guennec -- 2023-05-03
%     First version
%   1.0.1 -- M. Le Guennec -- 2023-05-21
%     Clearing the code to be easier to understand
%     Add the documentation

%% Inform the user 

disp(" "); disp("Calculating the RMS of the EMG signals")
disp("  Please wait, it might take a few minutes")

%% Executing code 

rmsMovingWindow = 0.250;  % We want to apply the RMS over 250 ms window

% Prepare a list of the name of the muscles
nameMuscles = {'GaMe', 'GaLa', 'SeTe', 'BiFe', 'VaMe', 'ReFe', 'VaLa', 'GlMa', 'LuEx'};


for nSession = 1:length(Data.session)

    disp(" "); disp("  Session n°" + nSession)

    for nSubject = 1:length(Data.session(nSession).subject)

         disp("    -> Subject " + nSubject)
        
        % First thing we need to do is find the set named "PRE"
        % We want to normalize the RMS by the maximal value of the RMS
        % enveloppe of the PRE set, so we need to calculate the RMS of this
        % set first, before calculating the RMS enveloppe of the rest of
        % the sets performed
        nSetsOfSubject = size(Data.session(nSession).subject(nSubject).set, 2);
        
        for nSet = 1:nSetsOfSubject
            setName = Data.session(nSession).subject(nSubject).set(nSet).name;
            if isempty(setName) == 0  % Some sets are missing
                if setName == "v1mspre"
                    iPre = nSet;  % iPre now contains the row number of PRE
                elseif setName == "PRE"
                    iPre = nSet;
                end
            end
        end

        % Now we can take the EMG data from the set and calculate its RMS
        EMG_pre = table2array(Data.session(nSession).subject(nSubject).set(iPre).EMG.data);
        timeEMG_Pre = EMG_pre(:, 1);
        EMG_pre = EMG_pre(:, 2:end);

        % rms function only take one measurement at a time
        for nMuscle = 1:size(EMG_pre, 2)
            EMG_preHp   = buttHighPass(SAMP_FREQ_EMG, 20, EMG_pre(:, nMuscle));
            EMG_preHpLp = buttLowPass(SAMP_FREQ_EMG, 450, EMG_preHp);

            [RMS, ~] = rms(EMG_preHpLp, SAMP_FREQ_EMG, rmsMovingWindow);
            maxRMSpre(nMuscle) = max(RMS);
        end
        
        % Now, we can loop through the EMG measurement of the subject
        for nSet = 1 : length(Data.session(nSession).subject(nSubject).set)

            setName = Data.session(nSession).subject(nSubject).set(nSet).name;
            if isempty(setName) == 0 && setName ~= "baseline" 

                EMG_set = table2array(Data.session(nSession).subject(nSubject).set(nSet).EMG.data);
                EMG_set = EMG_set(:, 2:end);

                % rms function only take one measurement at a time
                for nMuscle = 1:size(EMG_pre, 2)

                    % We begin by filtering the EMG signal
                    % we filter it between 20 and 450 Hz to follow the
                    % recommendations of Kamaruddin et al. (2015)
                    % The low pass filter is designed to get rid of 
                    % movement artifacts whereas the high pass is designed
                    % to get rid of signal aliasing
                    EMG_setHp   = buttHighPass(SAMP_FREQ_EMG, 20, EMG_set(:, nMuscle));
                    EMG_setHpLp = buttLowPass(SAMP_FREQ_EMG, 450, EMG_setHp);

                    % Load the matrix containing the squat phases indices
                    repetitionPhases = Data.session(nSession).subject(nSubject).set(nSet).EMG.phases;

                    for nRepetition = 1:size(repetitionPhases, 1)
                        % Prepare the indices from the beginning and ending
                        % of the repetition
                        iBegRepetition = repetitionPhases(nRepetition, 1);
                        iEndRepetition = repetitionPhases(nRepetition, 3);
                
                        % Cut the signal then claculate its RMS
                        EMG_rep = EMG_setHpLp(iBegRepetition : iEndRepetition);
                        [RMS, timeRMS] = rms(EMG_rep, SAMP_FREQ_EMG, rmsMovingWindow);
                        RMS = RMS./ maxRMSpre(nMuscle);  % The RMS is normalized with regards to the pre values calculated before

                        %%%%%%%%%%%%%%%% Extract variables %%%%%%%%%%%%%%%%

                        % For the area under the curve, we want to
                        % calculate the integer then sum all the values
                        aucRMS = cumtrapz(timeRMS, RMS);
                        aucRMS = @(a,b) max(aucRMS(timeRMS <= b)) - min(aucRMS(timeRMS >= a)); 
                        aucRMS = aucRMS(timeRMS(1), timeRMS(end));

                        % The max and mean have a function already existing
                        maxRMS = max(RMS);
                        meanRMS = mean(RMS);

                        % The time to maximum value is also easy to find
                        t2max = timeRMS(find(RMS == maxRMS));

                        % Save them 
                        Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(nMuscle).muscle = nameMuscles(nMuscle);
                        Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(nMuscle).meanRMS = meanRMS;
                        Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(nMuscle).maxRMS = maxRMS;
                        Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(nMuscle).t2maxRMS = t2max;
                        Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(nMuscle).aucRMS = aucRMS;
                    end
                end
            end
        end
    end
end

%% Clear the unused variables to keep a clear workspace 
clear EMG_pre EMG_preHp EMG_preHpLp EMG_rep EMG_set EMG_setHp EMG_setHpLp
clear timeEMG_pre timeEMG_set timeRMS
clear iBegRepetition iEndRepetition iPre
clear maxRMS maxRMSpre meanRMS nameMuscles t2max
clear nMuscle nRepetition nSession nSet nSetsOfSubject nSubject setName
clear repetitionPhases RMS rmsMovingWindow

%% Inform the user
disp(" "); disp("RMS calculated for all files"); disp(" ")


%% Function rms

function [RMS, time_RMS] = rms(signal, sampFreq, movingWindow)
    % rms : rectify and smooth a signal with a RMS filter
    %
    % Calling Sequence
    %   [RMS, time_RMS] = rms(signal, time, sampFreq, movingWindow)
    %
    % Parameters
    %   signal        : vector,  the signal to be filtered with the rms
    %   sampFreq      : number,  the sampling frequency of the signal
    %   movingWindow  : number,  the time epoch (in s) over which the rms
    %                            the RMS enveloppe will be calculated
    %
    % Output
    %   RMS           : vector,  the filtered signal
    %   time_RMS      : vector,  same size as RMS vector, time
    %                            corresponding to RMS output
    %
    % Description
    %   rms : calculates the rms over movingWindow * sampFreq points to
    %   rectify and smooth the signal
    
    % Authors
    %  Martin Le Guennec - Univ. Montpellier - France
    %
    % Versions
    %  Version 1.0.0 -- M. Le Guennec -- 2023-05-19
    %    First version

    % Verify that signal is a line
    if size(signal, 1) > size(signal, 2)
        signal = signal';
    end
    
    % Determine the number of points over which the RMS will be calculated
    nPointMovingWindow = movingWindow .* sampFreq;
    
    % Create time signal from the sampling frequency
    time = linspace(0, numel(signal)/sampFreq, numel(signal));

    % Create empty vectors 
    sizeRMS = numel(signal) - nPointMovingWindow;
    RMS = zeros((sizeRMS - nPointMovingWindow), 1);
    time_RMS = RMS;
    
    for nPoint = 1 : sizeRMS
    
        iBeg = nPoint;
        iEnd = iBeg + nPointMovingWindow;
    
        x = signal(iBeg : iEnd);
        rms = sqrt((sum(x.^2))./length(x));
        RMS(nPoint) = rms;
        time_RMS(nPoint) = mean(time(iBeg : iEnd));
    end

end

