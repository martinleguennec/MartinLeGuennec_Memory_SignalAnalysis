%% Short description
%   Calculates the instantaneous mean frequency of the EMG signals
%
% Calling sequence
%   run("waveletTransform.m")
%
% Parameters
%   Data  : structure,  contains all the EMG values, created with 
%                       loadFiles.m function
%
% Description
%   waveletTransform loops through the EMG files of sets of all subjects
%   then perform a continuous wavelet transform on each repetition signal,
%   thanks to this wavelet transform it then calculates the instantaneous
%   mean frequency of each repetition

% Authors
%   Martin Le Guennec - Univ. Montpellier - France
%
% Versions
%   1.0.0 -- M. Le Guennec -- 2023-05-10
%     First version
%   1.0.1 -- M. Le Guennec -- 2023-05-21
%     Clearing the code to be easier to understand
%     Add the documentation

%% Inform the user 
disp(" ")
disp("Calculating IMNF")
disp("  Please wait, it might take a few minutes")

%% Execute the code

figures = "N";

% Prepare a list of the name of the muscles
nameMuscles = {'GaMe', 'GaLa', 'SeTe', 'BiFe', 'VaMe', 'ReFe', 'VaLa', 'GlMa', 'LuEx'};


for nSession = 1:length(Data.session)

    disp(" "); disp("  Session n°" + nSession)

    for nSubject = 1:length(Data.session(nSession).subject)

        disp("    -> Subject " + nSubject)

        for nSet = 1 : length(Data.session(nSession).subject(nSubject).set)


            % Verify if the set exists and is not the baseline, we don't
            % interest ourselves in this measurement because there is no
            % movement during the baseline
            setName = Data.session(nSession).subject(nSubject).set(nSet).name;
            if isempty(setName) == 0 && setName ~= "baseline" 

                EMG_set = table2array(Data.session(nSession).subject(nSubject).set(nSet).EMG.data);
                timeEMG_set = EMG_set(:, 1);
                EMG_set = EMG_set(:, 2:end);

                % rms function only take one measurement at a time
                for nMuscle = 1:size(EMG_set, 2)

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
                
                        % Cut the signal then apply a wavelet transform
                        EMG_rep = EMG_setHpLp(iBegRepetition : iEndRepetition);
                        [wt, frequencies] = cwt(EMG_rep, 'amor', SAMP_FREQ_EMG, 'VoicesPerOctave', 48);

                        % Calculate the instantaneous mean frequency of the
                        % repetition
                        IMNF = imnf(wt, frequencies);

                        if figures == "Y"
                            figure(1); clf;
                            cwt(EMG_rep, 'amor', SAMP_FREQ_EMG, 'VoicesPerOctave', 48)
                            hold on
                            yline(IMNF)
                        end
                        
                        % Save the variable in the structure
                        Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).IMNF(nMuscle).muscle = nameMuscles(nMuscle);
                        Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).IMNF(nMuscle).value = IMNF;
                    end
                end
            end
        end
    end
end

%% Clear the workspace
clear nameMuscles nMuscle nRepetiton nSession nSet nSubject setName ...
    EMG_set EMG_setHp EMG_setHpLp EMG_rep timeEMG_Pre timeEMG_set ...
    IMNF figures repetitionsPhases ...
    frequencies wt repetitionPhases iBegRepetition iEndRepetition

%% Inform the user 
disp("  Al signal were treated")

%% imnf function
function IMNF = imnf(wt, f)
    % imnf : calculates the instantaneous median frequency of signal 
    %
    % Calling Sequence
    %   IMNF = imnf(wt, f)
    %
    % Parameters
    %   wt    : matrix,  dimension Na x N, where Na is the number of scales
    %                    and N is the number of samples in x, contains the
    %                    values calculated with the cwt function on the
    %                    signal
    %   f     : number,  dimension Na x 1, Scale-to-frequency (in Hz) 
    %                    conversions of the CWT
    %
    % Output
    %   IMNF  : number,  the instantaneous mean frequency calculated
    %
    % Description
    %   imnf : calculates the instantaneous mean frequency of signal
    
    % Authors
    %  Martin Le Guennec - Univ. Montpellier - France
    %
    % Versions
    %  Version 1.0.0 -- M. Le Guennec -- 2023-05-19
    %    First version

    wta = abs(wt);
    
    s_max = length(f);
    s = 1:s_max;
    w0 = f(1);
    w = w0./(s(2:end));
    
    for i = 1:size(wta, 2)
        allIMNF(i) = trapz(w', w'.*wta(2:end, i))./trapz(w', wta(2:end,i));
    end

    IMNF = mean(allIMNF);
end
