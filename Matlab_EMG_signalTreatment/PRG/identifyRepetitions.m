%% Short description
%   Find indices of the beginning and ending of concentric and eccentric
%   phases on the position them EMG data from the Data structure
%
% Calling sequence
%   run("identifyRepetitions")
%
% Parameters
%   Data  : structure,  contains all the LPT and EMG values, created with
%                       loadFiles.m function
%
% Description
%   identifyRepetitions loops through the position and EMG files of sets of 
%   all subjects and identify the beginning and ending of the concentric
%   and eccentric phases
%   Once the repetitions are identified on the position data, we multiply
%   them with the sampling frequency values of both LPT and EMG to find the
%   indices on the EMG data
%   This script is organized as follow :
%     I  - Execution of the code on the files from the Data structure
%     II - Definition of the identifyReps and findPeaks variables used in
%          the part I

% Authors
%   Martin Le Guennec - Univ. Montpellier - France
%
% Versions
%   1.0.0 -- M. Le Guennec -- 2023-05-03
%     First version
%   1.0.1 -- M. Le Guennec -- 2023-05-21
%     Clearing the code to be easier to understand
%     Add the documentation


%% Execute the code

% There is possibility to plot the figures, "Y" allows to plot figures
figures = "N";


for nSession = 1:length(Data.session)
    for nSubject = 1:length(Data.session(nSession).subject)
        for nSet = 1 : length(Data.session(nSession).subject(nSubject).set)

            setName = Data.session(nSession).subject(nSubject).set(nSet).name;

            % There won't be any repetition to find in the baseline
            if isempty(setName) == 0 && setName ~= "baseline" 
                position = Data.session(nSession).subject(nSubject).set(nSet).LPT.data.LPT;
                timePosition = Data.session(nSession).subject(nSubject).set(nSet).LPT.data.Time;

                % Filter the position signal from lpt at 10 Hz to get rid
                % of noise due to electrical furniture of the lpt device
                % (50 Hz) and the high frequency noise
                % We filter the position because we will derivate it to get
                % the velocity afterwards
                cutoffFreq = 10;
                positionLp = buttLowPass(SAMP_FREQ_LPT, cutoffFreq, position);

                % Find the repetitions in the position signal
                repetitions = identifyReps(positionLp, SAMP_FREQ_LPT, figures);
                
                % Save in Data structure
                Data.session(nSession).subject(nSubject).set(nSet).LPT.phases = repetitions;
                
                for nRepetition = 1:size(repetitions, 1)
                    
                    % Cut position and time signal, then create velocity
                    % vector for the repetition
                    positionRep = positionLp(repetitions(nRepetition, 1) : repetitions(nRepetition, 3));
                    timeRep = timePosition(repetitions(nRepetition, 1) : repetitions(nRepetition, 3));
                    velocityRep = diff(positionRep) ./ diff(timeRep);

                    % Prepare indices for end and begin of phases
                    iBegEccentric  = 1;
                    iEndEccentric  = repetitions(nRepetition, 2) - repetitions(nRepetition, 1);
                    iEndConcentric = repetitions(nRepetition, 3) - repetitions(nRepetition, 1);

                    % Extract the variables
                    repDuration           = timeRep(end) - timeRep(1);
                    repConcentricDuration = timeRep(end) - timeRep(iEndEccentric);
                    ROM                   = positionRep(iBegEccentric) - positionRep(iEndEccentric);
                    MVC                   = mean(velocityRep(iEndEccentric : iEndConcentric));
                    peakVelocity          = max(velocityRep(iEndEccentric : iEndConcentric));
                    t2p                   = timeRep(find(velocityRep == peakVelocity)) - timeRep(iEndEccentric);

                    % Save them in the structure
                    Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).repDuration = repDuration;
                    Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).repConcentricDuration = repConcentricDuration;
                    Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).ROM = ROM;
                    Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).MVC = MVC;
                    Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).peakVelocity = peakVelocity;
                    Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).t2p = t2p;
                end
                
                % We also want to have the begin and end of each phase but
                % for the EMG data
                % We will use the indices we get from the the LPT and
                % multiply them given that we have both sampling frequency
                % this should be easy
                repetitionsEMG = round((repetitions .* SAMP_FREQ_EMG) ./ SAMP_FREQ_LPT);

                % Save in the Data structure
                Data.session(nSession).subject(nSubject).set(nSet).EMG.phases = repetitionsEMG;
                
                % This allows a  quick visualization to verify if all the
                % graphs are good
                if figures == "Y"
                    % Prepare some data
                    EMG = (Data.session(nSession).subject(nSubject).set(nSet).EMG.data.RGLUTEUSMAXIMUS_EMG14_V_);
                    time =  (Data.session(nSession).subject(nSubject).set(nSet).EMG.data.X_s_);
    
                    % Make the graph
                    figure(1); clf;
                    sgtitle("Identification of eccentric (in green) and concentric (in red) phases")
                    plot(time, EMG); hold on
                    for rep = 1:size(repetitionsEMG, 1)
                        plot(time(repetitionsEMG(rep,1) : repetitionsEMG(rep,2)), EMG(repetitionsEMG(rep,1) : repetitionsEMG(rep,2)), "g")
                        plot(time(repetitionsEMG(rep,2) : repetitionsEMG(rep,3)), EMG(repetitionsEMG(rep,2) : repetitionsEMG(rep,3)), "r")
                    end
                    hold off
                    % Make the graph look good
                    title("On the position signal");
                    xlabel("Time (s)"); ylabel ("Position (m)"); 
                end
            end
        end
    end
end

%% Clear the unused variables to keep a clear workspace 
clear nRepetition nSession nSet nSubject rep setName                         % The number of the loops
clear figures cutoffFreq iBegEccentric iEndEccentric iEndConcentric          % The different indices 
clear time position positionLp positionRep timePosition timeRep velocityRep  % The different position, velocity and time vectors
clear MVC repConcentricDuration repDuration ROM t2p peakVelocity             % The calculated variables 
clear repetitions repetitionsEMG                                             % The matrix with repetitions indices

%% Inform the user 

disp("Repetitions' phases identified for all files")
disp(" ")

%% Function identifyReps
function phases = identifyReps(position, sampFreq, figures)
    % identifyReps : identifies the beginning and ending of eccentric and
    %                concentric phases of squat movements on a signal of
    %                linear position transducer
    %
    % Calling Sequence
    %   phases = identifyReps(position, sampFreq, figures)
    %
    % Parameters
    %   position  : vector,     the position signal to analyse 
    %   sampFreq  : number,     the sampling frequency of the position
    %   figures   : character,  if "Y" displays the figures
    %
    % Output
    %   phases    : matrix,     of dimension n x 3 where n is the number of
    %                           repetitions identified, gives in the first
    %                           column the beginning of the nth eccentric
    %                           phase, in the second column the ending of
    %                           the nth eccentric phase, and in the third
    %                           column the end of the nth concentric phase.
    %                           The beginning of the concentric phase
    %                           corresponds to the end of the eccentric
    %                           phase
    %
    % Description
    %   identifyReps identifies the beginning and ending of the concentric
    %   and eccentric phases of each squat performed on the position
    %   signal. It uses the position signal and the velocity, calculated at
    %   the beginning of the function, to find the phases relative to the
    %   peaks corresponding to the minimal position (low position of the
    %   squat)
    
    % Authors
    %  Martin Le Guennec - Univ. Montpellier - France
    %
    % Versions
    %  Version 1.0.0 -- M. Le Guennec -- 2023-02-13
    %    First version
    %  Version 1.0.1 -- M. Le Guennec -- 2023-05-20
    %    Clear the code
    %    Add documentation

    % The position signal must be in line
    if size(position, 1) > size(position, 2)
        position = position';
    end

    % Creation of time vectors
    timePosition = linspace (0, numel(position) ./ sampFreq, numel(position));
    timeVelocity = timePosition(2:end);

    % We will need the velocity to identify the phases 
    velocity = diff(position) ./ diff(timePosition);
    
    % Identification of the minimum position peaks 
    iMinPeaks = findPeaks(position); 
    
    % Identification of the concentric and eccentric phases
    for nPeak = 1:numel(iMinPeaks)
    
        % The ending of the eccentric phase simply corresponds to the peak
        % The beginning of the concentric phase corresponds to the last
        % positive velocity value before the peak, afterwards the subject
        % is lowering, therefore velocity is negative
        eccentricEnd = iMinPeaks(nPeak); 
        eccentricBeg = find(velocity(1 : eccentricEnd-10) >= 0);  % We add 10 to be sure that we don't take an artifact
        eccentricBeg = eccentricBeg(end) + 1; 
      
        % The ending of the concentric phase is the first velocity value
        % after the peak that is negative (amortization phase)
        concentricEnd = find(velocity(eccentricEnd + 1 : end) <= 0); 
        concentricEnd = eccentricEnd + concentricEnd(1);
    
        phases(nPeak, :) = [eccentricBeg, eccentricEnd, concentricEnd];
    end
    
    
    % If desired, make a plot of the position signals with eccentric and
    % concentric phase of each squat
    if figures == "Y"
        figure(1); clf;
        sgtitle("Identification of eccentric (in green) and concentric (in red) phases")
            
        subplot(2, 1, 1); % First with the position values
        plot(timePosition, position); hold on
        for rep = 1:numel(iMinPeaks)
            plot(timePosition(phases(rep,1) : phases(rep,2)), position(phases(rep,1) : phases(rep,2)), "g")
            plot(timePosition(phases(rep,2) : phases(rep,3)), position(phases(rep,2) : phases(rep,3)), "r")
        end
        hold off
        % Make the graph look good
        title("On the position signal");
        xlabel("Time (s)"); ylabel ("Position (m)"); 
            
        subplot(2, 1, 2); % Then with the velocity values
        plot(timeVelocity, velocity); hold on
        yline(0, "k")
        for rep = 1:numel(iMinPeaks)
            plot(timeVelocity(phases(rep,1) : phases(rep,2)), velocity(phases(rep,1) : phases(rep,2)), "g")
            plot(timeVelocity(phases(rep,2) : phases(rep,3)), velocity(phases(rep,2) : phases(rep,3)), "r")
        end
        hold off
        % Make the graph look good
        title("On the velocity signal");
        xlabel("Time (s)"); ylabel ("Velocity (mV.s^{-1})"); 
    end

end

%% subfunction findPeaks

function iMinPeaks = findPeaks(S)
    % findPeaks : find the peaks of a signal S under a defined baseline
    %
    % Calling Sequence
    %   iMinPeaks = findPeaks(S)
    %
    % Parameters
    %   S          : vector,  the signal whom peaks we want to find
    %
    % Output
    %   iMinPeaks  : vector,  the indices of the peaks
    %
    % Description
    %   findPeaks : find the negative peaks of a signal S, the value of the
    %   peaks are inferior than a threshold fixed at the middle between a
    %   baseline and the minimal value of the signal
    
    % Authors
    %  Martin Le Guennec - Univ. Montpellier - France
    %
    % Versions
    %  Version 1.0.0 -- M. Le Guennec -- 2023-02-13
    %    First version

    % We will need the derivate of the signal to find the peaks, then we
    % will compare each value n to the n+1 value, so we create these two
    % vectors too
    dS = diff(S);
    dS1 = dS(1 : end-1);
    dS2 = dS(2 : end);
    
    % The peaks of a signal can be detected by a change of sign of the
    % derivate : if the sign becomes positive it is a negative peak
    iNegativePeaks = find(sign(dS1) < sign(dS2)) + 1;
    
    % Now we keep only the lowest peaks
    % find the index of the n lowest peaks corresponding to the n squats
    baseline = mean(S(end-1000:end));                 % The signal was more stable at the end
    threshold = baseline - (baseline - min(S)) ./ 2;  % Define a threshold at the middle between baseline and minimum value
    iMinPeaks = find(S(iNegativePeaks) < threshold);   % Only keep peaks that are under these peaks
    iMinPeaks = iNegativePeaks(iMinPeaks);
end