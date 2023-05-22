%% Short description
%   Exports the variables as csv file
%
% Calling sequence
%   run("createTableForR.m")
%
% Parameters
%   Data  : structure,  contains all the variables calculated from EMG
%                       values, created with loadFiles function, then the
%                       variables were calculated with rmsCalculation and
%                       waveletTransform
%
% Description
%   createTableForR fetches the variables calculated in the Data structure
%   and stores them in a table named data_R. This table is then exported in
%   csv format

% Authors
%   Martin Le Guennec - Univ. Montpellier - France
%
% Versions
%   1.0.0 -- M. Le Guennec -- 2023-05-11
%     First version


%% Prepare empty table that we will fill afterwards
data_R = cell2table(cell(0, 50), 'VariableNames', [ ...
    "session", "subjectName", "setName", "repetition", "MVC", ...
    "IMNF_GaMe", "IMNF_GaLa", "IMNF_SeTe", "IMNF_BiFe", "IMNF_VaMe", "IMNF_ReFe", "IMNF_VaLa", "IMNF_GlMa", "IMNF_ExLo", ...
    "meanRMS_GaMe", "meanRMS_GaLa", "meanRMS_SeTe", "meanRMS_BiFe", "meanRMS_VaMe", "meanRMS_ReFe", "meanRMS_VaLa", "meanRMS_GlMa", "meanRMS_ExLo", ...
    "aucRMS_GaMe", "aucRMS_GaLa", "aucRMS_SeTe", "aucRMS_BiFe", "aucRMS_VaMe", "aucRMS_ReFe", "aucRMS_VaLa", "aucRMS_GlMa", "aucRMS_ExLo", ...
    "maxRMS_GaMe", "maxRMS_GaLa", "maxRMS_SeTe", "maxRMS_BiFe", "maxRMS_VaMe", "maxRMS_ReFe", "maxRMS_VaLa", "maxRMS_GlMa", "maxRMS_ExLo", ...
    "t2maxRMS_GaMe", "t2maxRMS_GaLa", "t2maxRMS_SeTe", "t2maxRMS_BiFe", "t2maxRMS_VaMe", "t2maxRMS_ReFe", "t2maxRMS_VaLa", "t2maxRMS_GlMa", "t2maxRMS_ExLo"]);

%% Fill the table repetition by repetition

for nSession = 1:length(Data.session)
    for nSubject = 1 :length(Data.session(nSession).subject)
        for nSet = 1 : length(Data.session(nSession).subject(nSubject).set)

            % Verify if the set exists and is not the baseline measurement 
            setName = Data.session(nSession).subject(nSubject).set(nSet).name;
            
            if isempty(setName) == 0 && setName ~= "baseline"
                for nRepetition = 1 : length(Data.session(nSession).subject(nSubject).set(nSet).repetitions)
                    
                    % Informtations about the session, set, subject,
                    % repetition
                    session       = Data.session(nSession).name;
                    subjectName   = Data.session(nSession).subject(nSubject).name;
                    setName       = Data.session(nSession).subject(nSubject).set(nSet).name;
                    repetition    = nRepetition;

                    % Motion analysis variables 
                    MVC           = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).MVC;

                    % Frequential analysis of the EMG
                    IMNF_GaMe     = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).IMNF(1).value;
                    IMNF_GaLa     = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).IMNF(2).value;
                    IMNF_SeTe     = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).IMNF(3).value;
                    IMNF_BiFe     = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).IMNF(4).value;
                    IMNF_VaMe     = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).IMNF(5).value;
                    IMNF_ReFe     = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).IMNF(6).value;
                    IMNF_VaLa     = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).IMNF(7).value;
                    IMNF_GlMa     = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).IMNF(8).value;
                    IMNF_ExLo     = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).IMNF(9).value;

                    % Variables about the RMS enveloppe
                    meanRMS_GaMe  = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(1).meanRMS;
                    meanRMS_GaLa  = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(2).meanRMS;
                    meanRMS_SeTe  = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(3).meanRMS;
                    meanRMS_BiFe  = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(4).meanRMS;
                    meanRMS_VaMe  = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(5).meanRMS;
                    meanRMS_ReFe  = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(6).meanRMS;
                    meanRMS_VaLa  = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(7).meanRMS;
                    meanRMS_GlMa  = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(8).meanRMS;
                    meanRMS_ExLo  = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(9).meanRMS;

                    aucRMS_GaMe   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(1).aucRMS;
                    aucRMS_GaLa   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(2).aucRMS;
                    aucRMS_SeTe   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(3).aucRMS;
                    aucRMS_BiFe   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(4).aucRMS;
                    aucRMS_VaMe   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(5).aucRMS;
                    aucRMS_ReFe   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(6).aucRMS;
                    aucRMS_VaLa   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(7).aucRMS;
                    aucRMS_GlMa   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(8).aucRMS;
                    aucRMS_ExLo   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(9).aucRMS;

                    maxRMS_GaMe   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(1).maxRMS;
                    maxRMS_GaLa   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(2).maxRMS;
                    maxRMS_SeTe   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(3).maxRMS;
                    maxRMS_BiFe   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(4).maxRMS;
                    maxRMS_VaMe   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(5).maxRMS;
                    maxRMS_ReFe   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(6).maxRMS;
                    maxRMS_VaLa   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(7).maxRMS;
                    maxRMS_GlMa   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(8).maxRMS;
                    maxRMS_ExLo   = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(9).maxRMS;

                    t2maxRMS_GaMe = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(1).t2maxRMS;
                    t2maxRMS_GaLa = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(2).t2maxRMS;
                    t2maxRMS_SeTe = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(3).t2maxRMS;
                    t2maxRMS_BiFe = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(4).t2maxRMS;
                    t2maxRMS_VaMe = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(5).t2maxRMS;
                    t2maxRMS_ReFe = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(6).t2maxRMS;
                    t2maxRMS_VaLa = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(7).t2maxRMS;
                    t2maxRMS_GlMa = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(8).t2maxRMS;
                    t2maxRMS_ExLo = Data.session(nSession).subject(nSubject).set(nSet).repetitions(nRepetition).RMS(9).t2maxRMS;

            
                    % We add all these values in the table
                    add_data = cell2table({ ...
                        session, subjectName, setName, repetition, MVC, ...
                        IMNF_GaMe, IMNF_GaLa, IMNF_SeTe, IMNF_BiFe, IMNF_VaMe, IMNF_ReFe, IMNF_VaLa, IMNF_GlMa, IMNF_ExLo, ...
                        meanRMS_GaMe, meanRMS_GaLa, meanRMS_SeTe, meanRMS_BiFe, meanRMS_VaMe, meanRMS_ReFe, meanRMS_VaLa, meanRMS_GlMa, meanRMS_ExLo, ...
                        aucRMS_GaMe, aucRMS_GaLa, aucRMS_SeTe, aucRMS_BiFe, aucRMS_VaMe, aucRMS_ReFe, aucRMS_VaLa, aucRMS_GlMa, aucRMS_ExLo, ...
                        maxRMS_GaMe, maxRMS_GaLa, maxRMS_SeTe, maxRMS_BiFe, maxRMS_VaMe, maxRMS_ReFe, maxRMS_VaLa, maxRMS_GlMa, maxRMS_ExLo, ...
                        t2maxRMS_GaMe, t2maxRMS_GaLa, t2maxRMS_SeTe, t2maxRMS_BiFe, t2maxRMS_VaMe, t2maxRMS_ReFe, t2maxRMS_VaLa, t2maxRMS_GlMa, t2maxRMS_ExLo}, ...
                        'VariableNames', [ ...
                        "session", "subjectName", "setName", "repetition", "MVC", ...
                        "IMNF_GaMe", "IMNF_GaLa", "IMNF_SeTe", "IMNF_BiFe", "IMNF_VaMe", "IMNF_ReFe", "IMNF_VaLa", "IMNF_GlMa", "IMNF_ExLo", ...
                        "meanRMS_GaMe", "meanRMS_GaLa", "meanRMS_SeTe", "meanRMS_BiFe", "meanRMS_VaMe", "meanRMS_ReFe", "meanRMS_VaLa", "meanRMS_GlMa", "meanRMS_ExLo", ...
                        "aucRMS_GaMe", "aucRMS_GaLa", "aucRMS_SeTe", "aucRMS_BiFe", "aucRMS_VaMe", "aucRMS_ReFe", "aucRMS_VaLa", "aucRMS_GlMa", "aucRMS_ExLo", ...
                        "maxRMS_GaMe", "maxRMS_GaLa", "maxRMS_SeTe", "maxRMS_BiFe", "maxRMS_VaMe", "maxRMS_ReFe", "maxRMS_VaLa", "maxRMS_GlMa", "maxRMS_ExLo", ...
                        "t2maxRMS_GaMe", "t2maxRMS_GaLa", "t2maxRMS_SeTe", "t2maxRMS_BiFe", "t2maxRMS_VaMe", "t2maxRMS_ReFe", "t2maxRMS_VaLa", "t2maxRMS_GlMa", "t2maxRMS_ExLo"] ...
                        );
                    
                    data_R = [data_R ; add_data];
                end
            end
        end
    end
end

clear add_data session subjectName setName repetition MVC ...
    nRepetition nSession nSet nSubject ...
    IMNF_GaMe IMNF_GaLa IMNF_SeTe IMNF_BiFe IMNF_VaMe IMNF_ReFe IMNF_VaLa IMNF_GlMa IMNF_ExLo ...
    meanRMS_GaMe meanRMS_GaLa meanRMS_SeTe meanRMS_BiFe meanRMS_VaMe meanRMS_ReFe meanRMS_VaLa meanRMS_GlMa meanRMS_ExLo ...
    aucRMS_GaMe aucRMS_GaLa aucRMS_SeTe aucRMS_BiFe aucRMS_VaMe aucRMS_ReFe aucRMS_VaLa aucRMS_GlMa aucRMS_ExLo ...
    maxRMS_GaMe maxRMS_GaLa maxRMS_SeTe maxRMS_BiFe maxRMS_VaMe maxRMS_ReFe maxRMS_VaLa maxRMS_GlMa maxRMS_ExLo ...
    t2maxRMS_GaMe t2maxRMS_GaLa t2maxRMS_SeTe t2maxRMS_BiFe t2maxRMS_VaMe t2maxRMS_ReFe t2maxRMS_VaLa t2maxRMS_GlMa t2maxRMS_ExLo

%% Save the table

fullFileName = fullfile(RES_PATH, "data_R.csv");

writetable(data_R, fullFileName)

%% Inform the user
disp(" "); disp("The variables are successfully export as csv file")
disp(" ")
disp("All the analysis was succesfully conducted, you can now close Matlab")

