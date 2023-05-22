function Data = loadFiles(dataDirPath)
    % loadFiles : load EMG and LPT data into a structure
    %
    % Calling Sequence
    %   listOfDir = listDirInDir(folderPath)
    %
    % Parameters
    %   dataDirPath  : string,     the path of the directory where the
    %                              data are stored
    %
    % Output
    %   Data         : structure,  containing all the data from the EMG 
    %                              files in csv format stored in dataDirPath
    %
    % Description
    %   loadFiles : loops through the directories and sub-directories of 
    %   data and loads the EMG files in csv format in a structure
    %   The data directory should be organized as follow :
    %     - Session number (01 or 02)
    %       - Subject ID
    %         - Files
    
    % Authors
    %  Martin Le Guennec - Univ. Montpellier - France
    %
    % Versions
    %  Version 1.0.0 -- M. Le Guennec -- 2023-05-03
    %    First version
    %  Version 1.1.0 -- M. Le Guennec -- 2023-05-19
    %    Save Data as .mat file to load it more quickly

    
    % Inform the user
    disp(" ")
    disp("Loading data ...")
    disp("  Please wait, it might take several minutes")
    disp(" ")

    % Verify if a file data.mat already exists
    % It is very time consumming to read the csv file, load them and
    % prepare them to be used
    % Once they were loaded once, we save the data in a .mat file that
    % we will load if it is already created
    fullPathDataStructure = fullfile(dataDirPath, "data.mat");

    if exist(fullPathDataStructure) ~= 0

        Data = load(fullPathDataStructure);

    else
    
        % List the different session files
        listSessionDir = listDirInDir(dataDirPath);
        
        % Find subject files in the session directories -----------------------
        for nSession = 1 : size(listSessionDir, 1)
        
            % Find the session name
            sessionDirFullPath = listSessionDir(nSession);
            % Split the full path and keep only the last part (session number)
            nameSession = strsplit(sessionDirFullPath, "/");
            nameSession = str2num(nameSession(end));
        
            % Create a new session in the structure
            Data.session(nSession).name = nameSession;
        
            % List the different subject directory in the session file 
            listSubjectDir = listDirInDir(listSessionDir(nSession));
        
            % Loop through the subject directories ----------------------------
            for nSubject = 1 : size(listSubjectDir, 1)
        
                % Find the subject name
                subjectDirFullPath = listSubjectDir(nSubject);
                % Split the full path and keep only the last part (subject ID)
                nameSubject = strsplit(subjectDirFullPath, "/");
                nameSubject = nameSubject(end);
        
                % Create a new subject in the structure's session
                Data.session(nSession).subject(nSubject).name = nameSubject;
        
                % List the content of the directory
                contentSubjectFolder = dir(subjectDirFullPath);
        
                % Loop through the files and find the EMG files ---------------
                for nFile = 1 : size(contentSubjectFolder, 1)
        
                    fullNameFile = contentSubjectFolder(nFile).name;  % File name with extension
                    [~, nameFile, extension] = fileparts(fullNameFile);  % File name and extension separately
                    nameFileParts = strsplit(nameFile, "_");
        
                    % Only keep EMG file in .csv format
                    if size(nameFileParts, 2) >= 3
                        if nameFileParts(3) == "EMG" && extension == ".csv"
        
                            % Make the full path name
                            currentFileFullPath = fullfile(subjectDirFullPath, fullNameFile);
        
                            % Extract the LPT and EMG data from the file
                            [currentDataLPT, currentDataEMG] = loadDelsysData(currentFileFullPath);
        
                            setName = cell2mat(nameFileParts(5));
                            % The set name can either be a number or a string
                            setNumber = str2num(setName);
                            % If setNumber is empty, it means that it is a 
                            % string
                            % Take the length of the structure Data already 
                            % existing  and add one to get an index
                            % It will be used afterwards to store the data in 
                            % the structure
                            if isempty(setNumber)
                                setNumber = length(Data.session(nSession).subject(nSubject).set) + 1;
                            end
        
                            % Save the data in the structure Data
                            Data.session(nSession).subject(nSubject).set(setNumber).name = setName;
                            Data.session(nSession).subject(nSubject).set(setNumber).EMG.data = currentDataEMG;
                            Data.session(nSession).subject(nSubject).set(setNumber).LPT.data = currentDataLPT;
                        end
                    end
                end
            end
        end

        % Save the data in a .mat file 
        save(fullPathDataStructure, "Data", "-v7.3", "-nocompression")
    end
    

    % Inform the user
    disp("  Data loaded !")
    disp(" ")

end


%% listDirInDir subfunction
function listOfDir = listDirInDir(folderPath)
    % listOfDir : find directories into directory
    %
    % Calling Sequence
    %   listOfDir = listDirInDir(folderPath)
    %
    % Parameters
    %   folderPath  : string,        the path of the directory where the
    %                                function searches the directories
    %
    % Output
    %   listOfDir   : string array,  array containing full path of all folders
    %                                contained in folderPath
    %
    % Description
    %  listOfDir : find directories into the specified directory
    
    % Authors
    %  Martin Le Guennec - Univ. Montpellier - France
    %
    % Versions
    %  Version 1.0.0 -- M. Le Guennec -- 2023-05-03
    %    First version
    
    % List folder content
    folderContent = dir(folderPath);
    
    % Loop through the folder content and identify the folders
    listOfDir = [];
    for nFolderContent = 1:size(folderContent)
    
        % Only take the folders
        if folderContent(nFolderContent).isdir == 1
    
            % Don't consider hidden folders (begin with ".")
            if folderContent(nFolderContent).name(1) ~= "."
    
                % Recreate the full path of the folder
                nameFolder = folderContent(nFolderContent).name;
                fullPathFolder = fullfile(folderPath, nameFolder);
    
                % Add the full path to the list
                listOfDir = [listOfDir;fullPathFolder];
            end
        end
    end
end


%% loadDelsysData subfunction
function [dataLPT, dataEMG] = loadDelsysData(fileName)
    % loadDelsysData : reads and sorts the data from a Delsys csv file
    %
    % Calling Sequence
    %   [dataLPT, dataEMG] = loadDelsysData(fileName)
    %
    % Parameters
    %   fileName  : string,  the full path of the delsys file, the file
    %                        should be of dimension n x 20 containing data 
    %                        from the 9 EMG sensors and the LPT, each column 
    %                        being associated with a column representing 
    %                        the sampling time
    %
    % Output
    %   dataLPT   : table,   dimension n x 2, contains time (column 1) and 
    %                        data (column 2) from the LPT
    %   dataEMG   : table,   dimension n x 10, contains time (column 1) and 
    %                        data (columns 2 to 10) from the 9 EMG sensors
    %
    % Description
    %  listOfDir : this function loads the data measured with the delsys
    %  from a csv file then creates two table : 
    %    - dataLPT : contains the raw LPT values with the first values
    %                corrected
    %    - dataEMG : contains the raw EMG values
    
    % Authors
    %  Martin Le Guennec - Univ. Montpellier - France
    %
    % Versions
    %  Version 1.0.0 -- M. Le Guennec -- 2023-05-03
    %    First version
    


    warning off
    dataDelsys = readtable(fileName);

    % ############################# LPT data #############################

    % The LPT data were recolted on the Analog Input Adapter 4
    % The associated time values are in the column X_s__2
    dataLPT = dataDelsys(:,{'X_s__2', 'AnalogInputAdapter4_Analog_A4_V_'});
    dataLPT.Properties.VariableNames = {'Time', 'LPT'};
    
    % Sometimes, last values are Nan ; don't take them into account
    iEndLPT = find(isnan(dataLPT.LPT) == 1); 
    iEndLPT = iEndLPT(1) - 500;
    dataLPT = dataLPT(1 : iEndLPT, :);

    % The first values are erroneus (they are worth 0)
    % Replace them with the first real values
    iBegLPT = find(table2array(dataLPT(:, 2)) > 1); 
    iBegLPT = iBegLPT(1);
    dataLPT(1 : iBegLPT, 2) = dataLPT(iBegLPT + 1, 2); 
    
    % Convert from mV to mm
    dataLPT.LPT = ((dataLPT.LPT * 286.3915) + 338.2548) ./ 1000;
    
    
    % ############################# EMG data #############################
    
    % Muscles from the bottom first and then upwards
    dataEMG = dataDelsys(:, { ...
        'X_s_'
        'RGASTROCNEMIUSMEDIALHEAD_EMG1_V_'
        'RGASTROCNEMIUSLATERALHEAD_EMG2_V_'
        'RSEMITENDINOSUS_EMG5_V_'
        'RBICEPSFEMORIS_EMG6_V_'
        'RVASTUSMEDIALIS_EMG10_V_'
        'RRECTUSFEMORIS_EMG11_V_'
        'RVASTUSLATERALIS_EMG13_V_'
        'RGLUTEUSMAXIMUS_EMG14_V_'
        'RTHORACOLUMBARFASCIA_EMG15_V_'
        });
end

