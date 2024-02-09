%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                      Function used to set the file and folder names of the images used for the analysis.                     %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       Settings [struct]: structure containing the different parameters for the analysis.                                     %
%       File [struct]: structure containing the file name parameters for the analysis.                                         %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       File [struct]: structure containing the file name parameters for the analysis.                                         %
%                                                                                                                              %
%   Last Revison Date: 09/02/2024                                                                                              %
%   Author: Xavier Trepat's group.                                                                                             %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars                                                                 %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function File = xOpenFiles(Settings, File)

    % Create necessary folders.
    if ~exist(File.pathname, 'dir')
        mkdir(File.pathname);
    end

    File.CroppPath = [File.pathname, filesep, 'Croppeddata'];       % Folder where the cropped images will be stored.
    File.DispPath  = [File.pathname, filesep, 'Displacements'];     % Folder where the 3D PIV results will be stored.

    if ~exist(File.CroppPath, 'dir')
        mkdir(File.CroppPath);
    end

    if ~exist(File.DispPath, 'dir')
        mkdir(File.DispPath);
    end    

    % Create or select the name of one file of each type.
    file_type = {'Beads', 'Trypsin', 'BF', 'Fluo', 'Fluo_1', 'Fluo_2'};
    file_type_long = {'beads', 'trypsin', 'phase contrast', 'fluorescence reconstruction', 'Fluo_1', 'Fluo_2'};

    for ii = 1:numel(file_type)

        % Open the "file_type" file:
        if Settings.Do.Open.(file_type{ii}) && ...
           (~isfield(File, [file_type{ii}, 'Name']) || (exist(File.([file_type{ii}, 'Name']), 'file') ~= 2))

            % If no name for the "file_type" file is provided, or if it does not correspond to a valid file, pop-up a dialog 
            % asking to select a "file_type" file.
            [Temp.([file_type{ii}, 'Name']), File.([file_type{ii}, 'Path'])] = ...
                    uigetfile([File.pathname, filesep, '*.', Settings.Imgfmt], ['Open ', file_type_long{ii}, ' image']);

            % Establish the bases of the names.
            File.Base.(file_type{ii}) = Temp.([file_type{ii}, 'Name'])(1:end-length(Settings.Imgfmt)-1);
            File.Base.(file_type{ii}) = File.Base.(file_type{ii})(1:find(isletter(File.Base.(file_type{ii})), 1, 'last'));

        elseif Settings.Do.Open.(file_type{ii}) && ...
               isfield(File, [file_type{ii}, 'Name']) && (exist(File.([file_type{ii}, 'Name']), 'file') == 2)

            % If the name for the "file_type" file is provided, and it points to a valid file, use that "file_type" file.
            [File.([file_type{ii}, 'Path']), Temp.([file_type{ii}, 'Name']), Temp.ext] = ...
                fileparts(File.([file_type{ii}, 'Name']));

            % There might be problems with non-standard extensions, such as ".ome.tif". This is fixed with these two lines.
            Temp.([file_type{ii}, 'Name']) = [Temp.([file_type{ii}, 'Name']), Temp.ext];
            Temp.([file_type{ii}, 'Name']) = Temp.([file_type{ii}, 'Name'])(1:end-length(Settings.Imgfmt)-1);

            % uigetfile ends PathName with filesep, while fileparts does not.
            File.([file_type{ii}, 'Path']) = [File.([file_type{ii}, 'Path']), filesep];

            % Establish the "file_type" of the names.
            File.Base.(file_type{ii}) = ...
                Temp.([file_type{ii}, 'Name'])(1:find(isletter(Temp.([file_type{ii}, 'Name'])), 1, 'last'));

        end

    end

    % Find all of the files for the timepoints.
    file_type = {'Beads', 'BF', 'Fluo', 'Fluo_1', 'Fluo_2'};

    for ii = 1:numel(file_type)

        % How many images form the "file_type" time series?
        if Settings.Do.Open.(file_type{ii})

            % Find all the files that start with File.Base.(file_type{ii}), have the extension Settings.Imgfmt, 
            % and sort them in natural order.
            dirTrac = dir([File.([file_type{ii}, 'Path']), File.Base.(file_type{ii}), '*.', Settings.Imgfmt]);
            dirTrac = natsortfiles({dirTrac.name});

            imNum = length(dirTrac);

            if Settings.Do.LimitNumberOfFrames
                imNum = min(imNum, Settings.Do.LimitNumberOfFrames);
            end

            File.NFiles.(file_type{ii}) = imNum;

            % Now build the names of the files:
            for jj = 1:File.NFiles.(file_type{ii})
                File.Name(jj).(file_type{ii}) = [ File.([file_type{ii}, 'Path']), dirTrac{jj}];
            end

        end

    end

    % Trypsin is different. Only use one trypsin image, the one selecter or prescribed.
    if Settings.Do.Open.Trypsin

        File.NFiles.Trypsin = 1;

        % Now build the names of the files:
        File.Name(1).Trypsin = [File.TrypsinPath, Temp.TrypsinName, '.', Settings.Imgfmt];

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%