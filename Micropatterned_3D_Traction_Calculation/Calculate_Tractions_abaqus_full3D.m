%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                   Calculate the 3D tractions exerted by a cell inside a micro-patterned well on a soft gel.                  %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%   Last Revison Date: 23/01/2024                                                                                              %
%   Author: Manuel Gomez Gonzalez                                                                                              %
%                                                                                                                              %
%   References:                                                                                                                %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% The parameter "exp_type" is appended to some results files. It can be used to set and group different parameters to different 
% experiment types.
exp_type = 'Wells';

% Number of threads for the parfors and abaqus.
n_th = 8;
n_th_abaqus = 8;

modify_workers_local_cluster(n_th, n_th);

% Initial time.
t_0 = 0;

% List of positions to analyze.
f_0 = 5;
f_f = f_0;

f_list = f_0:f_f;

% File.pathname is created by Matlab's PIV program, and it states where the PIV analysis data was saved. However, we might have 
% moved folders before performing the Abaqus analysis. If so, write in fold the folder where the PIV data currently resides.
fold = '/mnt/Data/Laura/2022-07-23_TFM_15kPa_15-19well_DMSO-Bleb/cells/220723_Results-1/15well_gel2_T0-Try/f';

% Settings needed for abaqus caltulation.
Settings_abaqus.abaqus_path = '/usr/simulia/abaqus/6.14-5/code/bin/abq6145';    % Path of the abaqus executable.
Settings_abaqus.abaqus_calculate = 1;           % 1 to calculate with Abaqus. 0 otherwise.
Settings_abaqus.abaqus_ODB_results = 1;         % 1 to extract data from Abaqus' ODB file into a matlab file.

% Padding used for the position and time names.
Settings_new.pad = 1;

% Different experiments might have different parameters. Define and group them here.
if strcmp(exp_type, 'Wells')

    % Offset in Z used to move the reference system from the Abaqus data to Matlab's data. Abaqus reference system is located 
    % at the lowest plane, under the center of the well. It has to be moved to the surface plane, under the center of the well.
    Settings_new.zoffset = -70;

    % The units off the PIV results are in pixels and planes. In the PIV, the settings PixelSizeXY and PixelSizeZ are set to 
    % indicate the pixel size in xy and the plane step of the image stacks. In case it was not specify in the PIV, or if it was 
    % wrong, they can be overwritten here.
    % Correct the pixel size of the Settings file, if needed.
    Settings_new.PixelSizeXY = 0.0995;
    Settings_new.PixelSizeZ = 0.17;

    % If the image stack was acquired inverted, meaning the bottom part of the gel is at the top part of the stack, the PIV 
    % data has to be inverted in Z. It is controlled by the parameter Settings_abaqus.dir_inv)='Inverted'.

    % Name of the set of nodes/elements where we impose the displacements.
    Settings_abaqus.set_name = 'Set-2' ;
    Settings_abaqus.set_bottom_name = 'Set-1';

    % Constrain the tractions to the location of cells or not. Options: 'Unconstrained' or 'Constrained'.
%     Settings_new.constrained = 'Unconstrained' ;
    Settings_new.constrained = 'Constrained' ;

    % Parameters to choose the correct abaqus model. The model is stored in a directory 
    % [fold_abaqus_models, filesep, 'well_R_', num2str(Settings_new.R), '_h_', num2str(Settings_new.h)].
    fold_abaqus_models = '/home/manu/Desktop/Abaqus_Models/15kPa/D_13/Ra1_Rb0.5';   % Directory where the models are stored.
    Settings_new.R = 6;         % Radius of the well model.
    Settings_new.h = 12;        % Depth of the well model.

    abaqus_model_loc = [fold_abaqus_models, filesep, 'well_R_', num2str(Settings_new.R), '_h_', num2str(Settings_new.h)];

    % Correct the RI missmatch artifact (=1) or not (=0).
    Settings_new.Correct_Disp_RI = 1;
    % Size of the well (15 or 19). The corection at the surface is different.
    Settings_new.Well_Size = 15;

    % When finding the location of the upper plane of the gel. Multiple cases can be controlled with two parameters:
    %   Settings.auto_surf_finding: if =1 or if it is not defined, the surface will be found automatically. 
    %                               If =0, the surface plane is selected manually.
    %   Settings.beads_on_gel_surface: if =1 or if it is not defined, assume that the beads are on the surface of the gel. 
    %                                  If =0, assume that the beads are everyehere in the gel.

end

% Job created in Abaqus to analyze the FEM problem.
Settings_abaqus.jobabaqus_file = 'Job-111.inp';

% Scratch directory used by abaqus. Set it if you have any problem with the default one (e.g. if /tmp is mapped to ram, and you 
% run out of ram memory).
Settings_abaqus.abaqus_scratch_dir = '/tmp';
% Settings_abaqus.abaqus_scratch_dir = '/mnt/Data/abaqus_temp';

% Maximum memory allowed to be used by abaqus, in megabytes.
Settings_abaqus.abaqus_mem_use = num2str(90*1024);

% Loop over a list of positions.
for ii = 1:numel(f_list)

    ff = f_list(ii);        % Position to analyze.

    close all; clc;

    disp(['Analyzing position ', num2str(ff), '...']);

    folder = [fold, sprintf(['%0', num2str(Settings_new.pad), 'd'], ff)];       % Folder to analyze.

    % Load the PIV's File.mat and Settings.mat.
    clear('File');
    clear('Settings');

    load([folder, '/File.mat']);
    load([folder, '/Settings.mat']);

    % Update the paths, in case the PIV data has changed location.
    File.pathname = folder;

    File.AbaqusPath = [File.pathname, filesep, 'abaqusdata'];
    File.BoundaryPath = [File.pathname, filesep, 'Boundary'];
    File.CroppPath = [File.pathname, filesep, 'Croppeddata'];
    File.DispPath = [File.pathname, filesep, 'Displacements'];
    File.ImPath = [File.pathname, filesep, 'Images'];
    File.MaskPath = [File.pathname, filesep, 'Mask'];
    File.TracPath = [File.pathname, filesep, 'Tractions_Abaqus'];

    % Update the Settings, if we specify new values.
    for fld = fieldnames(Settings_new)'
        Settings.(fld{1}) = Settings_new.(fld{1});
    end

    % Save the new File.mat, Settings.mat and Seeettings_abaqus.mat.
    save([File.pathname, filesep, 'File.mat'], 'File', '-mat');
    save([File.pathname, filesep, 'Settings.mat'], 'Settings', '-mat');
    save([File.pathname, filesep, 'Settings_abaqus.mat'], 'Settings_abaqus', '-mat');

    % Create the new directories.
    for new_dir = {File.AbaqusPath, File.BoundaryPath, File.CroppPath, File.DispPath, File.ImPath, File.MaskPath, File.TracPath}

        if ~exist(new_dir{1}, 'dir')
            mkdir(new_dir{1});
        end

    end

    copyfile([abaqus_model_loc, filesep, '*'], File.AbaqusPath);

    % Loop over all the timepoints.
    for tt = t_0:File.NFiles.Beads-1

        File = Tractions_abaqus_full3D(Settings, File, exp_type, Settings_abaqus, ff, tt, t_0, n_th, n_th_abaqus);

    end

end