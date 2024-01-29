%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%          Script used to calculate the 3D displacement of a soft substrate through Particle Image Velocimetry (PIV).          %
%                                                                                                                              %
%   Last Revison Date: 24/01/2024                                                                                              %
%   Based on codes created by Prof. Xavier Trepat's group.                                                                     %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%------------------------------------------------------------------------------------------------------------------------------%
%                                 Parameters regarding the file names, threads and compilation.                                %
%------------------------------------------------------------------------------------------------------------------------------%

n_th_disp = 8;                      % Number of parallel execution threads in displacement finder.
% By default, each worker of a parfor is started in singlethreaded mode. However, some of the functions that we use greatly 
% benefit from paralelism. We allow each worker of the parfor to use "n_th_per_worker" threads. Use with caution.
n_th_per_worker = 8;

% Paths of the image files used. If they exist, prompt windows to select them will not appear.
% These files names have a consistent structure, e.g. date_Illumination_magnification_location_time.tif.
pad = 1;                            % Zero-padding used in the tags. For example, pad=4 renders _t0000, _t0001, _t4204, etc.
f_0 = 1;                            % Initial location analyzed.
f_end = 1;                          % Final location analyzed.

f_list = f_0:f_end ;                % Vector that contains all the locations analyzed (they don't need to be consecutive). E.g.
                                    % f_list = [ 13, 15, 16, 19, 21 ];

ending = '_t0.';                    % Time tag of any of the images, e.g. '_t0015.'.

% Time tag of the actual trypsine image. The user might have captured several trypsin images for the same position 
% (e.g. _t0.tif, _t1.tif, etc.) and want to use a specific one for each
% position. E.g. ending_tryp = { '_t0002.', '_t0000.', '_t0003.'}.
ending_tryp = cell(1, f_end+1);

for ii=(f_0:f_end)+1
    ending_tryp{ii} = '_t0.';
end

dat = '210407';                     % Date of the experiment.
dat_tryp = '210407';                % Date of the trypsin acquisition.
mag = '40x';                        % Magnification.
fold = '/home/manu/Documents/Laura/Abaqus_Current';

% beads_col = 'Far_Red';            % Parameters of the illumination name used.
beads_col = 'Green';
% beads_col = 'Red';
% fluo_col = 'Green';
% fluo_col = 'Red';
fluo_col = 'Far_Red';
% BF = 'Bright_Field_Stack_';
% BF = 'Bright_Field_';
% BF = '_Organoids_Bright_Field_';

% cells_name = 'Organoids';
cells_name = 'Cells_';
modif = '';
% modif = 'Calyculin_';
modif_tryp = '_Try';
% modif_tryp = 'Tripsin_';

% Roots for the names of the files and folders used.
pathname = [fold, '/', dat, '_Results/gel2_t0-Try/f'];
BeadsName = [fold, filesep, dat, '/',  dat, '_15Kpa_15well_beads_gel2', '_f'];
%BFName = [fold, '/', dat, '/', dat, '_', cells_name, BF, modif, mag, '/', dat, '_', cells_name, BF, modif, mag, '_f'];
TrypsinName = [fold, filesep, dat, '/',  dat, '_15Kpa_15well_beads', modif_tryp, '_gel2_f'];
Fluo_1Name = [fold, filesep, dat, '/', dat,  '_15Kpa_15well_cell-tracker_gel2', '_f'];

%------------------------------------------------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------------------------------------------%
%                                           Setup fftw, parallelization and compile.                                           %
%------------------------------------------------------------------------------------------------------------------------------%

% Optimize a little bit the fft algorithm (it will not matter much when all dimensions are powers of 2).
fftw('planner', 'exhaustive');

if exist('fftw_wisdom.mat', 'file')
    load('fftw_wisdom.mat');
    fftw('wisdom', fftw_wisdom);
end

% Compile the C code, if needed.
compile_code();

% Set the number of contiguous parallel threads of execution.
if max(n_th_disp, n_th_per_worker) > 1
    modify_workers_local_cluster(n_th_disp, n_th_per_worker);
end

%------------------------------------------------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------------------------------------------%
%                                            Calculate displacements and tractions.                                            %
%------------------------------------------------------------------------------------------------------------------------------%

total_time = tic;

for ii=f_list

    clear('File', 'Settings');

    Settings = {};
    File.pathname = [pathname, sprintf(['%0', num2str(pad), 'd'], ii)];

    % When the stacks are already cropped, we will skip them. However, there is some information in the "File.mat" and 
    % "Settings.mat", obtained when the images where actually cropped, that we need to access.
    if exist([File.pathname, filesep, 'File.mat'], 'file') && exist([File.pathname, filesep, 'Settings.mat'], 'file')

        Settings = load([File.pathname, filesep, 'Settings.mat']);
        File = load([File.pathname, filesep, 'File.mat']);

        if isfield(File, 'File')
            File = File.File;
        end

        if isfield(Settings, 'Settings')
            Settings = Settings.Settings;
        end

    end

    Settings = set_settings(Settings);

    % Paths of the image files used. If they exist, prompt windows to select them will not appear.
    File.pathname = [pathname, sprintf(['%0', num2str(pad), 'd'], ii)];
    if exist('BeadsName', 'var')==1
        File.BeadsName = [BeadsName, sprintf(['%0', num2str(pad), 'd'], ii), ending, Settings.Imgfmt];
    end
    if exist('BFName', 'var')==1
        File.BFName = [BFName, sprintf(['%0', num2str(pad), 'd'], ii), ending, Settings.Imgfmt];
    end
    if exist('TrypsinName', 'var')==1
        File.TrypsinName = [TrypsinName, sprintf(['%0', num2str(pad), 'd'], ii), ending_tryp{ii+1}, Settings.Imgfmt];
    end
    if exist('Fluo_1Name', 'var')==1
        File.Fluo_1Name = [Fluo_1Name, sprintf(['%0', num2str(pad), 'd'], ii), ending, Settings.Imgfmt];
    end

    % Open files:
    File = xOpenFiles(Settings, File);

    % Crop Images:
    [Settings, File] = cell_cropper_3D_full(Settings, File);

    % Find Displacement Fields:
    File = displacement_finder_3D_longdisp(Settings, File, n_th_disp);

end

total_time = datestr(toc(total_time)/(60*60*24), 'HH:MM:SS')

%------------------------------------------------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------------------------------------------%
%                                                       Save parameters.                                                       %
%------------------------------------------------------------------------------------------------------------------------------%

% Save the File and Settings structures.
save([File.pathname, filesep, 'Settings.mat'], 'Settings', '-mat');
save([File.pathname, filesep, 'File.mat'], 'File', '-mat');
    
% Save the fftw wisdom.
fftw_wisdom = fftw('wisdom');
save('fftw_wisdom.mat', 'fftw_wisdom');

%------------------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
