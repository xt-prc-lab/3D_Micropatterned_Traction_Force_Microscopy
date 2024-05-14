%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                                 Function used to apply a 3D PIV to a pair of 3D stack images.                                %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       Settings [structure]: structure containing the different parameters for the analysis.                                  %
%       File [structure]: structure containing the parameters related to paths and file names.                                 %
%       n_th_disp [int scalar]: maximum number of threads used in the parallelization of the PIV.                              %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       File [structure]: structure containing the parameters related to paths and file names.                                 %
%                                                                                                                              %
%   Last Revison Date: 03/04/2024                                                                                              %
%   Based on a 2D PIV created by Prof. Xavier Trepat's group: Raimon Sunyer and Xavier Trepat.                                 %
%   3D PIV by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                  %
%                                                                                                                              %
%   References:                                                                                                                %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function File = displacement_finder_3D_longdisp(Settings, File, n_th_disp)

    disp(' ');
    disp('       *** Start displacement finder ***       ');
    disp(' ');

    % Import local variables from File and Settings structures.

    nBeadsFile = File.NFiles.Beads;             % Number of timepoints to analyze.
    CroppPath = File.CroppPath;                 % Location of the cropped files.
    DispPath = File.DispPath;                   % Location of the displacement files.
    namebeads = {File.Name.Beads};              % Base of the beads image file names.
    resolution = Settings.Resolution;           % XY-resolution of the PIV.
    resolution_z = Settings.Resolution_z;       % Z-resolution of the PIV.
    blocksizedec = Settings.Blocksizedec;       % Amount by which each side of the interrogation window shrinks at every 
    blocksizedec_z = Settings.Blocksizedec_z;   % iteration, in xy and z respectively.
    overlap = Settings.Overlap;                 % Overlap of the PIV boxes. 0 means no overlap. 0.5 means overlap of half box.
    overlap_z = Settings.Overlap_z;             % In xy and z respectively.

    padding = Settings.padding;                 % Padding in the X and Y directions for the PIV boxes.
    paddingz = Settings.paddingz;               % Padding in the Z direction for the PIV boxes.

    fft_window = Settings.fft_window;           % Window for FFT calculation. 1 means hanning.
    corr_method = Settings.corr_method;         % Method to compute max of cross-correlation: 1 for Centroid, 2 for Gaussian fit.
    pk_window = Settings.pk_window;             % Size of the sub-window for the peak of cross-correlation subpixel location.
    iter_max = Settings.iter_max;               % Maximum number of iterations after decrease of box size.
    iter_epsilon = Settings.iter_epsilon;       % epsilon, maximum acceptable relative error of the iteration.
    iter_N_conver = Settings.N_conver;          % N_conver: stop after N_conver iterations without improving the convergence.
    longdisp_del = Settings.longdisp_neigh;	    % del: for each box where longdisp is required, apply longdisp to the 
    longdisp_del_z = Settings.longdisp_neigh_z; % neighbouring "del" boxes around it. In XY and Z respectively.
    dynamic_window = Settings.dynamic_window;   % 1 if the sizes of the PIV windows are determined dynamically; 0 otherwise.
    filt_size = Settings.filt_size;             % Size of the Predictor filter, in window lengths, i.e. 
    filt_size_z = Settings.filt_size_z;         % filter_size = filt_size*window_size. In XY and Z respectively.
    filt_pow = Settings.filt_pow;               % Power of the Predictor filter.

    % Load the reference (trypsin) 3D stack image. It is common for the whole experiment. Load it only once.
    im1 = double(TIFFStack([CroppPath, filesep, 'crop_trypsin_', num2str(1), '.tif']));

    % When the stacks are already analyzed, we will skip them. However, there is some information in the "File.mat", 
    % obtained when the images where actually analyzed, that we need to access.
    if exist([File.pathname, filesep, 'File_old.mat'], 'file') == 2

        File_old = load([File.pathname, filesep, 'File_old.mat']);

        if isfield(File_old, 'File_old')
            File_old = File_old.File_old;
        end

        n_elem_old = numel(File_old.Name);

    end

    re_calculate = ~exist([File.pathname, filesep, 'File_old.mat'], 'file');

    % Loop over all the timepoints.
    for ii = 1:nBeadsFile

        [~, temp, ext] = fileparts(namebeads{ii});

        % There might be problems with non-standard extensions, such as ".ome.tif". This is fixed with these two lines.
        temp = [temp, ext];
        temp = temp(1:end-length(Settings.Imgfmt)-1);

        beads_file = [DispPath, filesep, 'GelDisp_', temp, '.dat'];

        if exist([CroppPath, filesep, 'crop_Beads_', num2str(ii), '.tif'], 'file') && (~exist(beads_file, 'file') || ...
           re_calculate || (ii > n_elem_old))

        	disp(['Processing crop_Beads_', num2str(ii),'.tif', ' of ', num2str(nBeadsFile), ' images']);

            % 3D-stack to process (beads image).
            im2 = double(TIFFStack([CroppPath, filesep, 'crop_Beads_', num2str(ii), '.tif']));

            [x, y, z, dx, dy, dz, pkh, conver, iter_time, ws_xy0, ws_z0, Convergence_Mask] = disp_on_blocks_longdisp_3D_parfor(...
            	im1, ...                    % Reference image.
                im2, ...                    % Measurement image.
                resolution, ...             % Resolution, in the Y direction (rows), of the PIV grid.
                resolution, ...             % Resolution, in the X direction (columns), of the PIV grid.
                resolution_z, ...           % Resolution, in the Z direction (layers), of the PIV grid.
                overlap, ...                % Overlap between blocks, in XY.
                overlap_z, ...              % Overlap between blocks, in Z.
                padding, ...                % Padding for the PIV, in XY, usually 0.
                paddingz, ...               % Padding in the Z-direction.
                blocksizedec, ...           % Amount by which each side of the interrogation box shrinks at every iteration.
                blocksizedec_z, ...         % In XY and Z respectively.
                fft_window, ...             % Window for FFT calculation. 1 means hanning.
                corr_method, ...            % Method to compute max of cross-correlation. 1 for 3DCentroid, 2 for Gaussian fit
                pk_window, ...              % Size of the sub-window for the peak of cross-correlation subpixel location.
                iter_max, ...               % Maximum number of iterations after decrease of box size.
                iter_epsilon, ...           % Maximum acceptable relative error of the iteration.
                iter_N_conver, ...          % N_conver: stop after N_conver iterations without improving the convergence.
                longdisp_del, ...           % For each box where longdisp is required, apply longdisp to the neighbouring "del" 
                longdisp_del_z, ...         % boxes around it. In XY and Z respectively.
                dynamic_window, ...         % 1 if the sizes of the PIV windows are determined dynamically; 0 otherwise.
                filt_size, ...              % Size of the Predictor filter, in window lengths, i.e. 
                filt_size_z, ...            % filter_size = filt_size*window_size. In XY and Z respectively.
                filt_pow, ...               % Power of the Predictor filter.
                n_th_disp ...
                );

            dx = inpaint_nans3(dx);
            dy = inpaint_nans3(dy);
            dz = inpaint_nans3(dz);

            % Save displacements:
            f3 = fopen(beads_file, 'w');

            for jj = 1:numel(x)
            	fprintf(f3, '%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n', ...
                        x(jj), y(jj), z(jj), dx(jj), dy(jj), dz(jj), pkh(jj), ws_xy0(jj), ws_z0(jj), Convergence_Mask(jj));
            end

            fclose(f3);

            f3 = fopen([DispPath, filesep, 'DispConvergence_', temp, '.dat'], 'w');

            for jj = 1:numel(conver)
            	fprintf(f3, '%15.5e%15.5e\n', conver(jj), iter_time(jj));
            end

            fclose(f3);

            File.TractionSize(ii).i = size(x,1);
            File.TractionSize(ii).j = size(x,2);
            File.TractionSize(ii).k = size(x,3);  

        end

    end

    save([File.pathname, filesep, 'File.mat'], 'File', '-mat');

    if exist([File.pathname, filesep, 'File_old.mat'], 'file')

        delete([File.pathname, filesep, 'File_old.mat']);

    end

    disp(' ');
    disp('        *** End displacement finder ***        ');
    disp(' ');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
