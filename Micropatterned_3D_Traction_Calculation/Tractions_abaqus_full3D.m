%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                 Function that writes an abaqus Fortran-subroutine to impose displacement Boundary Conditions                 %
%         (coming from PIV analysis) at the top nodes of a 3D Boundary Value Problem to be solved with Abaqus FEM code.        %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       Settings [structure]: Settings structure coming from the 3D PIV.                                                       %
%       File [structure]: File structure coming from the 3D PIV.                                                               %
%       exp_type [str]: String used to identify the type of experiment.                                                        %
%       Settings_abaqus [structure]: Settings structure containing some Abaqus parameters.                                     %
%       ff [int]: position analyzed.                                                                                           %
%       tt [int]: timepoint analyzed.                                                                                          %
%       t_0 [int]: value of the first timepoint of the series.                                                                 %
%       n_th [int]: number of threads used in the parfor loops.                                                                %
%       n_th_abaqus [int]: number of threads used by abaqus.                                                                   %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       File [structure]: Modified File structure.                                                                             %
%                                                                                                                              %
%   Last Revison Date: 22/01/2024                                                                                              %
%   Created by Ernest Latorre-Ibars and Manuel Gomez Gonzalez                                                                  %
%                                                                                                                              %
%   References:                                                                                                                %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function File = Tractions_abaqus_full3D(Settings, File, exp_type, Settings_abaqus, ff, tt, t_0, n_th, n_th_abaqus)

    %==========================================================================================================================%
    %                                                                                                                          %
    % This script is a bit mesy in the sense that there are many, many coordinate systems and projections of vectors. This is  %
    % a small list of the coordinate systems used.                                                                             %
    %                                                                                                                          %
    %   [  X_PIV ,  Y_PIV ,  Z_PIV  ] (Matrices of size [ny, nx, nz]) (units of microns):                                      %
    %           coordinates of the PIV data (center of the PIV windows), with origin at the pillar's center (in x and y) and   %
    %           at its base (in z).                                                                                            %
    %                                                                                                                          %
    %   [ DX_PIV , DY_PIV , DZ_PIV  ] (Matrices of size [ny, nx, nz]) (units of microns):                                      %
    %           displacements obtained with PIV, for PIV boxes centered at [ X_PIV, Y_PIV, Z_PIV ].                            %
    %                                                                                                                          %
    %   [  x_fluo,  y_fluo,  z_fluo ] (Matrices of size [sy, sx, sz]) (units in piexels):                                      %
    %           coordinates of the microscope data stacks (fluorescence, beads, etc.), with origin at the pillar's center      %
    %           (in x and y) and at its base (in z).                                                                           %
    %                                                                                                                          %
    %   [  x_img ,  y_img ,  z_img  ] (Matrices of size [sy, sx, sz]) (units in piexels):                                      %
    %           coordinates of the microscope data stacks (fluorescence, beads, etc.), with origin at the lower left corner    %
    %           (in x and y) and at the bottom (in z).                                                                         %
    %                                                                                                                          %
    %   [  X_nodes_q ,  Y_nodes_q ,  Z_nodes_q ] (Vectors of size [:]) (units in microns):                                     %
    %           coordinates of the surface nodes of the mesh of the FEM model, where there is fluorescence intensity,          %
    %           with origin at the pillar's center (in x and y) and at its base (in z) (as opposed to coords).                 %
    %                                                                                                                          %
    %   [ DX_nodes_q , DY_nodes_q , DZ_nodes_q ] (Vectors of size [:]) (units in microns):                                     %
    %           displacements, interpolated from the PIV data, of the surface nodes of the mesh of the FEM model, where there  %
    %           is fluorescence intensity, located at [ X_nodes_q, Y_nodes_q, Z_nodes_q ].                                     %
    %                                                                                                                          %
    %   [ DX_PIV_q , DY_PIV_q , DZ_PIV_q ] (Matrices of size [sy, sx, sz]) (units in microns):                                 %
    %           displacements, interpolated from the FEM data, of the surface nodes of the mesh of the FEM model, where there  %
    %           is fluorescence intensity, interpolated at [ x_fluo, y_fluo, z_fluo ].                                         %
    %                                                                                                                          %
    %   [ Sxx_nodes  , Syy_nodes  , Szz_nodes  , Sxy_nodes  , Sxz_nodes  , Syz_nodes  , P_nodes ] (Vectors of size [:]):       %
    %           Cauchy stress tensor, given in the deformed configuration, at the nodes, located at                            %
    %           [ x_def_nodes, y_def_nodes, Z_def_nodes ].                                                                     %
    %                                                                                                                          %
    %   [ tzx_nodes, tzy_nodes, tzz_nodes,  tn_nodes, tt_x_nodes, tt_y_nodes, tt_z_nodes, tt_nodes ] (Vectors of size [:]):    %
    %           Tractions at the external surface of the pillar, decomposed in (x, y, z) and in tangential and normal. Given   %
    %           in the deformed configuration, at the nodes [ X_nodes_q, Y_nodes_q, Z_nodes_q ].                               %
    %                                                                                                                          %
    %   [ tzx_PIV_q, tzy_PIV_q, tzz_PIV_q,  tn_PIV_q, tt_x_PIV_q, tt_y_PIV_q, tt_z_PIV_q, tt_PIV_q ]                           %
    %           (Matrices of size [ny, nx, nz]):                                                                               %
    %           Tractions at the external surface of the pillar, decomposed in (x, y, z) and in tangential and normal. Given   %
    %           in the deformed configuration, at the center of the PIV windows [ X_PIV, Y_PIV, Z_PIV ].                       %
    %                                                                                                                          %
    %   [ tzx_fluo_q, tzy_fluo_q, tzz_fluo_q,  tn_fluo_q, tt_x_fluo_q, tt_y_fluo_q, tt_z_fluo_q, tt_fluo_q ]                   %
    %           (Matrices of size [sy, sx, sz]):                                                                               %
    %           Tractions at the external surface of the pillar, decomposed in (x, y, z) and in tangential and normal. Given   %
    %           in the deformed configuration, at the coordinates of the microscope data stack [ x_fluo, y_fluo, z_fluo ].     %
    %                                                                                                                          %
    %   [ X_def_nodes, Y_def_nodes, Z_def_nodes ] (Vectors of size [:]) (units in microns):                                    %
    %           coordinates of the deformed surface nodes of the mesh of the FEM model, where there is fluorescence intensity, %
    %           with origin at the pillar's center (in x and y) and at its base (in z) (as opposed to coords).                 %
    %                                                                                                                          %
    %   coords (Matrix of size [:, 13]): data of the surface, from abaqus, where there is fluorescence intensity.              %
    %           The 1st column is the node number, the 2nd, 3rd and 4th columns are the x, y and z coordinates of each node,   %
    %           in microns ( X_nodes_q, Y_nodes_q, Z_nodes_q-zoffset ). The zero is at the center of the lowest plane.         %
    %           The 5th, 6th and 7th columns contain the x, y and z displacements, in microns                                  %
    %           ( DX_nodes_q, DY_nodes_q, DZ_nodes_q ). The 8th, 9th and 10th columns contain tx, ty and tz.                   %
    %           The 11th, 12th and 13th columns contain the deformed coordinates of the surface (x_def, y_def, z_def-zoffset). %
    %                                                                                                                          %
    %==========================================================================================================================%

    % Choose the figures that will be plotted.
    %   fig1, fig2 and fig3: display the experimental well(beads channel), the mask, the PIV grid, the mesh nodes and the 
    %                        fluorescence of the cell.
    %   fig4: 3D scatter plot of the mesh nodes.
    %   fig5: 3D quiver of the tractions field, in the deformed configuration.
    %   fig6: 3D scatter of the normal tractions, in the deformed configuration.
    %   fig7; 3D scatter of the tangential tractions, in the deformed configuration.
    %   fig8: scatter plot of the mesh nodes, quiver of the normals at the node locations, and quiver of the normals at the PIV 
    %         grid, in the deformed configuration.
    %   fig9: scatter plot of the mesh nodes and quiver of the normals at the node locations, in the deformed configuration.
    %   fig10, fig11, fig12: scatter plot of the mesh nodes, colored by PIV displacements DX, DY and DZ, in microns.
    %   fig13: scatter plot of the element face (at the surface) centers and quiver of the normals.
    %
%   figs_to_plot = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 ] ;
    figs_to_plot = 1:13;

    % Offset in Z used to move the reference system from the Abaqus data to Matlab's data. Abaqus' reference system is located 
    % at the lowest plane, under the center of the model. It has to be moved to the center of the surface plane.
    zoffset = Settings.zoffset;

    % Name of the set of nodes/elements where we impose the displacements (as Boundary Conditions).
    set_name = Settings_abaqus.set_name;

    % Name of the set of nodes/elements at the bottom (encastre).
    set_bottom_name = Settings_abaqus.set_bottom_name;

    % Abaqus' scratch directory;
    scratch_dir = Settings_abaqus.abaqus_scratch_dir;

    % Maximum memory allowed to be used by abaqus, in megabytes.
    mem_use = Settings_abaqus.abaqus_mem_use;

    % Job created in Abaqus to analyze the FEM problem.
    jobabaqus_file = 'Job-111.inp';
    jobabaqus_file_orig = 'Job-111_orig.inp';       % Original file, from Abaqus, before we do our modifications.

    conv_th = 0.80 ;                % Use only the PIV data with a correlation peak value higher than conv_th.

    % Folder and file names used for the analysis.
    pathabaqus = File.pathname;     % Folders where the Abaqus scripts, temporary files and data will reside.

    % Job file that will be sent to Abaqus.
    jobabaqus = [File.AbaqusPath, filesep, jobabaqus_file];

    % Image of the trypsin.
    crop_name = strsplit(File.Name(tt-t_0+1).CropTrypsin, filesep);
    File.Name(tt-t_0+1).CropTrypsin = [File.CroppPath, filesep, crop_name{end}];

    crop_name = strsplit(File.Name(tt-t_0+1).CropFluo_1, filesep);
    File.Name(tt-t_0+1).CropFluo_1 = [File.CroppPath, filesep, crop_name{end}];

    % Size of the PIV displacement grid.
    nx = File.TractionSize(tt-t_0+1).i;
    ny = File.TractionSize(tt-t_0+1).j;
    nz = File.TractionSize(tt-t_0+1).k;

    % Parameters needed for filtering the displacements. Obtain them from the PIV data.
    n_filt_trac = Settings.filt_size_trac;
    power_filt_trac = Settings.filt_pow_trac;
    overlap = Settings.Overlap;

    % If the PIV did not define the filter size and overlap in Z (because it used the same as for XY), define them here.
    if ~isfield(Settings, 'filt_size_z')
        Settings.filt_size_z = Settings.filt_size;
    end

    if ~isfield(Settings, 'Overlap_z')
        Settings.Overlap_z = Settings.Overlap;
    end

    n_filt_trac_z = Settings.filt_size_z;
    overlap_z = Settings.Overlap_z;

    % Padding used for the position and time names.
    pad = Settings.pad;

    % Name of the Abaqus job.
    jobname = [exp_type, '_f', sprintf(['%0', num2str(pad), 'd'], ff), '_t', sprintf(['%0', num2str(pad), 'd'], tt)];

    for suffix = {'.com', '.dat', '.msg', '.prt', '.sim', '.sta', ...
                  '_deformed_coord.dat', '_reactions.dat', '_stress_tensor.dat', '.odb', '.lck'}
        try
            delete([jobname, suffix{1}]);
        catch
        end
    end

    % In order to impose the displacement B.C.s to Abaqus, we need to do it through a Fortran subroutine. Given that we need to
    % impose a list of nodes and displacements, we will dynamically create the Fortran subroutine. It will be stored in the
    % following .f file.

    % Fortran subroutine used used by Abaqus to impose the PIV displacements.
    node_disp_subroutine_file = ['node_disp_subroutine', '_f', sprintf(['%0', num2str(pad), 'd'], ff), ...
                                                         '_t', sprintf(['%0', num2str(pad), 'd'], tt), '.f'];
    node_disp_subroutine = [File.AbaqusPath, filesep, node_disp_subroutine_file] ;

    % Abaqus results will be read from an ODB file. In order to access it, we use a Python script. We need to create this
    % script dynamically. It will be stored in the following .py file.
    
    % Extract the calculated results for the surface nodes/elements
    ReadODB = [File.AbaqusPath, filesep, 'ReadODB', '_f', sprintf(['%0', num2str(pad), 'd'], ff), ...
                                                    '_t', sprintf(['%0', num2str(pad), 'd'], tt), '.py'];
    % Extract the calculated results for all the nodes/elements
    ReadODB_all_elements = [File.AbaqusPath, filesep, 'ReadODB_all_nodes', '_f', sprintf(['%0', num2str(pad), 'd'], ff), ...
                                                                           '_t', sprintf(['%0', num2str(pad), 'd'], tt), '.py'];
    % Import the displacements previously calculated with the PIV.
    TempFile = [File.DispPath, filesep, 'GelDisp_', File.Base.Beads, sprintf(['%0', num2str(pad), 'd'], tt), '.dat'];
    TempData = load(TempFile);

    % Import a list with the coordinates of all the nodes in the mesh, the node number of the nodes that belong to the boundary 
    % where the PIV displecements will be imposed, etc..
    if isfield(Settings, "constrained") && strcmp(Settings.constrained, 'Constrained')

        % Save a copy of the original abaqus job file as jobabaqus_file_orig.
        copyfile(jobabaqus, [File.AbaqusPath, filesep, jobabaqus_file_orig])

        % Check if jobabaqus is using hybrid elements. If so, modify it to non-hybrid elements. Our constrained method uses the 
        % stiffness matrix calculated by Abaqus. It will not work if Abaqus uses its hybrid formulation.
        modify_hybrid_elements([File.AbaqusPath, filesep, jobabaqus_file]);

        % Save a copy the jobabaqus as jobabaqus_file_cell_loads.
        jobabaqus_file_cell_loads = 'Job-111_cell_loads.inp';

        copyfile(jobabaqus, [File.AbaqusPath, filesep, jobabaqus_file_cell_loads])

        % Create and save the job that computes the stiffness matrix for the F.E.M. problem.
        jobabaqus_file_matrix_input = 'Job-111_mtx_matrix_input.inp' ;
        jobabaqus_matrix_input = [File.AbaqusPath, filesep, jobabaqus_file_matrix_input];

        create_stiffness_matrix_job(jobabaqus, jobabaqus_matrix_input);

        % Obtain the node and element information from the abaqus job.
        [nodes_coord, elements_nodes, nodes_to_impose_disp, elements_surface, nodes_fixed] = ...
                extract_abaqus_data(jobabaqus, set_name, set_bottom_name);

    else

        % Obtain the node and element information from the abaqus job.
        [nodes_coord, elements_nodes, nodes_to_impose_disp, elements_surface, ~]  = ...
                extract_abaqus_data( jobabaqus, set_name, set_bottom_name);

    end

    % Add the node numbers to the list of elements belonging to the surface. Add one column per node belonging to the element.
    elements_surface = sort(elements_surface(:));
    elements_surface = elements_surface(~isnan(elements_surface));
    elements_surface = [elements_surface, zeros(size(elements_surface, 1), size(elements_nodes, 2)-1)];

    % Keep only the nodes of the elements that belong to the surface. Add the nodes of each of those elements.
    for jj = 1:size(elements_surface, 1)
        elements_surface(jj, 2:end) = elements_nodes(elements_surface(jj, 1) == elements_nodes(:, 1), 2:end);
    end

    % For each element that is on the surface, one, two or three nodes belong to the surface, the other(s) do(es)n't. 
    % Unless the element is degenerated or in a corner. The following detection and selection of node elements will breake in 
    % that case. We need to check that in the model/mesh we use there are no or very few of those elements. If there are zero 
    % of those elements, there is no problem. If there are a few of them, a small local error is expected to be introduced. 
    % To do: for now it is enough to manually check the model/mesh before running the code. In the future, I need to implement 
    % a better authomatic case handling.

    % Substitute the node numbers that don't belong to the surface by zeros.
    elements_surface_all = elements_surface;        % Vector containing all the nodes of the element.
    elements_surface(:, 2:end) = ...
            elements_surface(:, 2:end).*ismember(elements_surface(:, 2:end), nodes_to_impose_disp);

    % Keep only the elements with three (or more) nodes on the surface, i.e. discard those with only one or two nodes on the 
    % surface.
    elements_surface_all = elements_surface_all(sum(elements_surface(:, 2:end)~=0, 2) > 2, :);
    elements_surface = elements_surface(sum(elements_surface(:, 2:end)~=0, 2) > 2, :);

    % Create a matrix for the nodes where the PIV displacements will be imposed, "coords". They are on the surface.
    %   coords(:, 1): node number.
    %   coords(:, 2:4): x, y and z coordinates of each node, in microns. The zero is at the center of the lowest plane.
    %   coords(:, 5:7): x, y and z displacements, in microns.
    %   coords(:, 8:10): x, y and z tractions, in Pa.
    %   coords(:, 11:13): deformed x, y and z coordinates, in microns.
    coords = NaN(numel(nodes_to_impose_disp), 13);
    coords(:, 1:size(nodes_coord,2)) = nodes_coord(nodes_to_impose_disp, :);

    % Importing 3D discrete displacement field from PIV analysis and scaling it from pixels to um.
    % Colums 1, 2 and 3 are point coordinates and colums 4, 5 and 6 are Dx, Dy and Dz respectively.
    % Coordinates start in the lowest plane, lower left corner.
    PIV_data = TempData(:, 1:6);
    PIV_data(:, [1:2, 4:5]) = PIV_data(:, [1:2, 4:5])*Settings.PixelSizeXY;
    PIV_data(:, [3, 6]) = PIV_data(:, [3, 6])*Settings.PixelSizeZ;

    % If the image stack was acquired inverted, meaning the bottom part of the gel is at the top part of the stack, the PIV 
    % data has to be inverted in Z. It is controlled by the parameter Settings_abaqus.dir_inv)='Inverted'.
    if isfield(Settings_abaqus, 'dir_inv') && strcmp(Settings_abaqus.dir_inv, 'Inverted')

        PIV_data(:, 3) =  PIV_data(end:-1:1, 3);
        PIV_data(:, 6) = -PIV_data(      : , 6);

    end

    % Center the PIV coordinate grid. Move the origin to the center of the lowest plane.
    X_PIV = PIV_data(:, 1) - (max(PIV_data(:, 1)) + min(PIV_data(:, 1)))/2;
    Y_PIV = PIV_data(:, 2) - (max(PIV_data(:, 2)) + min(PIV_data(:, 2)))/2;
    Z_PIV = PIV_data(:, 3);

    X_PIV = reshape(X_PIV, [nx, ny, nz]);
    Y_PIV = reshape(Y_PIV, [nx, ny, nz]);
    Z_PIV = reshape(Z_PIV, [nx, ny, nz]);

    % Find the location of the top plane of the gel.
    % Load the trypsin stack, find the median of the intensity of each plane, locate the maximum of this intensity profile.
    im3 = double(TIFFStack(File.Name(tt-t_0+1).CropTrypsin));
    im_fluo = double(TIFFStack(File.Name(tt-t_0+1).CropFluo_1));

    if isfield(Settings_abaqus, 'dir_inv') && strcmp(Settings_abaqus.dir_inv, 'Inverted')

        im3 = im3(:, :, end:-1:1);
        im_fluo = im_fluo(:, :, end:-1:1);

    end

    [sy, sx, sz] = size(im3);

    % Find the location of the upper plane of the gel, max_int, and the central plane of the well, cent_int. Multiple cases can 
    % be controlled with two parameters:
    %   Settings.auto_surf_finding: if =1 or if it is not defined, the surface will be found automatically. 
    %                               If =0, the surface plane is selected manually.
    %   Settings.beads_on_gel_surface: if =1 or if it is not defined, assume that the beads are on the surface of the gel. 
    %                                  If =0, assume that the beads are everyehere in the gel.
    %
    % We distinguish three cases to find max_int:
    %   1) If the gel has beads only on its surface, we calculate the average intensity of each plane, and the plane with the 
    %      higher intensity, max_int is the surface.
    %   2) If the gel has beads everywhere, we calculate the average intensity of each plane, and then the peaks of its profile. 
    %      The upper surface, max_int, is the last peak, which might be lower than the previous peaks.
    %   3) Manual plane selection if neither of those work.
    %
    % We distinguish two cases to find cent_int:
    %   1) If the surface of the gel was selected manually, or if the beads are everywhere inside the gel, the center is 
    %      calculated with the depth of the well, because we don't have a reliable way to find the bottom of the gel. We could 
    %      ask the user for the bottom plane, but the shape of the well should not vary much, so the detection is not crucial.
    %   2) If the beads are at the surface of the gel, it is not possible to find the bottom plane of the well by looking at 
    %      the intensity profile, its peak is too shallow and wide. It can be reliably found, howerver, with the profile of the 
    %      standard deviation of the intensity at each plane. The accuracy is lower for the top plane of a PAA gel, but much 
    %      higher for the bottom plane of the well.

%     z_profile = squeeze(median(im3, [1, 2]));
    z_profile = squeeze(mean(im3, [1, 2]));     % Average intensity of each plane of the trypsin stack.

    if ~isfield(Settings, 'auto_surf_finding')||Settings.auto_surf_finding

        if ~isfield(Settings, 'beads_on_gel_surface')||Settings.beads_on_gel_surface

            [~, max_int] = max(z_profile);      % max_int is the plane where the surface is located.

            std_inte = zeros(1, sz);

            % Calculate the std intensity of each plane.
            for jj=1:sz
                std_inte(jj) = std(im3(:, :, jj), 0, 'all');
            end

            [pks, locs] = findpeaks(std_inte);
            [~, ind] = maxk(pks, 2);

            cent_int = round(sum(locs(ind))/2); % cent_int is the central plane of the well.

        else

            [~, locs_mean] = findpeaks(z_profile);
            max_int = locs_mean(end);

            cent_int = round(locs_mean(end) - Settings.h/Settings.PixelSizeZ/2);

        end

    else

        max_int = manually_select_plane(File, Settings, im3);

        cent_int = round(max_int - Settings.h/Settings.PixelSizeZ/2);

    end

    % Find the center and radius of the well. First find a binary mask of the central plane, then fit it to a circle.
    im = im3(:, :, cent_int);
    im = imadjust(mat2gray(im), [0 0.2], []);

    % Create a polygonal mask of the pillar base.
    mask = mask_beads_manual(double(im));

    % Find the center and mean radius of the mask, i.e. the pillar base.
    [Index_row_p, Index_col_p] = find(bwperim(mask));

    [z, R_mean] = fitcircle([Index_row_p, Index_col_p]');

    cylinder.x_cent = z(2);
    cylinder.y_cent = z(1);
    cylinder.R = R_mean;

    TempFile = [File.BoundaryPath, filesep, 'mask_well_shape.mat'];
    save(TempFile, 'cylinder');

    % Shift the center of the coordinate system to the center of the well in (x, y) and the upper plane of the gel.
    x_cent = z(2);
    y_cent = z(1);
    disp(['Radius of the experimental well: ', num2str(R_mean*Settings.PixelSizeXY), ' microns.']);

    X_PIV = X_PIV - (x_cent - sx/2)*Settings.PixelSizeXY;
    Y_PIV = Y_PIV - (y_cent - sy/2)*Settings.PixelSizeXY;
    Z_PIV = Z_PIV - max_int*Settings.PixelSizeZ;

    x_fluo = (1:sx) - x_cent;
    y_fluo = (1:sy) - y_cent;
    z_fluo = (1:sz) - max_int;

    [y_fluo, x_fluo, z_fluo] = ndgrid(y_fluo, x_fluo, z_fluo);
    [y_img , x_img , z_img ] = ndgrid(1:sy, 1:sx, 1:sz);

    % Create a 3D stack mask of the ideal well. Use the center of the experimental well and the radius of the mesh.
    [~, R_FEM] = fitcircle([coords((coords(:, 4)<=-zoffset-1)&(coords(:, 4)>=-zoffset-Settings.h+1), 2), ...
                            coords((coords(:, 4)<=-zoffset-1)&(coords(:, 4)>=-zoffset-Settings.h+1), 3)]');
    disp(['Radius of the model well: ', num2str(R_FEM), ' microns.']);

    mask_circ = zeros(sy, sx);      % Mask of the ideal circular well.
    mask_circ(x_fluo(:, :, 1).^2 + y_fluo(:, :, 1).^2 <= (R_FEM/Settings.PixelSizeXY).^2) = 1;

    [Index_row_circ, Index_col_circ] = find(bwperim(mask_circ));        % Circumference of the ideal well.

    % Mask in the image grid.
    mask_well_fluo = zeros(sy, sx, sz);

    for jj=1:numel(Index_row_circ)
        mask_well_fluo(Index_row_circ(jj), Index_col_circ(jj), 1:max_int) = 1;
    end

    % Save the mask stack.
    mask_fluo_name = 'Mask_fluo_Well.tif';

    imwrite(uint16(mask_well_fluo(:, :, 1)), [File.MaskPath, filesep, mask_fluo_name], 'compression', 'none');
    for kk=2:sz
        imwrite(uint16(mask_well_fluo(:, :, kk)), [File.MaskPath, filesep, mask_fluo_name], 'WriteMode', 'append', ...
                'compression', 'none')
    end

    % Same mask in the PIV grid.
    mask_well_PIV = zeros(ny, nx, nz);

    dx = X_PIV( 1, 2, 1 ) - X_PIV( 1, 1, 1 ) ;
    dy = Y_PIV( 2, 1, 1 ) - Y_PIV( 1, 1, 1 ) ;
    dz = Z_PIV( 1, 1, 2 ) - Z_PIV( 1, 1, 1 ) ;

    % This is a little bit tricky. When comparing both grids, there might be a very small offset, due to 
    % truncation errors. Add a very small offset to the limits.
    eps_grid = 1e-12 ;

    % Go through each point in the PIV grid. Define a voxel around it. Find the points of the image grid that lie inside this 
    % voxel. Check if any of those points have a value of 1 in the masks stack. If so, assign a value of 1 in the PIV grid 
    % mask. The result is slightly different than creating the mask with the same method as mask_well_fluo but with the coarser 
    % PIV grid.
    for ll=1:nx*ny*nz
%     parfor( ll=1:nx*ny*nz, n_th )

        % Matlab complains about the use of the find() function here, and recommends using logical indexing. However, my tests 
        % show that, for this specific loop, logical indexing is 50% slower.
        ind_x = find((squeeze(x_fluo(1, :, 1))*Settings.PixelSizeXY>=X_PIV(ll)-dx/2-eps_grid)& ...
                     (squeeze(x_fluo(1, :, 1))*Settings.PixelSizeXY<=X_PIV(ll)+dx/2+eps_grid));
        ind_y = find((squeeze(y_fluo(:, 1, 1))*Settings.PixelSizeXY>=Y_PIV(ll)-dy/2-eps_grid)& ...
                     (squeeze(y_fluo(:, 1, 1))*Settings.PixelSizeXY<=Y_PIV(ll)+dy/2+eps_grid));
        ind_z = find((squeeze(z_fluo(1, 1, :))*Settings.PixelSizeZ >=Z_PIV(ll)-dz/2-eps_grid)& ...
                     (squeeze(z_fluo(1, 1, :))*Settings.PixelSizeZ <=Z_PIV(ll)+dz/2+eps_grid));

        if any(mask_well_fluo(ind_y, ind_x, ind_z), 'all')
            mask_well_PIV(ll) = 1;
        end

    end

%     delete(gcp('nocreate'));

    % Save the PIV mask as a text file
    f = fopen([File.MaskPath, filesep, 'Mask_PIV_Well.dat'], 'w');

    ind = find(mask_well_PIV == 1);

    for jj = 1:numel(ind)
        fprintf(f, '%15.5e%15.5e%15.5e%15.5e\n', X_PIV(ind(jj)), Y_PIV(ind(jj)), Z_PIV(ind(jj)), mask_well_PIV(ind(jj)));
    end

    fclose(f) ;

    % Keep only the displacements at the points where the fluorescence signal is higher than threshold_fluo.
    im_fluo_filt = imgaussfilt(im_fluo, 2.0);

    im_fluo_filt2 = NaN(size(im_fluo_filt));
    mask2 = double(mask) ;
    mask2(~mask2) = NaN;

    for kk=1:sz
        im_fluo_filt2(:, :, kk) = im_fluo_filt(:, :, kk).*mask2;
    end

%     threshold_fluo = nanmean(im_fluo_filt(:));
    threshold_fluo = nanmean(im_fluo_filt2(:));
    mask_fluo_0 = (im_fluo_filt >= threshold_fluo);

    % Dilate the mask.
    se = strel('disk', 11);
    mask_fluo = imdilate(mask_fluo_0, se);

    % Save the mask stack.
    imwrite(uint16(mask_fluo(:, :, 1)), [File.MaskPath, filesep, 'Mask_fluo_CAFs.tif'], 'compression', 'none')
    for kk=2:sz
        imwrite(uint16(mask_fluo(:, :, kk)), [File.MaskPath, filesep, 'Mask_fluo_CAFs.tif'], ...
                'WriteMode', 'append', 'compression', 'none')
    end

    % Calculate the intersection between the well and the fluorescence masks.
    se = strel('disk', 2);
    mask_cell_on_well = imdilate(mask_well_fluo, se).*mask_fluo;

    % Save the mask stack.
    imwrite(uint16(mask_cell_on_well(:, :, 1)), [File.MaskPath, filesep, 'Mask_fluo_cell_on_well.tif'], 'compression', 'none')

    for kk=2:sz
        imwrite(uint16(mask_cell_on_well(:, :, kk)), [File.MaskPath, filesep, 'Mask_fluo_cell_on_well.tif'], ...
                'WriteMode', 'append', 'compression', 'none')
    end

    mask_cell_on_well_filled = zeros(size(mask_well_fluo));
    mask_cell_on_well_filled(:, :, max_int) = 1-mask_circ;

    pillar_pixels = find(mask_well_fluo);
    [row_pillar, col_pillar, h_pillar] = ind2sub(size(mask_well_fluo), pillar_pixels);

    for mm=1:numel(pillar_pixels)

        if any(mask_cell_on_well(row_pillar(mm), col_pillar(mm), h_pillar(mm):end))
            mask_cell_on_well_filled(row_pillar(mm), col_pillar(mm), h_pillar(mm)) = 1;
        end

    end

    % Save the mask stack.
    imwrite(uint16(mask_cell_on_well_filled(:, :, 1)), ...
            [File.MaskPath, filesep, 'Mask_fluo_cell_on_well_filled.tif'], 'compression', 'none')
    for kk=2:sz
        imwrite(uint16(mask_cell_on_well_filled(:, :, kk)), ...
                [File.MaskPath, filesep, 'Mask_fluo_cell_on_well_filled.tif'], 'WriteMode', 'append', 'compression', 'none')
    end

    % x and y coordinates of the edges of the PIV data.
    a = X_PIV(1,1,1)/Settings.PixelSizeXY;
    b = X_PIV(1,end,1)/Settings.PixelSizeXY;
    c = Y_PIV(1,1,1)/Settings.PixelSizeXY;
    d = Y_PIV(end,1,1)/Settings.PixelSizeXY;

    if any(figs_to_plot == 1)

        % Plot the plane of the beads at the surface of the gel; mark the PIV coordinates with circles; 
        % mark the abaqus nodes with circles; and mark the edge of the well mask.
        fig1 = figure(1);

        set(fig1, 'Visible', 'on');
        ax1 = axes('Parent', fig1);

        imagesc(ax1, X_PIV(1, [1, end], 1)/Settings.PixelSizeXY, Y_PIV([1, end], 1, 1)/Settings.PixelSizeXY, ...
                imadjust(mat2gray(im), [0 0.2], []));
        colormap(ax1, 'gray');
        hold(ax1, 'on');
        scatter(ax1, X_PIV(:)/Settings.PixelSizeXY, Y_PIV(:)/Settings.PixelSizeXY, 'DisplayName','PIV grid');
        axis(ax1, 'image');

        scatter(ax1, ...
                coords(coords(:, 4)<=-zoffset+1, 2)/Settings.PixelSizeXY, ...
                coords(coords(:, 4)<=-zoffset+1, 3)/Settings.PixelSizeXY, 'DisplayName','Mesh nodes');

        scatter(ax1, ((b-a)/(sx-1)*(Index_col_p-sx)+b), ((d-c)/(sy-1)*(Index_row_p-sy)+d), 'DisplayName','Mask');

        xlabel(ax1, 'X [px]');
        ylabel(ax1, 'Y [px]');

        legend(ax1, 'location', 'bestoutside');

    end

    if any(figs_to_plot == 2)

        % Plot a yz cut of the beads through the center of the well; and mark the abaqus nodes with circles.
        fig2 = figure(2);

        set(fig2, 'Visible', 'on');
        ax2 = axes('Parent', fig2);

        imagesc(ax2, ([1, sz] - max_int)*Settings.PixelSizeZ-zoffset, [min(Y_PIV(:)), max(Y_PIV(:))], ...
                imadjust(mat2gray(squeeze(im3(:, round(x_cent), :))), [0 0.2], []));
        colormap(ax2, 'gray');
        hold(ax2, 'on');
        scatter(ax2, coords(:, 4), coords(:, 3), 'DisplayName','Mesh nodes');
        axis(ax2, 'image');

        xlabel(ax2, 'Z [\mum]');
        ylabel(ax2, 'Y [\mum]');

        legend(ax2, 'location', 'bestoutside');

    end

    if any(figs_to_plot == 3)

        % Plot a yz projection of the fluorescence; and mark the abaqus nodes with circles.
        fig3 = figure(3);

        set(fig3, 'Visible', 'on');
        ax3 = axes('Parent', fig3);

        imagesc(ax3, ([1, sz] - max_int)*Settings.PixelSizeZ-zoffset, [min(Y_PIV(:)), max(Y_PIV(:))], ...
                imadjust(mat2gray(squeeze(sum(im_fluo, 2))), [0 0.2], []));
        colormap(ax3, 'gray');
        hold(ax3, 'on');
        axis(ax3, 'image');
        scatter(ax3, coords(:, 4), coords(:, 3), 'DisplayName','Mesh nodes');

        xlabel(ax3, 'Z [\mum]');
        ylabel(ax3, 'Y [\mum]');

        legend(ax3, 'location', 'bestoutside');

    end

    clear("im3");

    % Reshaping and filter DX_PIV,DY_PIV and DZ_PIV.
    DX_PIV = reshape(PIV_data(:,4), [nx ny nz]);
    DY_PIV = reshape(PIV_data(:,5), [nx ny nz]);
    DZ_PIV = reshape(PIV_data(:,6), [nx ny nz]);

    % Keep only the PIV points which have converged and that have a high correlation value.
    pkh = reshape(TempData(:,7), [nx ny nz]);
    conv_mask = reshape(TempData(:,end), [nx ny nz]);
    mask = (conv_mask)&(pkh >= conv_th);

%     mask(:) = 1;        % Use this line to test how the convergence criterion performs.

    DX_PIV = removeOutliers(DX_PIV);
    DY_PIV = removeOutliers(DY_PIV);
    DZ_PIV = removeOutliers(DZ_PIV);

    DX_PIV = DX_PIV{1};
    DY_PIV = DY_PIV{1};
    DZ_PIV = DZ_PIV{1};

    DX_PIV = filter_3D_Disp(DX_PIV, ...
                            [n_filt_trac/(1 - overlap), n_filt_trac/(1 - overlap), n_filt_trac_z/(1 - overlap_z)], ...
                            power_filt_trac);
    DY_PIV = filter_3D_Disp(DY_PIV, ...
                            [n_filt_trac/(1 - overlap), n_filt_trac/(1 - overlap), n_filt_trac_z/(1 - overlap_z)], ...
                            power_filt_trac);
    DZ_PIV = filter_3D_Disp(DZ_PIV, ...
                            [n_filt_trac/(1 - overlap), n_filt_trac/(1 - overlap), n_filt_trac_z/(1 - overlap_z)], ...
                            power_filt_trac);

    % If the FOV is in equilibrium, the mean of the displacements can be substracted here.
%     % Substracting the mean (is the FOV in equilibrium?).
%     DX_PIV = DX_PIV - nanmean(DX_PIV(:));
%     DY_PIV = DY_PIV - nanmean(DY_PIV(:));
%     DZ_PIV = DZ_PIV - nanmean(DZ_PIV(:));

    % The PIV grid and the FEM mesh will, in general, be very different. Keep only the mesh nodes inside of the PIV grid 
    % x and y limits.
    mask2 = (coords(:,2)>=max(min(X_PIV(:)), min(X_PIV(:))-(File.Drift.shx_centroid+1)*Settings.PixelSizeXY)) & ...
            (coords(:,2)<=min(max(X_PIV(:)), max(X_PIV(:))-(File.Drift.shx_centroid+1)*Settings.PixelSizeXY)) & ...
            (coords(:,3)>=max(min(Y_PIV(:)), min(Y_PIV(:))-(File.Drift.shy_centroid+1)*Settings.PixelSizeXY)) & ...
            (coords(:,3)<=min(max(Y_PIV(:)), max(Y_PIV(:))-(File.Drift.shy_centroid+1)*Settings.PixelSizeXY));
    coords = coords(mask2, :);

    % Will interpolate the PIV displacement in the FEM mesh grid nodes.
    X_nodes_q = coords(:,2);
    Y_nodes_q = coords(:,3);
    Z_nodes_q = coords(:,4) + zoffset;

    if any( figs_to_plot == 4 )

        % Mark in 3D the abaqus nodes with circles.
        fig4 = figure(4);

        set(fig4, 'Visible', 'on');
        ax4 = axes('Parent', fig4);

        scatter3(ax4, X_nodes_q, Y_nodes_q, Z_nodes_q, 'DisplayName','Mesh nodes');
        axis(ax4, 'image');

        xlabel(ax4, 'X [\mum]');
        ylabel(ax4, 'Y [\mum]');
        zlabel(ax4, 'Z [\mum]');

        legend(ax4, 'location', 'northeast');

    end

    % Interpolate the PIV data fields DX_PIV, DY_PIV and DZ_PIV, which are applied on (X_PIV, Y_PIV, Z_PIV),
    % to the query points (X_nodes_q, Y_nodes_q, Z_nodes_q).
    FDX_PIV = scatteredInterpolant(X_PIV(mask), Y_PIV(mask), Z_PIV(mask), DX_PIV(mask), 'natural');
    DX_nodes_q = FDX_PIV(X_nodes_q, Y_nodes_q, Z_nodes_q);
    FDY_PIV = scatteredInterpolant(X_PIV(mask), Y_PIV(mask), Z_PIV(mask), DY_PIV(mask), 'natural');
    DY_nodes_q = FDY_PIV(X_nodes_q, Y_nodes_q, Z_nodes_q);
    FDZ_PIV = scatteredInterpolant(X_PIV(mask), Y_PIV(mask), Z_PIV(mask), DZ_PIV(mask), 'natural');
    DZ_nodes_q = FDZ_PIV(X_nodes_q, Y_nodes_q, Z_nodes_q);

    % Apply the RI mismatch correction.
    if isfield(Settings, "Correct_Disp_RI") && Settings.Correct_Disp_RI

        [DX_nodes_q_corrected, DY_nodes_q_corrected, DZ_nodes_q_corrected] = ...
                        correct_disp_RI_matched(DX_nodes_q, DY_nodes_q, DZ_nodes_q, X_nodes_q, Y_nodes_q, Z_nodes_q, Settings);

        DX_nodes_q = DX_nodes_q_corrected;
        DY_nodes_q = DY_nodes_q_corrected;
        DZ_nodes_q = DZ_nodes_q_corrected;

    end

    if any(figs_to_plot == 10)

        % Mark in 3D the abaqus nodes with circles, colored by displacements.
        fig10 = figure(10);

        set(fig10, 'Visible', 'on');
        ax10 = axes('Parent', fig10);

        scatter3(ax10, X_nodes_q, Y_nodes_q, Z_nodes_q, [], DX_nodes_q, 'filled');
        axis(ax10, 'image');
        cb10 = colorbar(ax10);

        xlabel(ax10, 'X [\mum]');
        ylabel(ax10, 'Y [\mum]');
        zlabel(ax10, 'Z [\mum]');

        ylabel(cb10, 'PIV u_x [\mum]');

    end

    if any(figs_to_plot == 11)

        % Mark in 3D the abaqus nodes with circles, colored by displacements.
        fig11 = figure(11);

        set(fig11, 'Visible', 'on');
        ax11 = axes('Parent', fig11);

        scatter3(ax11, X_nodes_q, Y_nodes_q, Z_nodes_q, [], DY_nodes_q, 'filled');
        axis(ax11, 'image');
        cb11 = colorbar(ax11);

        xlabel(ax11, 'X [\mum]');
        ylabel(ax11, 'Y [\mum]');
        zlabel(ax11, 'Z [\mum]');

        ylabel(cb11, 'PIV u_y [\mum]');

    end

    if any(figs_to_plot == 12)

        % Mark in 3D the abaqus nodes with circles, colored by displacements.
        fig12 = figure(12);

        set(fig12, 'Visible', 'on');
        ax12 = axes('Parent', fig12);

        scatter3(ax12, X_nodes_q, Y_nodes_q, Z_nodes_q, [], DZ_nodes_q, 'filled');
        axis(ax12, 'image');
        cb12 = colorbar(ax12);

        xlabel(ax12, 'X [\mum]');
        ylabel(ax12, 'Y [\mum]');
        zlabel(ax12, 'Z [\mum]');

        ylabel(cb12, 'PIV u_z [\mum]');

    end

    % Coords contain, now, in the first column, the node number; in the second, third and fourth columns, 
    % the x, y and z coordinates of each node, in microns; and in the fifth, sixth and seventh columns, the x, y and z 
    % displacements, in microns.
    coords(:, 5) = DX_nodes_q;
    coords(:, 6) = DY_nodes_q;
    coords(:, 7) = DZ_nodes_q;

    % Write the Fortran subroutine used to impose the PIV calculated displacements.
    f = fopen(node_disp_subroutine , 'w');
    fprintf(f, '      SUBROUTINE  DISP(U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)\nC\n      INCLUDE ''ABA_PARAM.INC''\nC\n      DIMENSION U(3),TIME(2),COORDS(3)\nC\n      integer :: i\n      real, dimension(%d) :: n\n      real, dimension(%d) :: dx\n      real, dimension(%d) :: dy\n      real, dimension(%d) :: dz\n', length(coords), length(coords), length(coords), length(coords));
    for k = 1:length(coords)

        node_str = num2str(coords(k, 1));
        fprintf(f, '      n(%d)=', k);
        fprintf(f, node_str);
        fprintf(f, '\n');

        node_str = num2str(coords(k, 5));
        fprintf(f, '      dx(%d)=', k);
        fprintf(f, node_str);
        fprintf(f, '\n');

        node_str = num2str(coords(k, 6));
        fprintf(f, '      dy(%d)=', k);
        fprintf(f, node_str);
        fprintf(f, '\n');

        node_str = num2str(coords(k, 7));
        fprintf(f, '      dz(%d)=', k);
        fprintf(f, node_str);
        fprintf(f, '\n');

    end

    fprintf(f, '      DO i=1,%d,1', length(coords));
    fprintf(f, '\n      IF(NODE.EQ.INT(n(i)).AND.JDOF.EQ.1) THEN\n         U(1)=dx(i)\n      ELSEIF(NODE.EQ.INT(n(i)).AND.JDOF.EQ.2) THEN\n         U(1)=dy(i)\n      ELSEIF(NODE.EQ.INT(n(i)).AND.JDOF.EQ.3) THEN\n         U(1)=dz(i)\n      ENDIF\n      ENDDO\n      RETURN\n      END');
    fclose(f);

    % Results from Abaqus are stored in a ODB file. The Reaction Forces are extracted from this file with a python script. 
    % Here, we dynamically create this script and launch it.   
    f = fopen(ReadODB, 'w');

    fprintf(f, 'import odbAccess\nnames=[''');
    fprintf(f, jobname);
    fprintf(f, ''']\nnameOfStep=''Step-1''\nnodeNumber=[');

    for k = 1:length(coords)-1
        node_str = num2str(coords(k, 1));
        fprintf(f, node_str);
        fprintf(f, ',');
    end

    node_str = num2str(coords(length(coords), 1));
    fprintf(f, node_str);
    fprintf(f, ']\n');
    fprintf(f, 'for y in range(len(names)):\n	NameOfFileR=names[y]+''_reactions.dat''\n	NameOfFileC=names[y]+''_deformed_coord.dat''\n	NameOfFileS=names[y]+''_stress_tensor.dat''\n	FileResultsR=open(NameOfFileR,''w'')\n	FileResultsC=open(NameOfFileC,''w'')\n	FileResultsS=open(NameOfFileS,''w'')\n	Name=names[y]+''.odb''\n	myOdb = odbAccess.openOdb(path=Name)\n	lastStep=myOdb.steps[nameOfStep]\n	for z in range(1, len(lastStep.frames)):\n		lastFrame = myOdb.steps[nameOfStep].frames[z]\n		Times=lastFrame.frameValue\n		ReactionForce = lastFrame.fieldOutputs[''RF'']\n		Coords = lastFrame.fieldOutputs[''COORD'']\n		Stress = lastFrame.fieldOutputs[''S'']\n		for v in ReactionForce.values:\n			if v.nodeLabel in nodeNumber:\n				FileResultsR.write(''%%d\\t %%10.8E\\t %%10.8E\\t %%10.8E\\t \\n'' %% (v.nodeLabel,v.data[0],v.data[1],v.data[2]))\n		for v in Coords.values:\n			if v.nodeLabel in nodeNumber:\n				FileResultsC.write(''%%d\\t %%10.8E\\t %%10.8E\\t %%10.8E\\t \\n'' %% (v.nodeLabel,v.data[0],v.data[1],v.data[2]))\n		for v in Stress.values:\n			if v.nodeLabel in nodeNumber:\n				FileResultsS.write(''%%d\\t %%10.8E\\t %%10.8E\\t %%10.8E\\t %%10.8E\\t %%10.8E\\t %%10.8E\\t %%10.8E\\t \\n'' %% (v.nodeLabel,v.data[0],v.data[1],v.data[2],v.data[3],v.data[4],v.data[5],v.press))\n	myOdb.close()\n	FileResultsR.close()\n	FileResultsC.close()\n	FileResultsS.close()');
    fclose(f);

    % Do the same to obtain the deformed coordinates of all the nodes, rather than only the surface ones.   
    f = fopen(ReadODB_all_elements, 'w');

    fprintf(f, 'import odbAccess\nnames=[''');
    fprintf(f, jobname);
    fprintf(f, ''']\nnameOfStep=''Step-1''\nnodeNumber=[');

    for k = 1:length(nodes_coord)-1
        node_str = num2str(nodes_coord(k, 1));
        fprintf(f, node_str);
        fprintf(f, ',');
    end

    node_str = num2str(nodes_coord(length(nodes_coord), 1));
    fprintf(f, node_str);
    fprintf(f, ']\n');
    fprintf(f, 'for y in range(len(names)):\n	NameOfFileC=names[y]+''_deformed_coord_all_elements.dat''\n	FileResultsC=open(NameOfFileC,''w'')\n	Name=names[y]+''.odb''\n	myOdb = odbAccess.openOdb(path=Name)\n	lastStep=myOdb.steps[nameOfStep]\n	for z in range(1, len(lastStep.frames)):\n		lastFrame = myOdb.steps[nameOfStep].frames[z]\n		Times=lastFrame.frameValue\n		Coords = lastFrame.fieldOutputs[''COORD'']\n		for v in Coords.values:\n			if v.nodeLabel in nodeNumber:\n				FileResultsC.write(''%%d\\t %%10.8E\\t %%10.8E\\t %%10.8E\\t \\n'' %% (v.nodeLabel,v.data[0],v.data[1],v.data[2]))\n		myOdb.close()\n	FileResultsC.close()');
    fclose(f);

    % Calculate the FEM solution in Abaqus.
    
    % There is a bug in Abaqus for Linux where "standard" would not stop when the calculation is completed. The script 
    % stop_abaqus.m checks if the simulation is completed and kills "standard". With the following line, stop_abaqus.m is run in 
    % the background while the rest of the matlab script is run.
    delete(gcp('nocreate'));
    job_stop = batch(@stop_abaqus, 0, {exp_type});

    % If using the constrained method, obtain the striffness matrix.
    if isfield(Settings, "constrained") && strcmp(Settings.constrained, 'Constrained')

        system([Settings_abaqus.abaqus_path, ' job="', jobname, '"', ' input="', jobabaqus_matrix_input, '"', ...
                ' scratch="', scratch_dir, '"', ' memory="', mem_use, ' mb"', ' interactive', ' cpus=', num2str(n_th_abaqus)])

        movefile([jobname, '_STIF2.mtx'], [pathabaqus, filesep, 'abaqusdata', filesep, jobname, '_STIF2_mtx_matrix_input.mtx']);

        for suffix = {'_X2.sim', '.1.SMABulk', '.com', '.dat', '.msg', '.odb', '.prt', '.sim', '.sta'}
            try
                delete([jobname, suffix{1}]);
            catch
            end
        end

    end

    % Launch abaqus analysis.
    system([Settings_abaqus.abaqus_path, ' job="', jobname, '"', ' input="', jobabaqus, '"', ...
            ' scratch="', scratch_dir, '"', ' memory="', mem_use, ' mb"', ' user="', node_disp_subroutine, '"', ...
            ' interactive', ' cpus=', num2str(n_th_abaqus)])

    cancel(job_stop);
    delete(job_stop);
    clear job_stop;

    modify_workers_local_cluster(n_th, n_th);

    % Extract the results from the .odb file to text files.
    system([Settings_abaqus.abaqus_path, ' cae nogui="' ReadODB, '"'])
    system([Settings_abaqus.abaqus_path, ' cae nogui="' ReadODB_all_elements, '"'])

    if isfield(Settings, "constrained") && strcmp(Settings.constrained, 'Constrained')

        % Solve the constrained problem.
        calculate_constrained_tractions(File.AbaqusPath, jobname, coords(:, 1), [DX_nodes_q, DY_nodes_q, DZ_nodes_q], ...
                                        nodes_to_impose_disp, nodes_fixed, nodes_coord, elements_surface, elements_surface_all);

        for suffix = {'_deformed_coord.dat', '_deformed_coord_all_elements.dat', '_reactions.dat', '_stress_tensor.dat', '.odb'}
            suffix_2 = suffix{1};
            suffix_2 = [suffix_2(1:end-4), '_unconstrained', suffix_2(end-3:end)];
            movefile([jobname, suffix{1}], [File.AbaqusPath, filesep, jobname, suffix_2]);
        end

    else

        for suffix = {'_deformed_coord.dat', '_deformed_coord_all_elements.dat', '_reactions.dat', '_stress_tensor.dat', '.odb'}
            movefile([jobname, suffix{1}], File.AbaqusPath);
        end

    end

    for suffix = {'.com', '.dat', '.msg', '.prt', '.sim', '.sta'}
        try
            delete([jobname, suffix{1}]);
        catch
        end
    end

   % Load the tractions, and keep only the tractions calculated in the nodes of the surface. Add tx, ty and tz to 
   % the columns 8, 9 and 10 of coords.
   tractions = importdata([File.AbaqusPath, filesep, jobname, '_reactions.dat']);
   coords(coords(:, 1) == tractions(:, 1), 8:10) = tractions(tractions(:, 1) == coords(:, 1), 2:4);

   % Load the deformed coordinates, and keep only the values calculated at the nodes of the surface. 
   % Add x_def, y_def and z_def to the columns 11, 12 and 13 of coords.
   deformed_coord = importdata([File.AbaqusPath, filesep, jobname, '_deformed_coord.dat']);
   coords(coords(:, 1) == deformed_coord(:, 1), 11:13) = deformed_coord(deformed_coord(:, 1) == coords(:, 1), 2:4);

   % Load the deformed coordinates of the nodes of all the elements. They are used to calculate the normals.
   deformed_coord_all_elements = importdata([File.AbaqusPath, filesep, jobname, '_deformed_coord_all_elements.dat']);

   if ~isfield(Settings, "constrained") || ~strcmp(Settings.constrained, 'Constrained')

        % Load the components of the stress tensor at each point. The stress tensor is calculated at each element's face, 
        % not at each node. So each node has several values of the stress tensor. We average those values in each node.
        stress_tensor_components = importdata([File.AbaqusPath, filesep, jobname, '_stress_tensor.dat']);
        stress_tensor_components = sortrows(stress_tensor_components);

        tmp = zeros(size(coords, 1), 8);
        for jj = 1:size(coords, 1)
            tmp(jj, :) = nanmean(stress_tensor_components(stress_tensor_components(:,1)==coords(jj,1), :), 1);
        end

        % Cauchy stress tensor, given in the deformed configuration, values at the nodes.
        Sxx_nodes = tmp(:, 2);
        Syy_nodes = tmp(:, 3);
        Szz_nodes = tmp(:, 4);
        Sxy_nodes = tmp(:, 5);
        Sxz_nodes = tmp(:, 6);
        Syz_nodes = tmp(:, 7);
        P_nodes   = tmp(:, 8);

    end

    X_def_nodes = coords(:, 11);
    Y_def_nodes = coords(:, 12);
    Z_def_nodes = coords(:, 13) + zoffset;

    % Similar to the stress tensor, the normals are defined at each face of each element, not at each node. Here, we compute 
    % the normals at each node.

    % Calculate the normals based on the deformed configuration. In the matrx "normals", each row represents a surface face of 
    % an element that is on the surface. Columns 1, 2 and 3 represent the center of the face (location of the normal). Columns 
    % 4, 5 and 6 represent the components of the normal vector at that face.

    % Typically, one would calculate the normal to the surface as the cross product of the two vectors defined by three points 
    % on the surface. But, we could have elements with faces defined by four points. In that case, if the element is distorted, 
    % the four points might not define a plane anymore. To take this case into account, we calculate the normal vector with a 
    % least square fit of the points of the face to a plane.

    normals = zeros(size(elements_surface, 1), 6);

    elements_surface_area = zeros(size(elements_surface, 1), 2);
    elements_surface_area(:, 1) = elements_surface(:, 1);

    nodes_surface_area = zeros(size(coords, 1), 2);
    nodes_surface_area(:, 1) = coords(:, 1);

    for jj = 1:size(elements_surface, 1)

        % Nodes of an element at the surface.
        tmp = elements_surface(jj, 2:end);
        tmp = tmp(tmp ~= 0);

        % Check if at least three of the nodes are on the surface and are on the subset of nodes where the displacements 
        % have been imposed.
        good_face = 1;
        for kk=1:numel(tmp)
            if isempty(coords(coords(:, 1) == tmp(kk), 11:13))
                good_face = 0;
            end
        end

        if good_face

            % Calculate the area of the face of the element on the surface.
            x1 = [coords(coords(:, 1)==tmp(1), 11), coords(coords(:, 1)==tmp(1), 12), coords(coords(:, 1)==tmp(1), 13)];
            x2 = [coords(coords(:, 1)==tmp(2), 11), coords(coords(:, 1)==tmp(2), 12), coords(coords(:, 1)==tmp(2), 13)];
            x3 = [coords(coords(:, 1)==tmp(3), 11), coords(coords(:, 1)==tmp(3), 12), coords(coords(:, 1)==tmp(3), 13)];

            elements_surface_area(jj, 2) = norm(cross((x3-x1), (x2-x1)))/2;

            % For all the nodes of the element, calculate the coordinate sums.
            X = 0;
            Y = 0;
            Z = 0;
            XY = 0;
            XZ = 0;
            YZ = 0;
            X2 = 0;
            Y2 = 0;
            Z2 = 0;

            for kk=1:numel(tmp)

                X = X + coords(coords(:, 1) == tmp(kk), 11);
                Y = Y + coords(coords(:, 1) == tmp(kk), 12);
                Z = Z + coords(coords(:, 1) == tmp(kk), 13);

                XY = XY + coords(coords(:, 1) == tmp(kk), 11)*coords(coords(:, 1) == tmp(kk), 12);
                XZ = XZ + coords(coords(:, 1) == tmp(kk), 11)*coords(coords(:, 1) == tmp(kk), 13);
                YZ = YZ + coords(coords(:, 1) == tmp(kk), 12)*coords(coords(:, 1) == tmp(kk), 13);

                X2 = X2 + coords(coords(:, 1) == tmp(kk), 11)^2;
                Y2 = Y2 + coords(coords(:, 1) == tmp(kk), 12)^2;
                Z2 = Z2 + coords(coords(:, 1) == tmp(kk), 13)^2;

                % Calculate the center of the face corresponding to the normal. 
                normals(jj, 1:3) = normals(jj, 1:3) + coords(coords(:, 1) == tmp(kk), 11:13);

                % Assign the correspoinding area to each node of the surface.
                nodes_surface_area(nodes_surface_area(:, 1)==tmp(kk), 2) = ...
                    nodes_surface_area(nodes_surface_area(:, 1)==tmp(kk), 2) + elements_surface_area(jj, 2)/3;

            end

            normals(jj, 1:3) = normals(jj, 1:3)/numel(tmp);

            % The coordinates of the normal vectors are given by the following least squares expressions.
            normals(jj, 4) = (-X*Y2*Z2 + X*YZ^2 + XY*Y*Z2 - XY*YZ*Z - XZ*Y*YZ + XZ*Y2*Z)/...
                             (X2*Y2*Z2 - X2*YZ^2 - XY^2*Z2 + 2*XY*XZ*YZ - XZ^2*Y2);
            normals(jj, 5) = ( X*XY*Z2 - X*XZ*YZ - X2*Y*Z2 + X2*YZ*Z - XY*XZ*Z + XZ^2*Y)/...
                             (X2*Y2*Z2 - X2*YZ^2 - XY^2*Z2 + 2*XY*XZ*YZ - XZ^2*Y2);
            normals(jj, 6) = (-X*XY*YZ + X*XZ*Y2 + X2*Y*YZ - X2*Y2*Z + XY^2*Z - XY*XZ*Y)/...
                             (X2*Y2*Z2 - X2*YZ^2 - XY^2*Z2 + 2*XY*XZ*YZ - XZ^2*Y2);

            normals(jj, 4:6) = normals(jj, 4:6)/sqrt(normals(jj, 4)^2 + normals(jj, 5)^2 + normals(jj, 6)^2);

            % Calculate the sign that makes the normal point outwards.
            % Lisf of all the nodes of the element:
            tmp = elements_surface_all(jj, 2:end);
            x_cm = [0, 0, 0];                   % Coordinates of the center of mass of the element.

            for kk=1:numel(tmp)

                % Calculate the center of the face corresponding to the normal. 
                x_cm = x_cm + deformed_coord_all_elements(deformed_coord_all_elements(:, 1) == tmp(kk), 2:4);

            end

            x_cm = x_cm/numel(tmp);

            v_cf_cm = x_cm - normals(jj, 1:3);
            signnormal = -sign(dot(v_cf_cm, normals(jj, 4:6)));
            normals(jj, 4:6) = normals(jj, 4:6)*signnormal;

        end

    end

    % The normals have been calculated at the center of the deformed faces. We interpolate the normals to the nodes.
%     ind = all(normals, 2);
    ind = any(normals(:, 4:6), 2);
    normals = normals(ind, :);

    if any(figs_to_plot == 13)

        % Plot a scatter3 of the face center locations and a quiver3 of the normals.
        fig13 = figure(13);

        set(fig13, 'Visible', 'on');
        ax13 = axes('Parent', fig13);

        scatter3(ax13, normals(:, 1), normals(:, 2), normals(:, 3), 10, 'cr', 'filled', 'DisplayName', "Element's face center");
        hold(ax13, 'on');
        quiver3(ax13, normals(:, 1), normals(:, 2), normals(:, 3), normals(:, 4), normals(:, 5), normals(:, 6), 3, ...
                'DisplayName', 'Local normal');
        axis(ax13, 'image');

        xlabel(ax13, 'X [\mum]');
        ylabel(ax13, 'Y [\mum]');
        zlabel(ax13, 'Z [\mum]');

        legend(ax13, 'location', 'northeast');

    end

    % Interpolate at the node deformed locations.
    FNx_nodes = scatteredInterpolant(normals(:,1), normals(:,2), normals(:,3), normals(:,4), 'natural', 'linear');
    FNy_nodes = scatteredInterpolant(normals(:,1), normals(:,2), normals(:,3), normals(:,5), 'natural', 'linear');
    FNz_nodes = scatteredInterpolant(normals(:,1), normals(:,2), normals(:,3), normals(:,6), 'natural', 'linear');

    Nx_nodes = FNx_nodes(X_def_nodes, Y_def_nodes, Z_def_nodes-zoffset);
    Ny_nodes = FNy_nodes(X_def_nodes, Y_def_nodes, Z_def_nodes-zoffset);
    Nz_nodes = FNz_nodes(X_def_nodes, Y_def_nodes, Z_def_nodes-zoffset);

    normN = sqrt(Nx_nodes.*Nx_nodes + Ny_nodes.*Ny_nodes + Nz_nodes.*Nz_nodes);
    Nx_nodes = Nx_nodes./normN;
    Ny_nodes = Ny_nodes./normN;
    Nz_nodes = Nz_nodes./normN;

    if isfield(Settings, "constrained") && strcmp(Settings.constrained, 'Constrained')

        % In the constrained method, the tractions are readitly calculated.
        tzx_nodes = tractions(:, 2)./nodes_surface_area(:, 2);
        tzy_nodes = tractions(:, 3)./nodes_surface_area(:, 2);
        tzz_nodes = tractions(:, 4)./nodes_surface_area(:, 2);

    else

        % In the unconstrained method, calculate the tractions normal to the external surface from the stress tensor.
        tzx_nodes = Sxx_nodes.*Nx_nodes + Sxy_nodes.*Ny_nodes + Sxz_nodes.*Nz_nodes;
        tzy_nodes = Sxy_nodes.*Nx_nodes + Syy_nodes.*Ny_nodes + Syz_nodes.*Nz_nodes;
        tzz_nodes = Sxz_nodes.*Nx_nodes + Syz_nodes.*Ny_nodes + Szz_nodes.*Nz_nodes;

    end

    % Equilibrate the tractions, if the field is in equilibrium.
    tzx_nodes = tzx_nodes - nanmean(tzx_nodes(:));
    tzy_nodes = tzy_nodes - nanmean(tzy_nodes(:));
    tzz_nodes = tzz_nodes - nanmean(tzz_nodes(:));

    % Calculate the normal tractions.
    tn_nodes = Nx_nodes.*tzx_nodes + Ny_nodes.*tzy_nodes + Nz_nodes.*tzz_nodes;

    % Calculate the components of the tangential tractions.
    tt_x_nodes = tzx_nodes - tn_nodes.*Nx_nodes;
    tt_y_nodes = tzy_nodes - tn_nodes.*Ny_nodes;
    tt_z_nodes = tzz_nodes - tn_nodes.*Nz_nodes;

    tt_nodes = sqrt(tt_x_nodes.^2 + tt_y_nodes.^2 + tt_z_nodes.^2);
        
    % Save the tracton field.
    f = fopen([File.TracPath, filesep, 'Traction_abaqus_abaqus_grid_', File.Base.Beads, ...
               sprintf(['%0', num2str(pad), 'd'], tt), '.dat'], 'w');

    for jj = 1:size(DX_nodes_q, 1)
        fprintf(f, '%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n', ...
                coords(jj, 1), X_def_nodes(jj), Y_def_nodes(jj), Z_def_nodes(jj), ...
                DX_nodes_q(jj), DY_nodes_q(jj), DZ_nodes_q(jj), ...
                tzx_nodes(jj), tzy_nodes(jj), tzz_nodes(jj), tn_nodes(jj), tt_nodes(jj), ...
                tt_x_nodes(jj), tt_y_nodes(jj), tt_z_nodes(jj));
    end

    fclose(f);

    if ~isfield(Settings, "constrained") || ~strcmp(Settings.constrained, 'Constrained')

        f = fopen([File.TracPath, filesep, 'Stress_tensor_abaqus_abaqus_grid_', File.Base.Beads, ...
                   sprintf(['%0', num2str(pad), 'd'], tt), '.dat'], 'w');

        for jj = 1:size(DX_nodes_q, 1)
            fprintf(f, '%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n', ...
                    X_def_nodes(jj), Y_def_nodes(jj), Z_def_nodes(jj), ...
                    Sxx_nodes(jj), Syy_nodes(jj), Szz_nodes(jj), Sxy_nodes(jj), Sxz_nodes(jj), Syz_nodes(jj), P_nodes(jj));
        end

        fclose(f);

    end

    % Plot the traction field.
    if any(figs_to_plot == 5)

        % Plot a 3D quiver of the tractions, in the deformed configuration.
        fig5 = figure(5);

        set(fig5, 'Visible', 'on');
        ax5 = axes('Parent', fig5);

        quiver3(ax5, X_def_nodes, Y_def_nodes, Z_def_nodes, tzx_nodes, tzy_nodes, tzz_nodes, 3, ...
                'DisplayName', 'Local traction');
        axis(ax5, 'equal');

        xlabel(ax5, 'X [\mum]');
        ylabel(ax5, 'Y [\mum]');
        zlabel(ax5, 'Z [\mum]');

        legend(ax5, 'location', 'northeast');

    end

    if any(figs_to_plot == 6)

        % Plot a scatter3 of the normal tractions, in the deformed configuration.
        fig6 = figure(6);

        set(fig6, 'Visible', 'on');
        ax6 = axes('Parent', fig6);

        scatter3(ax6, X_def_nodes, Y_def_nodes, Z_def_nodes, 170, tn_nodes, 's', 'filled');
        set(ax6, 'Clim', [-1000, 1000]);
        cb6 = colorbar(ax6);
        axis(ax6, 'equal');

        xlabel(ax6, 'X [\mum]');
        ylabel(ax6, 'Y [\mum]');
        zlabel(ax6, 'Z [\mum]');

        ylabel(cb6, 't_n [Pa]');

    end

    if any(figs_to_plot == 7)

        % Plot a scatter3 of the tangential tractions, in the deformed configuration.
        fig7 = figure(7);

        set(fig7, 'Visible', 'on');
        ax7 = axes('Parent', fig7);

        scatter3(ax7, X_def_nodes, Y_def_nodes, Z_def_nodes, 170, tt_nodes, 's', 'filled');
        set(ax7, 'Clim', [-1000, 1000]);
        cb7 = colorbar(ax7);
        axis(ax7, 'equal');

        xlabel(ax7, 'X [\mum]');
        ylabel(ax7, 'Y [\mum]');
        zlabel(ax7, 'Z [\mum]');

        ylabel(cb7, 't_t [Pa]');

    end

    % For each voxel of the original fluorescence image, find if any mesh node lies inside it. If so, calculate the mean of the 
    % displacement of all the nodes inside the voxel. Use the PIV displacement not the F.E.M. displacement. For each voxel with 
    % nodes inside, save its coordinates and displacement.
    DX_PIV_q = NaN(sy, sx, sz);
    DY_PIV_q = NaN(sy, sx, sz);
    DZ_PIV_q = NaN(sy, sx, sz);

    surf_mask = zeros(sy, sx, sz);

    parfor(ll=1:sx*sy*sz, n_th)
%     for ll=1:sx*sy*sz

        ind = (X_def_nodes>=(x_fluo(ll) - 1/2)*Settings.PixelSizeXY)& ...
              (X_def_nodes<=(x_fluo(ll) + 1/2)*Settings.PixelSizeXY)& ...
              (Y_def_nodes>=(y_fluo(ll) - 1/2)*Settings.PixelSizeXY)& ...
              (Y_def_nodes<=(y_fluo(ll) + 1/2)*Settings.PixelSizeXY)& ...
              (Z_def_nodes>=(z_fluo(ll) - 1/2)*Settings.PixelSizeZ )& ...
              (Z_def_nodes<=(z_fluo(ll) + 1/2)*Settings.PixelSizeZ );

        if any(ind)

            DX_PIV_q(ll) = nanmean(DX_nodes_q(ind));
            DY_PIV_q(ll) = nanmean(DY_nodes_q(ind));
            DZ_PIV_q(ll) = nanmean(DZ_nodes_q(ind));

            surf_mask(ll) = 1;

        end

    end

    % Save the displacements of the nodes in the image grid.
    f = fopen([File.TracPath, filesep, 'Displacements_PIV_Image_grid_', File.Base.Beads, ...
               sprintf(['%0', num2str(pad), 'd'], tt), '.dat'], 'w');

    ind = find(surf_mask);

    for kk = 1:numel(ind)

        jj = ind(kk) ;

        fprintf(f, '%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n', ...
                x_img(jj), y_img(jj), z_img(jj), DX_PIV_q(jj), DY_PIV_q(jj), DZ_PIV_q(jj));

    end

    fclose(f);

    clear("DX_nodes_q", "DY_nodes_q", "DZ_nodes_q", "DX_PIV_q", "DY_PIV_q", "DZ_PIV_q", "surf_mask", "ind");

    % Interpolate the tractions and normals in the PIV grid. For each box of the PIV, check if there are any mesh nodes lying 
    % inside. If so, calculate the average of the tractions and normal vector of all the nodes inside.
    dx = abs(X_PIV(1, 2, 1) - X_PIV(1, 1, 1));
    dy = abs(Y_PIV(2, 1, 1) - Y_PIV(1, 1, 1));
    dz = abs(Z_PIV(1, 1, 2) - Z_PIV(1, 1, 1));
            
    tzx_PIV_q = NaN(ny, nx, nz);
    tzy_PIV_q = NaN(ny, nx, nz);
    tzz_PIV_q = NaN(ny, nx, nz);
    Nx_PIV_q  = NaN(ny, nx, nz);
    Ny_PIV_q  = NaN(ny, nx, nz);
    Nz_PIV_q  = NaN(ny, nx, nz);

    surf_mask = zeros(ny, nx, nz);

    parfor(ll=1:nx*ny*nz, n_th)
%     for ll=1:nx*ny*nz

        ind = (X_def_nodes>=X_PIV(ll)-dx/2)&(X_def_nodes<=X_PIV(ll)+dx/2)&...
              (Y_def_nodes>=Y_PIV(ll)-dy/2)&(Y_def_nodes<=Y_PIV(ll)+dy/2)&...
              (Z_def_nodes>=Z_PIV(ll)-dz/2)&(Z_def_nodes<=Z_PIV(ll)+dz/2);

        if any(ind)

            tzx_PIV_q(ll) = nanmean(tzx_nodes(ind));
            tzy_PIV_q(ll) = nanmean(tzy_nodes(ind));
            tzz_PIV_q(ll) = nanmean(tzz_nodes(ind));
            Nx_PIV_q(ll)  = nanmean(Nx_nodes(ind));
            Ny_PIV_q(ll)  = nanmean(Ny_nodes(ind));
            Nz_PIV_q(ll)  = nanmean(Nz_nodes(ind));

            surf_mask(ll) = 1;

        end

    end

    normN = sqrt(Nx_PIV_q.*Nx_PIV_q + Ny_PIV_q.*Ny_PIV_q + Nz_PIV_q.*Nz_PIV_q);
    Nx_PIV_q = Nx_PIV_q./normN;
    Ny_PIV_q = Ny_PIV_q./normN;
    Nz_PIV_q = Nz_PIV_q./normN;

    tzx_PIV_q = tzx_PIV_q - nanmean(tzx_PIV_q(:));
    tzy_PIV_q = tzy_PIV_q - nanmean(tzy_PIV_q(:));
    tzz_PIV_q = tzz_PIV_q - nanmean(tzz_PIV_q(:));

    tn_PIV_q = Nx_PIV_q.*tzx_PIV_q + Ny_PIV_q.*tzy_PIV_q + Nz_PIV_q.*tzz_PIV_q;

    tt_x_PIV_q = tzx_PIV_q - tn_PIV_q.*Nx_PIV_q;
    tt_y_PIV_q = tzy_PIV_q - tn_PIV_q.*Ny_PIV_q;
    tt_z_PIV_q = tzz_PIV_q - tn_PIV_q.*Nz_PIV_q;

    tt_PIV_q = sqrt(tt_x_PIV_q.^2 + tt_y_PIV_q.^2 + tt_z_PIV_q.^2);

    % Save the tracton in the PIV grid.
    f = fopen([File.TracPath, filesep, 'Traction_abaqus_PIV_grid_', File.Base.Beads, ...
               sprintf(['%0', num2str(pad), 'd'], tt), '.dat'], 'w');

    ind = find(surf_mask);

    for kk = 1:numel(ind)

        jj = ind(kk);

        fprintf(f, '%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n', ...
                X_PIV(jj), Y_PIV(jj), Z_PIV(jj), tzx_PIV_q(jj), tzy_PIV_q(jj), tzz_PIV_q(jj), ...
                tt_x_PIV_q(jj), tt_y_PIV_q(jj), tt_z_PIV_q(jj), tt_PIV_q(jj), tn_PIV_q(jj), ...
                Nx_PIV_q(jj), Ny_PIV_q(jj), Nz_PIV_q(jj));

    end

    fclose(f);

    if any(figs_to_plot == 8)

        % Plot a scatter3 of the node locations, a quiver3 of the normals at the node locations, and a quiver3 of the 
        % normals at the PIV grid, in the deformed configuration.
        fig8 = figure(8);

        set(fig8, 'Visible', 'on');
        ax8 = axes('Parent', fig8);

        scatter3(ax8, X_def_nodes, Y_def_nodes, Z_def_nodes, 10, 'cr', 'filled', 'DisplayName', 'Mesh nodes');
        hold(ax8, 'on');
        quiver3(ax8, X_def_nodes, Y_def_nodes, Z_def_nodes, Nx_nodes(:), Ny_nodes(:), Nz_nodes(:), 3, ...
                'DisplayName', 'Local normals at the nodes');
        quiver3(ax8, X_PIV(:), Y_PIV(:), Z_PIV(:), Nx_PIV_q(:), Ny_PIV_q(:), Nz_PIV_q(:), 3, ...
                'DisplayName', 'Local normals at the PIV grid');
        axis(ax8, 'image');

        xlabel(ax8, 'X [\mum]');
        ylabel(ax8, 'Y [\mum]');
        zlabel(ax8, 'Z [\mum]');

        legend(ax8, 'location', 'northeast');

    end

    clear( "X_PIV", "Y_PIV", "Z_PIV", "dx", "dy", "dz", "tzx_PIV_q", "tzy_PIV_q", "tzz_PIV_q", ...
        "Nx_PIV_q", "Ny_PIV_q", "Nz_PIV_q", "surf_mask", "normN", ...
        "tn_PIV_q", "tt_x_PIV_q", "tt_y_PIV_q", "tt_z_PIV_q", "tt_PIV_q", "ind" ) ;

    % For each voxel of the original fluorescence image, find if any mesh node lies inside it. If so, calculate the mean of the 
    % tractions and normals of all the nodes inside the voxel. For each voxel with nodes inside, save its coordinates, traction 
    % components and normal.
    tzx_fluo_q = NaN(sy, sx, sz);
    tzy_fluo_q = NaN(sy, sx, sz);
    tzz_fluo_q = NaN(sy, sx, sz);
    Nx_fluo_q = NaN(sy, sx, sz);
    Ny_fluo_q = NaN(sy, sx, sz);
    Nz_fluo_q = NaN(sy, sx, sz);

    surf_mask = zeros(sy, sx, sz);

%     parfor(ll=1:sx*sy*sz, n_th)
    for ll=1:sx*sy*sz

        ind = (X_def_nodes>=(x_fluo(ll) - 1/2)*Settings.PixelSizeXY)& ...
              (X_def_nodes<=(x_fluo(ll) + 1/2)*Settings.PixelSizeXY)& ...
              (Y_def_nodes>=(y_fluo(ll) - 1/2)*Settings.PixelSizeXY)& ...
              (Y_def_nodes<=(y_fluo(ll) + 1/2)*Settings.PixelSizeXY)& ...
              (Z_def_nodes>=(z_fluo(ll) - 1/2)*Settings.PixelSizeZ )& ...
              (Z_def_nodes<=(z_fluo(ll) + 1/2)*Settings.PixelSizeZ );

        if any( ind )

            tzx_fluo_q(ll) = nanmean(tzx_nodes(ind));
            tzy_fluo_q(ll) = nanmean(tzy_nodes(ind));
            tzz_fluo_q(ll) = nanmean(tzz_nodes(ind));
            Nx_fluo_q(ll) = nanmean(Nx_nodes(ind));
            Ny_fluo_q(ll) = nanmean(Ny_nodes(ind));
            Nz_fluo_q(ll) = nanmean(Nz_nodes(ind));

            surf_mask(ll) = 1;

        end

    end

    normN = sqrt(Nx_fluo_q.*Nx_fluo_q + Ny_fluo_q.*Ny_fluo_q + Nz_fluo_q.*Nz_fluo_q);
    Nx_fluo_q = Nx_fluo_q./normN;
    Ny_fluo_q = Ny_fluo_q./normN;
    Nz_fluo_q = Nz_fluo_q./normN;

    tzx_fluo_q = tzx_fluo_q - nanmean(tzx_fluo_q(:));
    tzy_fluo_q = tzy_fluo_q - nanmean(tzy_fluo_q(:));
    tzz_fluo_q = tzz_fluo_q - nanmean(tzz_fluo_q(:));

    tn_fluo_q = Nx_fluo_q.*tzx_fluo_q + Ny_fluo_q.*tzy_fluo_q + Nz_fluo_q.*tzz_fluo_q;

    tt_x_fluo_q = tzx_fluo_q - tn_fluo_q.*Nx_fluo_q;
    tt_y_fluo_q = tzy_fluo_q - tn_fluo_q.*Ny_fluo_q;
    tt_z_fluo_q = tzz_fluo_q - tn_fluo_q.*Nz_fluo_q;

    tt_fluo_q = sqrt(tt_x_fluo_q.^2 + tt_y_fluo_q.^2 + tt_z_fluo_q.^2);

    % Save the tracton and normals in the image grid.
    f = fopen([File.TracPath, filesep, 'Traction_abaqus_Image_grid_', File.Base.Beads, ...
               sprintf(['%0', num2str(pad), 'd'], tt), '.dat'], 'w');

    ind = find( surf_mask ) ;

    for kk = 1:numel(ind)

        jj = ind(kk) ;

        fprintf(f, '%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n', ...
                x_fluo(jj)*Settings.PixelSizeXY, y_fluo(jj)*Settings.PixelSizeXY, z_fluo(jj)*Settings.PixelSizeZ, ...
                tzx_fluo_q(jj), tzy_fluo_q(jj), tzz_fluo_q(jj), tt_x_fluo_q(jj), tt_y_fluo_q(jj), tt_z_fluo_q(jj), ...
                tt_fluo_q(jj), tn_fluo_q(jj), Nx_fluo_q(jj), Ny_fluo_q(jj), Nz_fluo_q(jj));

    end

    fclose(f);

    if any(figs_to_plot == 9)

        % Plot a scatter3 of the node locations, and a quiver3 of the normals at the node locations, in the deformed 
        % configuration.
        fig9 = figure(9);

        set(fig9, 'Visible', 'on');
        ax9 = axes('Parent', fig9);

        scatter3(ax9, X_def_nodes, Y_def_nodes, Z_def_nodes, 10, 'cr', 'filled', 'DisplayName', 'Mesh nodes');
        hold(ax9, 'on');
        quiver3(ax9, X_def_nodes, Y_def_nodes, Z_def_nodes, Nx_nodes(:), Ny_nodes(:), Nz_nodes(:), 3, 'DisplayName', 'Normals');
        axis(ax9, 'image') ;

        xlabel(ax9, 'X [\mum]');
        ylabel(ax9, 'Y [\mum]');
        zlabel(ax9, 'Z [\mum]');

        legend(ax9, 'location', 'northeast');

    end

    clear("x_fluo", "y_fluo", "z_fluo", "x_fluo", "y_fluo", "z_fluo", "tzx_fluo_q", "tzy_fluo_q", "tzz_fluo_q", ...
          "Nx_fluo_q", "Ny_fluo_q", "Nz_fluo_q", "surf_mask", "normN", "tn_fluo_q", ...
          "tt_x_fluo_q", "tt_y_fluo_q", "tt_z_fluo_q", "tt_fluo_q", "ind");

    disp(' ');
    disp(['          *** End timepoint ,', num2str(tt), ' ***        ']);
    disp(' ');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
