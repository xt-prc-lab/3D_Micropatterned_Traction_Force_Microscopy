%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                   Solve the F.E.M. problem where the displacements are prescribed at a subset of the nodes                   %
%                  at the surface of the gel (3D PIV displacements of the gel surface as Boundary Conditions),                 %
%    and the tractions are constrained to the walls of the well (prescribe traction-free at the upper surface of the gel).     %
%                    The stiffness matrix of the unconstrained problem was previously calculated in Abaqus.                    %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       folder [str]: location of the unconstrained problem's stiffness matrix file.                                           %
%       jobname [str]: name of the Abaqus job.                                                                                 %
%       nodes_surf_bc [matrix]: numbers of the nodes of the surface where we apply the boundary conditions (might not be the   %
%                               same as the whole upper surface).                                                              %
%       def_surf_bc [matrix]: boundary conditions, i.e. measured deformations of the nodes in nodes_surf_bc, in microns.       %
%                             The columns are the x, y and z displacements, respectively.                                      %
%       nodes_surf [matrix]: nodes at the upper surface of the gel.                                                            %
%       nodes_bottom [matrix]: nodes at the bottom surface of the model. They are fixed (encastre B.C.).                       %
%       nodes_coordinates [matrix]: [node_number, x_0, y_0, z_0], for all the nodes.                                           %
%       elements_surface [matrix]: for all the elements with a face on the surface, the first column is the element number,    %
%                                  the second, third and fourth columns are the node numbers. The nodes that are not on the    %
%                                  surface are substituted by zero.                                                            %
%       elements_surface_all [matrix]: for all the elements with a face on the surface, the first column is the element        %
%                                      number, the second, third and fourth columns are the node numbers.                      %
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

function calculate_constrained_tractions(folder, jobname, nodes_surf_bc, def_surf_bc, nodes_surf, nodes_bottom, ...
                                         nodes_coordinates, elements_surface, elements_surface_all)

    %% Import the stiffness matrix of the un-constrained problem.
    % Output from this block:
    %   * k_rigidity: rigidity matrix of the unconstrained problem.

    rigidity_matrix = load([folder, filesep, jobname, '_STIF2_mtx_matrix_input.mtx']);

    % Reference on how the data is stored in the stiffness matrix file: 
    % https://xizou.github.io/2016/04/07/how-to-export-stiffness-matrix-from-abaqus/
    ind_pos = (rigidity_matrix(:, 1)>=0)&(rigidity_matrix(:, 3)>=0);
    ind_row = (rigidity_matrix(:,1)-1)*3 + rigidity_matrix(:, 2);
    ind_col = (rigidity_matrix(:,3)-1)*3 + rigidity_matrix(:, 4);

    k_rigidity_table = [[ind_row(ind_pos); ind_col(ind_pos)], ...
                        [ind_col(ind_pos); ind_row(ind_pos)], ...
                        [rigidity_matrix(ind_pos, 5); rigidity_matrix(ind_pos, 5)]];
    k_rigidity_table = unique(k_rigidity_table, 'rows');
    k_rigidity = sparse(k_rigidity_table(:, 1), k_rigidity_table(:, 2), k_rigidity_table(:, 3));

    clear rigidity_matrix ind_pos ind_row ind_col k_rigidity_table;

    %% Create vectors with the node+degree_of_freedom locations:
    % Outputs from this block:
    %   * nodes_free_surf = [node_number] part of the upper surfaces that is force-free.
    %   * nodes_well_surf = [node_number] part of the upper surface where forces are constrained.

    n_nodes = size(nodes_coordinates(:, 1), 1);       % Number of nodes.

    % Separate the nodes of the upper surface and the well.
    surf_z_pos = max(nodes_coordinates(nodes_surf, 4));

    nodes_free_surf = nodes_coordinates(nodes_coordinates(:, 4)==surf_z_pos, 1);
    nodes_well_surf = setdiff(nodes_surf, nodes_free_surf);

    nodes_dof_well_surf = [3*(nodes_well_surf-1)+1; 3*(nodes_well_surf-1)+2; 3*(nodes_well_surf-1)+3];

    %% Remove the degrees of freedom of the bottom nodes.
    nodes_dof_bottom = [3*(nodes_bottom-1)+1; 3*(nodes_bottom-1)+2; 3*(nodes_bottom-1)+3];

    k_rigidity_nobottom = k_rigidity;
    k_rigidity_nobottom(nodes_dof_bottom, :) = [];
    k_rigidity_nobottom(:, nodes_dof_bottom) = [];

    % Matrices S and L.
    S_nobottom = sparse([3*(nodes_surf_bc-1)+1; 3*(nodes_surf_bc-1)+2; 3*(nodes_surf_bc-1)+3], ...
                        [3*(nodes_surf_bc-1)+1; 3*(nodes_surf_bc-1)+2; 3*(nodes_surf_bc-1)+3], ...
                        ones(3*size(nodes_surf_bc, 1), 1), 3*n_nodes, 3*n_nodes);
    S_nobottom(nodes_dof_bottom, :) = [];
    S_nobottom(:, nodes_dof_bottom) = [];

    L_nobottom = speye(3*n_nodes, 3*n_nodes);
    L_nobottom(3*(nodes_well_surf-1)+1, 3*(nodes_well_surf-1)+1) = 0;
    L_nobottom(3*(nodes_well_surf-1)+2, 3*(nodes_well_surf-1)+2) = 0;
    L_nobottom(3*(nodes_well_surf-1)+3, 3*(nodes_well_surf-1)+3) = 0;

    L_nobottom([nodes_dof_bottom; nodes_dof_well_surf], :) = [];
    L_nobottom(:, nodes_dof_bottom) = [];

    %% Ensemble the modified stiffness matrix.

    k = 1;  % Energy constant.

    % Number of nodes not belonging to the bottom.
    n_nobottom = n_nodes-size(nodes_bottom, 1);

    % Modified rigidity matrix excluding the nodes belonging to the bottom.
    k_mod_nobottom = sparse(4*3*n_nobottom-3*size(nodes_well_surf, 1), 4*3*n_nobottom-3*size(nodes_well_surf, 1));

    k_mod_nobottom(1:3*n_nobottom, 1:3*n_nobottom) = k*S_nobottom;
    k_mod_nobottom(1:3*n_nobottom, (2*3*n_nobottom+1):3*3*n_nobottom) = k_rigidity_nobottom';
    k_mod_nobottom((3*n_nobottom+1):2*3*n_nobottom, (2*3*n_nobottom+1):3*3*n_nobottom) = -speye(3*n_nobottom, 3*n_nobottom);
    k_mod_nobottom((3*n_nobottom+1):2*3*n_nobottom, (3*3*n_nobottom+1):end) = L_nobottom';
    k_mod_nobottom((2*3*n_nobottom+1):3*3*n_nobottom, 1:3*n_nobottom) = k_rigidity_nobottom;
    k_mod_nobottom((2*3*n_nobottom+1):3*3*n_nobottom, (3*n_nobottom+1):2*3*n_nobottom) = -speye(3*n_nobottom, 3*n_nobottom);
    k_mod_nobottom((3*3*n_nobottom+1):end, (3*n_nobottom+1):2*3*n_nobottom) = L_nobottom;

    % Forcing vector, b, excluding the nodes belonging to the bottom.
    b_nobottom = sparse(4*3*n_nodes, 1);
    b_nobottom(3*(nodes_surf_bc-1)+1) = k*def_surf_bc(:, 1);
    b_nobottom(3*(nodes_surf_bc-1)+2) = k*def_surf_bc(:, 2);
    b_nobottom(3*(nodes_surf_bc-1)+3) = k*def_surf_bc(:, 3);

    b_nobottom([nodes_dof_bottom; nodes_dof_bottom+3*n_nodes; nodes_dof_bottom+2*3*n_nodes; [nodes_dof_bottom; nodes_dof_well_surf]+3*3*n_nodes]) = [];

    %% Ensemble the force-gradient dumping matrix.
    alpha = 5e-10;      % Weight penality constant for the force gradients.

    k_alpha = sparse(3*n_nodes, 3*n_nodes);

    for jj=1:size(elements_surface, 1)

        % Nodes of an element at the surface.
        tmp = elements_surface( jj, 2:end ) ;
        tmp = tmp( tmp ~= 0 ) ;

        if numel(tmp)==4

            % Node of the element not on the surface.
            tmp4 = tmp(4);
            tmp = tmp(1:3);

        else

            % Node of the element not on the surface.
            tmp4 = setdiff(elements_surface_all( jj, 2:end ), tmp);

        end

        % Check if at least three of the nodes are on the surface and are on nodes_well_surf (surface of the well).
        good_face = 1 ;
        for kk=1:numel(tmp)
            if ~any(nodes_well_surf == tmp(kk))
                good_face = 0 ;
            end
        end

        if good_face

            % Calculate the cooirdinates of the three nodes of the face of the element on the surface.
            x1 = nodes_coordinates(nodes_coordinates(:, 1)==tmp(1), 2:4);
            x2 = nodes_coordinates(nodes_coordinates(:, 1)==tmp(2), 2:4);
            x3 = nodes_coordinates(nodes_coordinates(:, 1)==tmp(3), 2:4);

            % Calculate the coordinates of the other node.
            x4 = nodes_coordinates(nodes_coordinates(:, 1)==tmp4, 2:4);

            % Calculate the area of the face of the element on the surface.
            elements_surface_area = norm(cross((x3-x1), (x2-x1)))/2;

            % Calculate the derivatives of the shape functions for the element.
            [C11, C12, C13, C14, C21, C22, C23, C24, C31, C32, C33, C34, C41, C42, C43, C44] = ...
            	coeffs_penalization_function(x1(1), x2(1), x3(1), x4(1), x1(2), x2(2), x3(2), x4(2), x1(3), x2(3), x3(3), x4(3));

            % Fill in the terms of k_alpha:
            for ll=1:3

                k_alpha(3*(tmp(1)-1)+ll, 3*(tmp(1)-1)+ll) = k_alpha(3*(tmp(1)-1)+ll, 3*(tmp(1)-1)+ll) + ...
                                                            C11*2*elements_surface_area/3;
                k_alpha(3*(tmp(1)-1)+ll, 3*(tmp(2)-1)+ll) = k_alpha(3*(tmp(1)-1)+ll, 3*(tmp(2)-1)+ll) + ...
                                                            C12*2*elements_surface_area/3;
                k_alpha(3*(tmp(1)-1)+ll, 3*(tmp(3)-1)+ll) = k_alpha(3*(tmp(1)-1)+ll, 3*(tmp(3)-1)+ll) + ...
                                                            C13*2*elements_surface_area/3;
                k_alpha(3*(tmp(1)-1)+ll, 3*(tmp4  -1)+ll) = k_alpha(3*(tmp(1)-1)+ll, 3*(tmp4-1)+ll) + ...
                                                            C14*2*elements_surface_area/3;

                k_alpha(3*(tmp(2)-1)+ll, 3*(tmp(1)-1)+ll) = k_alpha(3*(tmp(2)-1)+ll, 3*(tmp(1)-1)+ll) + ...
                                                            C21*2*elements_surface_area/3;
                k_alpha(3*(tmp(2)-1)+ll, 3*(tmp(2)-1)+ll) = k_alpha(3*(tmp(2)-1)+ll, 3*(tmp(2)-1)+ll) + ...
                                                            C22*2*elements_surface_area/3;
                k_alpha(3*(tmp(2)-1)+ll, 3*(tmp(3)-1)+ll) = k_alpha(3*(tmp(2)-1)+ll, 3*(tmp(3)-1)+ll) + ...
                                                            C23*2*elements_surface_area/3;
                k_alpha(3*(tmp(2)-1)+ll, 3*(tmp4  -1)+ll) = k_alpha(3*(tmp(2)-1)+ll, 3*(tmp4-1)+ll) + ...
                                                            C24*2*elements_surface_area/3;

                k_alpha(3*(tmp(3)-1)+ll, 3*(tmp(1)-1)+ll) = k_alpha(3*(tmp(3)-1)+ll, 3*(tmp(1)-1)+ll) + ...
                                                            C31*2*elements_surface_area/3;
                k_alpha(3*(tmp(3)-1)+ll, 3*(tmp(2)-1)+ll) = k_alpha(3*(tmp(3)-1)+ll, 3*(tmp(2)-1)+ll) + ...
                                                            C32*2*elements_surface_area/3;
                k_alpha(3*(tmp(3)-1)+ll, 3*(tmp(3)-1)+ll) = k_alpha(3*(tmp(3)-1)+ll, 3*(tmp(3)-1)+ll) + ...
                                                            C33*2*elements_surface_area/3;
                k_alpha(3*(tmp(3)-1)+ll, 3*(tmp4  -1)+ll) = k_alpha(3*(tmp(3)-1)+ll, 3*(tmp4-1)+ll) + ...
                                                            C34*2*elements_surface_area/3;

                k_alpha(3*(tmp4  -1)+ll, 3*(tmp(1)-1)+ll) = k_alpha(3*(tmp4-1)+ll, 3*(tmp(1)-1)+ll) + ...
                                                            C41*2*elements_surface_area/3;
                k_alpha(3*(tmp4  -1)+ll, 3*(tmp(2)-1)+ll) = k_alpha(3*(tmp4-1)+ll, 3*(tmp(2)-1)+ll) + ...
                                                            C42*2*elements_surface_area/3;
                k_alpha(3*(tmp4  -1)+ll, 3*(tmp(3)-1)+ll) = k_alpha(3*(tmp4-1)+ll, 3*(tmp(3)-1)+ll) + ...
                                                            C43*2*elements_surface_area/3;
                k_alpha(3*(tmp4  -1)+ll, 3*(tmp4  -1)+ll) = k_alpha(3*(tmp4-1)+ll, 3*(tmp4-1)+ll) + ...
                                                            C44*2*elements_surface_area/3;

            end

        end

    end

    k_alpha = alpha*k_alpha;

    k_alpha_nobottom = k_alpha;
    k_alpha_nobottom(nodes_dof_bottom, :) = [];
    k_alpha_nobottom(:, nodes_dof_bottom) = [];

    % Modified rigidity matrix excluding the nodes belonging to the bottom.
    k_mod_nobottom((3*n_nobottom+1):2*3*n_nobottom, (3*n_nobottom+1):2*3*n_nobottom) = k_alpha_nobottom;

    %% Solve the problem.
    x_nobottom = k_mod_nobottom\b_nobottom;

    %% Save the data.
    nodes_numbers_nobottom = nodes_coordinates(:, 1);
    nodes_numbers_nobottom(nodes_bottom) = [];

    nodes_def = zeros(size(nodes_coordinates));
    nodes_def(:, 1) = nodes_coordinates(:, 1);
    nodes_def(nodes_numbers_nobottom, 2:4) = ...
                [x_nobottom(1:3:3*n_nobottom), x_nobottom(2:3:3*n_nobottom), x_nobottom(3:3:3*n_nobottom)];

    nodes_trac = zeros(size(nodes_coordinates));
    nodes_trac(:, 1) = nodes_coordinates(:, 1);
    nodes_trac(nodes_numbers_nobottom, 2:4) = ...
                [x_nobottom((3*n_nobottom+1):3:6*n_nobottom), ...
                 x_nobottom((3*n_nobottom+2):3:6*n_nobottom), ...
                 x_nobottom((3*n_nobottom+3):3:6*n_nobottom)];

    % Save (...)_deformed_coord.dat
    nodes_def_coordinates = nodes_coordinates;
    nodes_def_coordinates(:, 2:end) = nodes_def_coordinates(:, 2:end) + nodes_def(:, 2:end);

    f = fopen([folder, filesep, jobname, '_deformed_coord.dat'], 'w');
    for jj = 1:size(nodes_surf_bc, 1)
        fprintf(f, '%d\t %10.8E\t %10.8E\t %10.8E\t \n', ...
            nodes_def_coordinates(nodes_surf_bc(jj), 1), nodes_def_coordinates(nodes_surf_bc(jj), 2), ...
            nodes_def_coordinates(nodes_surf_bc(jj), 3), nodes_def_coordinates(nodes_surf_bc(jj), 4));
    end
    fclose(f);

    % Save (...)_reactions.dat
    f = fopen([folder, filesep, jobname, '_reactions.dat'], 'w');
    for jj = 1:size(nodes_surf_bc, 1)
        fprintf(f, '%d\t %10.8E\t %10.8E\t %10.8E\t \n', ...
            nodes_trac(nodes_surf_bc(jj), 1), nodes_trac(nodes_surf_bc(jj), 2), ...
            nodes_trac(nodes_surf_bc(jj), 3), nodes_trac(nodes_surf_bc(jj), 4));
    end
    fclose(f);

    % Save (...)_deformed_coord_all_elements.dat
    f = fopen([folder, filesep, jobname, '_deformed_coord_all_elements.dat'], 'w');
    for jj = 1:n_nodes
        fprintf(f, '%d\t %10.8E\t %10.8E\t %10.8E\t \n', ...
            nodes_def_coordinates(jj, 1), nodes_def_coordinates(jj, 2), ...
            nodes_def_coordinates(jj, 3), nodes_def_coordinates(jj, 4));
    end
    fclose(f);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%