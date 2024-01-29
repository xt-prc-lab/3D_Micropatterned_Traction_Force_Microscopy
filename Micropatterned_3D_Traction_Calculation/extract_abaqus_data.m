%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                                   Extract the node and element data of an Abaqus job file.                                   %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       file [str]: Abaqus job file.                                                                                           %
%       set_name [str]: name of the set of nodes where the Boundary Conditions (PIV displacements) are applied.                %
%       set_fixed_name [str]: name of the set of fixed nodes (encastre).                                                       %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       nodes_coordinates [matrix]: nodes belonging to the whole model. The first column is the node number, the second, third %
%                                   and fourth columns are their x, y and z coordinates, respectivelly.                        %
%       elements_nodes [matrix]: elements belonging to the whole model. The first column is the element number, the following  %
%                                columns are the numbers of the nodes that form the element.                                   %
%       nodes_to_impose [matrix]: number of the nodes where the Boundary Conditions (PIV displacements) are applied.           %
%       elements_surface_nodes [matrix]: number of the elements belonging to the B.C. set.                                     %
%       nodes_fixed [matrix]: number of the fixed nodes (encastre).                                                            %
%                                                                                                                              %
%   Last Revison Date: 22/01/2024                                                                                              %
%   Author: Manuel Gomez Gonzalez                                                                                              %
%                                                                                                                              %
%   References:                                                                                                                %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nodes_coordinates, elements_nodes, nodes_to_impose, elements_surface_nodes, nodes_fixed] = ...
            extract_abaqus_data(file, set_name, set_fixed_name)

    % Open the Abaqus model and import all the text lines.
    fid = fopen(file,'r');
    MyTextFile = textscan(fid,'%s','delimiter','\n','whitespace', '');
    fclose(fid);
    MyTextFile = [MyTextFile{:}];

    % Find the first and last line of the definition of the nodes of the part model.
    str = '*Part, name=';
    ind = find(strncmpi(MyTextFile, str, numel(str)));
    ind0_nodes_part = ind(1) + 2;

    str = '*Element, type=';
    ind = find(strncmpi(MyTextFile, str, numel(str)));
    ind = ind(ind>ind0_nodes_part);
    indf_nodes_part = ind(1) - 1;

    % Find the first and last line of the definition of the elements of the part model.
    ind0_elements_part = indf_nodes_part + 2;

    str = '*Nset, nset=';
    ind = find(strncmpi(MyTextFile, str, numel(str)));
    ind = ind(ind>ind0_elements_part);
    indf_elements_part = ind(1) - 1;

    % Find the first and last line of the definition of the nodes of the fixed B.C. set.
    str = ['*Nset, nset=', set_fixed_name, ', instance='];
    ind = find(strncmpi(MyTextFile, str, numel(str)));
    ind0_fixed_nodes_set = ind(1) + 1;

    str = ['*Elset, elset=', set_fixed_name, ', instance='];
    ind = find(strncmpi(MyTextFile, str, numel(str)));
    ind = ind(ind>ind0_fixed_nodes_set);
    indf_fixed_nodes_set = ind(1) - 1;

    % Find the first and last line of the definition of the nodes of the B.C. set.
    str = ['*Nset, nset=', set_name];
    ind = find(strncmpi(MyTextFile, str, numel(str)));
    ind0_nodes_set = ind(1) + 1;

    str = ['*Elset, elset=', set_name];
    ind = find(strncmpi(MyTextFile, str, numel(str)));
    ind = ind(ind>ind0_nodes_set);
    indf_nodes_set = ind(1) - 1;

    % Find the first and last line of the definition of the elements of the B.C. set.
    ind0_elements_set = indf_nodes_set + 2;

    str = '*';
    ind = find(strncmpi(MyTextFile, str, numel(str)));
    ind = ind(ind>ind0_elements_set);
    indf_elements_set = ind(1) - 1;

    % Create a matrix containing the nodes belonging to the whole model. 
    % The first column is the node number, the second, third and fourth columns
    % are the x, y and z coordinates, respectivelly.
    temp = textscan(MyTextFile{ind0_nodes_part}, '%s', 'Delimiter', ',');
    temp = str2double(temp{1}');
    nodes_coordinates = zeros(indf_nodes_part-ind0_nodes_part+1, numel(temp));

    for ii=1:indf_nodes_part-ind0_nodes_part+1

        temp = textscan(MyTextFile{ii+ind0_nodes_part-1}, '%s', 'Delimiter', ',');
        nodes_coordinates(ii, :) = str2double(temp{1}');

    end

    % Create a matrix containing the elements belonging to the whole model. 
    % The first column is the element number, 
    % the following columns are the numbers of the nodes that form the element.
    temp = textscan(MyTextFile{ind0_elements_part}, '%s');
    temp = str2double(temp{1}');
    elements_nodes = zeros(indf_elements_part-ind0_elements_part+1, numel(temp));

    for ii=1:indf_elements_part-ind0_elements_part+1

        temp = textscan(MyTextFile{ii+ind0_elements_part-1}, '%s', 'Delimiter', ',');
        elements_nodes(ii, :) = str2double(temp{1}');

    end

    % Create a matrix containing the number of the nodes belonging to the fixed B.C. set.
    temp = textscan(MyTextFile{ind0_fixed_nodes_set}, '%s', 'Delimiter', ',');
    temp = str2double(temp{1}');
    temp2 = numel(temp);
    nodes_fixed = zeros(indf_fixed_nodes_set-ind0_fixed_nodes_set+1, temp2);

    for ii=1:indf_fixed_nodes_set-ind0_fixed_nodes_set+1

        temp = textscan(MyTextFile{ii+ind0_fixed_nodes_set-1}, '%s', 'Delimiter', ',');
        temp = str2double(temp{1}');
        nodes_fixed(ii, 1:numel(temp)) = temp;

    end

    nodes_fixed = nodes_fixed';
    nodes_fixed = nodes_fixed(:);
    nodes_fixed = nodes_fixed(1:end-(temp2-numel(temp)));

    % Create a matrix containing the number of the nodes belonging to the B.C. set.
    temp = textscan(MyTextFile{ind0_nodes_set}, '%s', 'Delimiter', ',');
    temp = str2double(temp{1}');
    temp2 = numel(temp);
    nodes_to_impose = zeros(indf_nodes_set-ind0_nodes_set+1, temp2);

    for ii=1:indf_nodes_set-ind0_nodes_set+1

        temp = textscan(MyTextFile{ii+ind0_nodes_set-1}, '%s', 'Delimiter', ',');
        temp = str2double(temp{1}');
        nodes_to_impose(ii, 1:numel(temp)) = temp;

    end

    nodes_to_impose = nodes_to_impose';
    nodes_to_impose = nodes_to_impose(:);
    nodes_to_impose = nodes_to_impose(1:end-(temp2-numel(temp)));

    % Create a matrix containing the number of the elements belonging to the B.C. set.
    temp = textscan(MyTextFile{ind0_elements_set}, '%s', 'Delimiter', ',');
    temp = str2double(temp{1}');
    temp2 = numel(temp);
    elements_surface_nodes = zeros(indf_elements_set-ind0_elements_set+1, temp2);

    for ii=1:indf_elements_set-ind0_elements_set+1

        temp = textscan(MyTextFile{ii+ind0_elements_set-1}, '%s', 'Delimiter', ',');
        temp = str2double(temp{1}');
        elements_surface_nodes(ii, 1:numel(temp)) = temp;

    end

    elements_surface_nodes = elements_surface_nodes';
    elements_surface_nodes = elements_surface_nodes(:);
    elements_surface_nodes = elements_surface_nodes(1:end-(temp2-numel(temp)));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%