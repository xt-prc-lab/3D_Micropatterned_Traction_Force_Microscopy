%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%             Check if the Abaqus model outputs the stiffness matrix. If not, modify it so the matrix is exported.             %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       file [str]: Abaqus job file.                                                                                           %
%       file_save [str]: modified Abaqus job file.                                                                             %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%   Last Revison Date: 22/01/2024                                                                                              %
%   Author: Manuel Gomez Gonzalez                                                                                              %
%                                                                                                                              %
%   References:                                                                                                                %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_stiffness_matrix_job(file, file_save)

    % Open the Abaqus model and import all the text lines.
    fid = fopen(file,'r');
    MyTextFile = textscan(fid,'%s','delimiter','\n','whitespace', '');
    fclose(fid);
    MyTextFile = [MyTextFile{:}];

    % Find the first and last lines of the definition of the upper-surface boundary conditions.
    str = '*Boundary, user';
    ind = find(strncmpi(MyTextFile, str, numel(str)));
    ind0_BC = ind(1);

    str = '**';
    ind = find(strncmpi(MyTextFile((ind0_BC+1):end), str, numel(str)));
    indf_BC = ind0_BC+ind(1);

    MyTextFile((ind0_BC+1):(indf_BC-1)) = [];

    % Find if the abaqus job outputs the stiffness matrix.
    str = '*MATRIX OUTPUT, STIFFNESS, LOAD, FORMAT=MATRIX INPUT';
    ind = find(strncmpi(MyTextFile, str, numel(str)));

    % If the job does not output the stiffness matrix, modify it so it does.
    if isempty(ind)

        MyTextFile{end+1} = '**';
        MyTextFile{end+1} = '** Output Global Stiffness Matrix';
        MyTextFile{end+1} = '**';
        MyTextFile{end+1} = '*Step, name=Global_Stiffness_Matrix';
        MyTextFile{end+1} = '*MATRIX GENERATE, STIFFNESS, LOAD';
        MyTextFile{end+1} = '*MATRIX OUTPUT, STIFFNESS, LOAD, FORMAT=MATRIX INPUT';
        MyTextFile{end+1} = '*End Step';

    end

    % Save the modified abaqus model.
    fid = fopen(file_save, 'w');
    MyTextFile = MyTextFile.';
    fprintf(fid,'%s\n', MyTextFile{:});
    fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%