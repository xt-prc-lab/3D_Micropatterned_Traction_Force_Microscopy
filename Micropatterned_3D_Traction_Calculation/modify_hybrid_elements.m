%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                  Check if the Abaqus model uses Hybrid elements. If so, modify the job file to NOT use them.                 %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       file [str]: Abaqus job file.                                                                                           %
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

function modify_hybrid_elements(file)

    % Open the Abaqus model and import all the text lines.
    fid = fopen(file,'r');
    MyTextFile = textscan(fid,'%s','delimiter','\n','whitespace', '');
    fclose(fid);
    MyTextFile = [MyTextFile{:}];

    % Find the line of the definition of the elements of the part model.
    str = '*Element, type=';
    ind = find(strncmpi(MyTextFile, str, numel(str)));
    ind_elements_part_def = ind(1);

    % Check if the element type is hybrid. If so, change it so hybrid elements are not used.
    if MyTextFile{ind_elements_part_def}(end)=='H'
        MyTextFile{ind_elements_part_def} = MyTextFile{ind_elements_part_def}(1:end-1);
    end

    % Save the modified abaqus model.
    fid = fopen(file, 'w');
    MyTextFile = MyTextFile.';
    fprintf(fid,'%s\n', MyTextFile{:});
    fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%