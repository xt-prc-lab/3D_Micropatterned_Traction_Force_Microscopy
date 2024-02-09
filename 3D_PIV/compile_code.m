%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                                          Function used to compile the c functions.                                           %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%   Last Revison Date: 09/02/2024                                                                                              %
%   Author: Manuel Gomez Gonzalez                                                                                              %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compile_code()

    % Compile the functions only if any of the mex files does not exist.
    if ((exist('xsubpix_shift_3D_mex', 'file') ~= 3) || (exist('xwindow3D_mex', 'file') ~= 3)) && isunix
        % If the platform is unix, use some compiler optimizations.
        mex COPTIMFLAGS='-O2 -march=native -fomit-frame-pointer -pipe -DNDEBUG' -output xsubpix_shift_3D_mex xsubpix_shift_3D.c
        mex COPTIMFLAGS='-O2 -march=native -fomit-frame-pointer -pipe -DNDEBUG' -output xwindow3D_mex xwindow3D.c

    elseif ((exist('xsubpix_shift_3D_mex', 'file') ~= 3) || (exist('xwindow3D_mex', 'file') ~= 3)) && ~isunix
        % I am not able to check the compiler optimizations for non-unix platforms.

        mex -output xsubpix_shift_3D_mex xsubpix_shift_3D.c
        mex -output xwindow3D_mex xwindow3D.c

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%