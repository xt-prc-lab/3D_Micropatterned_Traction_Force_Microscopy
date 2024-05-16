%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                                 Function used to apply a 3D filter to the displacement field.                                %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       u0 [matrix]: 3D field to filter.                                                                                       %
%       filter_length [matrix]: size of the applied filter.                                                                    %
%       z [scalar]: power of the applied filter.                                                                               %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       u [matrix]: filtered 2D field.                                                                                         %
%                                                                                                                              %
%   Last Revison Date: 22/01/2024                                                                                              %
%   Author: Manuel Gomez Gonzalez                                                                                              %
%                                                                                                                              %
%   References:                                                                                                                %
%       "Effect of predictor corrector filtering on the stability and spatial resolution of iterative PIV interrogation";      %
%       F. F. J. Schrijer and F. Scarano; Experiments in Fluids, November 2008, Volume 45, Issue 5, pp 927–941.                %
%       https://doi.org/10.1007/s00348-008-0511-7                                                                              %
%                                                                                                                              %
%       "A Fast Iterative Digital Volume Correlation Algorithm for Large Deformations"; E. Bar-Kochba, J. Toyjanova,           %
%       E. Andrews, K.-S. Kim and C. Franck; Experimental Mechanics, January 2015, Volume 55, Issue 1, pp 261–274.             %
%       https://doi.org/10.1007/s11340-014-9874-2                                                                              %
%                                                                                                                              %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------------------------------------------------------------------------------------------------%
%                                            Check the inputs and apply the filter.                                            %
%------------------------------------------------------------------------------------------------------------------------------%

function u = filter_3D_Disp(u0, filter_length, z)

    % Check inputs and assign default values.
    if nargin==0
        error('Wrong number of inputs');

    elseif nargin==1
        filter_length = [1, 1, 1];
        z = 0.0075;

    elseif nargin==2
        z = 0.0075;

    end

    % If only one filter length is provided, use if for the three dimmensions.
    if isscalar(filter_length)
        filter_length(2, 3) = filter_length;
    end

    % Keep only positive filter lengths.
    if any(filter_length<0)
        filter_length = abs(filter_length);

    elseif any(filter_length==0)
        filter_length(filter_length==0) = 1;

    end

    if z==0
        % Impulse filter. Don't filter the 2D field.
        u = double(u0);

    else
        % Define the kernel of the filter, based on the filter size and power.
        kern = filter_kernel(filter_length, z);

        % Convolve the 2D field with the filter kernel.
        u = double(convn(u0, kern, 'same'));

    end

end

%------------------------------------------------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------------------------------------------%
%                                           Generate the kernel used for the filter.                                           %
%------------------------------------------------------------------------------------------------------------------------------%

function kern = filter_kernel(filter_length, z)

    % We create a symmetric filter. The half-length of the filter must be an integer.
    l = 2.*round(filter_length./2);

    [X, Y, Z] = meshgrid(-l(1)/2:l(1)/2, -l(2)/2:l(2)/2, -l(3)/2:l(3)/2);

    X = abs(X);
    Y = abs(Y);
    Z = abs(Z);

    octant_1 = (X >= Y)&(X >= Z);
    octant_2 = (X <  Y)&(Y >= Z);
    octant_3 = (Z >  Y)&(Z >  X);

    % Kernel.
    kern = ((l(1)/2)^z - X.^z).*octant_1 + ((l(2)/2)^z - Y.^z).*octant_2 + ((l(3)/2)^z - Z.^z).*octant_3;
    kern = kern/sum(kern(:));

end

%------------------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%