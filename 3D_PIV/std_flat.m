%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                        Function used to calculate the standard deviation of all the terms of an array.                       %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       a [n-D matrix]: matrix whose elements we use to calculate the standard deviation.                                      %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       st [double scalar]: standard deviation of all the terms of the matrix a.                                               %
%                                                                                                                              %
%   Last Revison Date: 15/04/2024                                                                                              %
%   Author: Manuel Gomez Gonzalez                                                                                              %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function st = std_flat( a )

    st = sqrt(mean((a - mean(a, 'all')).^2, 'all'));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
