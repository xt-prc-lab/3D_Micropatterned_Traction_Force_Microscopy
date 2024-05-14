%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%           Function used to repeat a matrix an integer number of times in the i, j and k directions. It also adds a           %
%                                        specified padding in the i, j and k directions.                                       %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       A [matrix]: matrix to expand.                                                                                          %
%       DimRow [int]: number of repetitions along the rows.                                                                    %
%       DimCol [int]: number of repetitions along the columns.                                                                 %
%       DimZ   [int]: number of repetitions along the layers.                                                                  %
%       padYTop [int]: padding added at the top of the rows.                                                                   %
%       padYBottom [int]: padding added at the bottom of the rows.                                                             %
%       padXLeft [int]: padding added at the left of the columns.                                                              %
%       padXRight [int]: padding added at the right of the columns.                                                            %
%       padZTop [int]: padding added at the top of the layers.                                                                 %
%       padZBottom [int]: padding added at the bottom of the layers.                                                           %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       E [3D matrix]: expanded matrix.                                                                                        %
%                                                                                                                              %
%   Last Revison Date: 26/02/2024                                                                                              %
%   Based on the code created by Prof. Xavier Trepat for 2D expansion (2008).                                                  %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                %
%                                                                                                                              %
%   References:                                                                                                                %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E = xExpandMatrix_3D(A,DimRow,DimCol,DimZ,padYTop,padYBottom,padXLeft,padXRight,padZTop,padZBottom,padding)

    % Dimension of the final matrix before padding:
    FinalDimRow = DimRow*size(A,1);
    FinalDimCol = DimCol*size(A,2);
    FinalDimZ   = DimZ*  size(A,3);

    % Initialize the final array.
    E = zeros(FinalDimRow+padYTop+padYBottom, FinalDimCol+padXLeft+padXRight, FinalDimZ+padZTop+padZBottom);
    E(:) = padding;

    % Expand the matrix:
    if (DimRow==1) && (DimCol==1) && (DimZ==1)
        E((1:FinalDimRow)+padYTop, (1:FinalDimCol)+padXLeft, (1:FinalDimZ)+padZTop ) = A;

    else
        E((1:FinalDimRow)+padYTop, (1:FinalDimCol)+padXLeft, (1:FinalDimZ)+padZTop) = ...
                    repmat(A, [DimRow, DimCol, DimZ]);

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%