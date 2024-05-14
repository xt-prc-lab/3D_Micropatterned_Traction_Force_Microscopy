%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%               Function that applies a fast PIV to calculate the displacement field of im1 with respect to im2.               %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       im1 [3D matrix]: reference 3D stack image.                                                                             %
%       im2 [3D matrix]: deformed 3D stack image on which the displacement field will be computed.                             %
%       iblocksize [int]: number of rows of the PIV box.                                                                       %
%       jblocksize [int]: number of columns of the PIV box.                                                                    %
%       kblocksize [int]: number of layers of the PIV box.                                                                     %
%       overlap_i [scalar]: overlap in rows of the PIV boxes.                                                                  %
%       overlap_j [scalar]: overlap in columns of the PIV boxes.                                                               %
%       overlap_k [scalar]: overlap in layers of the PIV boxes.                                                                %
%       padding_i [int]: zero-padding, in rows, around the PIV boxes.                                                          %
%       padding_j [int]: zero-padding, in columns, around the PIV boxes.                                                       %
%       padding_k [int]: zero-padding, in layers, around the PIV boxes.                                                        %
%       window [int]: selects the multiplicative window used on the PIV boxes.                                                 %
%                     1 means Hanning window, otherwise means no windowing.                                                    %
%       method [int]: Method used for the location of the peak of the 3D correlation.                                          %
%                     1 means the 3D centroid, 2 means Gaussian fit.                                                           %
%       N [int]: number of points around the pixel peak location used to locate it with subpixel accuracy.                     %
%       Ni [int]: minimum number of iterations. The algorithm will iterate until the number of iterations is larger than Ni    %
%                 and the desired accuracy is reached.                                                                         %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       x [3D matrix]: x-coordinate of the displacement field.                                                                 %
%       y [3D matrix]: y-coordinate of the displacement field.                                                                 %
%       z [3D matrix]: z-coordinate of the displacement field.                                                                 %
%       dx [3D matrix]: x-component of the displacement field.                                                                 %
%       dy [3D matrix]: y-component of the displacement field.                                                                 %
%       dz [3D matrix]: z-component of the displacement field.                                                                 %
%       pkh [scalar]: value of the cross-correlation at the peak.                                                              %
%                                                                                                                              %
%   Last Revison Date: 04/03/2024                                                                                              %
%   Based on the 2D PIV created by Xavier Trepat (03/2008).                                                                    %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                %
%                                                                                                                              %
%   References:                                                                                                                %
%       "A correlation-based continuous window-shift technique to reduce the peak-locking effect in digital PIV image          %
%       evaluation"; Lichuan Gui and Steven T. Wereley; Experiments in Fluids, April 2002, Volume 32, Issue 4, pp 506-517.     %
%       https://doi.org/10.1007/s00348-001-0396-1                                                                              %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, y, z, dx, dy, dz, pkh] = ...
                    disp_on_blocks_fast_3D(im1, im2, iblocksize, jblocksize, kblocksize, overlap_i, overlap_j, overlap_k, ...
                                           padding_i, padding_j, padding_k, window, method, N, Ni)

    sz = size(im1);             % Size of the image.

    % Check the number of PIV boxes along each direction. If it is negative in any direction, don't run the PIV.
    yaux = (sz(1)/iblocksize-1)*(1/(1-overlap_i)) + 1;
    xaux = (sz(2)/jblocksize-1)*(1/(1-overlap_j)) + 1;
    zaux = (sz(3)/kblocksize-1)*(1/(1-overlap_k)) + 1;

    if (xaux<0) || (yaux<0) || (zaux<0)

        x = NaN;
        y = NaN;
        z = NaN;
        dx = NaN;
        dy = NaN;
        dz = NaN;
        pkh = NaN;

    else

        % Padd the images, if needed.
        im1 = xExpandMatrix_3D(im1, 1, 1, 1, padding_i, padding_i, padding_j, padding_j, padding_k, padding_k, 0); 
        im2 = xExpandMatrix_3D(im2, 1, 1, 1, padding_i, padding_i, padding_j, padding_j, padding_k, padding_k, 0);

        inci = round(iblocksize*(1-overlap_i));             % Define increment in Y.
        incj = round(jblocksize*(1-overlap_j));             % Define increment in X.
        inck = round(kblocksize*(1-overlap_k));             % Define increment in Z.

        if (inci<1 || inci>iblocksize) || (incj<1 || incj>jblocksize) || (inck<1 || inck>kblocksize)
            error('Wrong Overlap in Correlation Algorithm')
        end

        % x, y and z coordinates of the PIV boxes.
        [x, y, z] = meshgrid((1:incj:sz(2)-jblocksize+1) + jblocksize/2, ...
                             (1:inci:sz(1)-iblocksize+1) + iblocksize/2, ...
                             (1:inck:sz(3)-kblocksize+1) + kblocksize/2);

        dx  = zeros(size(x));               % Displacements in x.
        dy  = zeros(size(x));               % Displacements in y.
        dz  = zeros(size(x));               % Displacements in z.
        pkh = zeros(size(x));               % Height of the xcorr peak.
        c   = zeros(size(x));               % Cross correlation of the images.

        % Loop along the rows.
        for ki = 1:inci:sz(1)-iblocksize+1

            fprintf('.');

            % Loop along the columns.
            for kj = 1:incj:sz(2)-jblocksize+1

                % Loop along the layers.
                for kk = 1:inck:sz(3)-kblocksize+1

                    % Initialize iterative process.
                    niterations = 0;
                    DX = inf;
                    DY = inf;
                    DZ = inf;

                    % Repeat until the number of iterations is larger than Ni AND the displacement calculated in two successive
                    % iterations is less than 0.02 pixels in every direction.
                    while niterations<=Ni && (abs(dx((ki+inci-1)/inci, (kj+incj-1)/incj, (kk+inck-1)/inck) - DX )>0.02 || ...
                                              abs(dy((ki+inci-1)/inci, (kj+incj-1)/incj, (kk+inck-1)/inck) - DY )>0.02 || ...
                                              abs(dz((ki+inci-1)/inci, (kj+incj-1)/incj, (kk+inck-1)/inck) - DZ )>0.02)

                        niterations = niterations+1;
                        DX = dx((ki+inci-1)/inci, (kj+incj-1)/incj, (kk+inck-1)/inck);
                        DY = dy((ki+inci-1)/inci, (kj+incj-1)/incj, (kk+inck-1)/inck);
                        DZ = dz((ki+inci-1)/inci, (kj+incj-1)/incj, (kk+inck-1)/inck);

                        im11 = im1(ki:ki+iblocksize+2*padding_i-1, ...
                                   kj:kj+jblocksize+2*padding_j-1, ...
                                   kk:kk+kblocksize+2*padding_k-1);
                        im22 = im2(ki:ki+iblocksize+2*padding_i-1, ...
                                   kj:kj+jblocksize+2*padding_j-1, ...
                                   kk:kk+kblocksize+2*padding_k-1);

                        % Subpixel-shift the images. Skip the first iteration.
                        % We shift an subimage of the target size plus the padding, so if the shift is large, we will have 
                        % actual data at the sides rather than zeros.
                        if DX || DY
                            im11 = xsubpix_shift_3D_mex(im11, DX/2, DY/2, DZ/2);
                            im22 = xsubpix_shift_3D_mex(im22,-DX/2,-DY/2,-DZ/2);
                        end

                        % The cross-correlation is performed with a subimage of the target size, without padding.
                        x_sect = (1:iblocksize)+padding_i;
                        y_sect = (1:jblocksize)+padding_j;
                        z_sect = (1:kblocksize)+padding_k;

                        % Subtract the mean
                        im11(x_sect, y_sect, z_sect) = im11(x_sect, y_sect, z_sect) - ...
                                                       mean(im11(x_sect, y_sect, z_sect), 'all', 'omitnan');
                        im22(x_sect, y_sect, z_sect) = im22(x_sect, y_sect, z_sect) - ...
                                                       mean(im22(x_sect, y_sect, z_sect), 'all', 'omitnan');            

                        % Multiply each block by a window.
                        im11(x_sect, y_sect, z_sect) = xwindow3D_mex(im11(x_sect, y_sect, z_sect), window);
                        im22(x_sect, y_sect, z_sect) = xwindow3D_mex(im22(x_sect, y_sect, z_sect), window);

                        % Cross-correlation.
                        c = real(xcorrf3(im11(x_sect, y_sect, z_sect), im22(x_sect, y_sect, z_sect), 'no'));

                        % Find the location of the cross-correlation's peak.
                        [xc, yc, zc] = x_cntr_3D(c, N, method);
                        dx((ki+inci-1)/inci, (kj+incj-1)/incj, (kk+inck-1)/inck) = real(DX+xc);
                        dy((ki+inci-1)/inci, (kj+incj-1)/incj, (kk+inck-1)/inck) = real(DY+yc);
                        dz((ki+inci-1)/inci, (kj+incj-1)/incj, (kk+inck-1)/inck) = real(DZ+zc); 

                    end

                end

                % Height of the cross-correlation's peak.
                pkh((ki+inci-1)/inci, (kj+incj-1)/incj, (kk+inck-1)/inck) = max(abs(c(:)));     
 
            end

        end

    end

    fprintf('\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%