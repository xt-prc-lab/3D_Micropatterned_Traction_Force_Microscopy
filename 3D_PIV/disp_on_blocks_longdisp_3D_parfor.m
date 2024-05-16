%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%        Function that applies an accurate (slower) PIV to calculate the displacement field of im1 with respect to im2.        %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       im1 [3D matrix]: reference 3D stack image.                                                                             %
%       im2 [3D matrix]: deformed 3D stack image on which the displacement field will be computed.                             %
%       iblocksize [int scalar]: number of rows of the PIV box.                                                                %
%       jblocksize [int scalar]: number of columns of the PIV box.                                                             %
%       kblocksize [int scalar]: number of layers of the PIV box.                                                              %
%       overlap_i [double scalar]: overlap in rows of the PIV boxes.                                                           %
%       overlap_j [double scalar]: overlap in columns of the PIV boxes.                                                        %
%       overlap_k [double scalar]: overlap in layers of the PIV boxes.                                                         %
%       padding_i [int scalar]: zero-padding, in rows, around the PIV boxes.                                                   %
%       padding_j [int scalar]: zero-padding, in columns, around the PIV boxes.                                                %
%       padding_k [int scalar]: zero-padding, in layers, around the PIV boxes.                                                 %
%       blocksizedec [int scalar]: amount by which each side of the interrogation box, in xy, shrinks at each iteration.       %
%       blocksizedec_z [int scalar]: amount by which each side of the interrogation box, in z, shrinks at each iteration.      %
%       window [int scalar]: selects the multiplicative window used on the PIV boxes.                                          %
%                            1 means Hanning window, otherwise means no windowing.                                             %
%       method [int scalar]: Method used for the location of the peak of the 3D correlation.                                   %
%                            1 means the 3D centroid, 2 means Gaussian fit.                                                    %
%       N [int scalar]: number of points around the pixel peak location used to locate it with subpixel accuracy.              %
%       N_lim [int scalar]: maximum number of iterations after decrease of box size. The algorithm will iterate until the      %
%                           number of iterations, after the initial box size decrease, is larger than N_lim or all boxes have  %
%                           converged, whatever happens first.                                                                 %
%       epsilon [double scalar]: maximum acceptable relative error of the iteration for each box to be considered as converged.%
%       N_conver [int scalar]: stop iterating after N_conver iterations without new boxes converging.                          %
%       del [int scalar]: For each box where longdisp is required, apply longdisp to the neighbouring "del" boxes around it.   %
%       del_z [int scalar]: In XY and Z respectively.                                                                          %
%       dynamic_window [int scalar]: when determining the initial window size, we could use different criteria for the         %
%                                    different boxes. 1 if the sizes of the PIV boxes are determined dynamically; 2 if all the %
%                                    boxes larger than blocksize have the maximum initial box size; 0 if all the boxes have    %
%                                    the maximum box size.                                                                     %
%       filt_size_i [int scalar]: Y-size of the Predictor filter, in window lengths, i.e. filter_size = filt_size*window_size. %
%       filt_size_j [int scalar]: X-size of the Predictor filter, in window lengths, i.e. filter_size = filt_size*window_size. %
%       filt_size_k [int scalar]: Z-size of the Predictor filter, in window lengths, i.e. filter_size = filt_size*window_size. %
%       filt_pow [double scalar]: power of the Predictor filter.                                                               %
%       n_th_disp [int scalar]: maximum number of threads used in the parallelization of the PIV.                              %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       x [3D matrix]: x-coordinate of the displacement field.                                                                 %
%       y [3D matrix]: y-coordinate of the displacement field.                                                                 %
%       z [3D matrix]: z-coordinate of the displacement field.                                                                 %
%       dx [3D matrix]: x-component of the displacement field.                                                                 %
%       dy [3D matrix]: y-component of the displacement field.                                                                 %
%       dz [3D matrix]: z-component of the displacement field.                                                                 %
%       pkh [3D matrix]: value of the cross-correlation at the peak.                                                           %
%       conver [vector]: percentage of PIV boxes that have reached convergence at each iteration.                              %
%       iter_time [vector]: time taken to run each iteration.                                                                  %
%       ws_xy0 [3D matrix]: initial window size, in xy, of each PIV box, before any shrinking.                                 %
%       ws_z0 [3D matrix]: initial window size, in z, of each PIV box, before any shrinking.                                   %
%       Convergence_Mask [3D matrix]: for each PIV box, 1 if it has reached convergence, 0 otherwise.                          %
%                                                                                                                              %
%   Last Revison Date: 16/05/2024                                                                                              %
%   Based on the 2D PIV created by Xavier Trepat (03/2008).                                                                    %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                %
%                                                                                                                              %
%   References:                                                                                                                %
%       "A correlation-based continuous window-shift technique to reduce the peak-locking effect in digital PIV image          %
%       evaluation"; Lichuan Gui and Steven T. Wereley; Experiments in Fluids, April 2002, Volume 32, Issue 4, pp 506-517.     %
%       https://doi.org/10.1007/s00348-001-0396-1                                                                              %
%                                                                                                                              %
%       "Particle Image Velocimetry, A Practical Guide"; Markus Raffel, Christian E. Willert, Steve T. Wereley,                %
%       Jürgen Kompenhans; Springer, 2007. https://doi.org/10.1007/978-3-540-72308-0                                           %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x, y, z, dx, dy, dz, pkh, conver, iter_time, ws_xy0, ws_z0, Convergence_Mask] = ...
    disp_on_blocks_longdisp_3D_parfor(im1, im2, iblocksize, jblocksize, kblocksize, overlap_i, overlap_j, overlap_k, ...
                                      padding_i, padding_j, padding_k, blocksizedec, blocksizedec_z, window, method, ...
                                      N, N_lim, epsilon, N_conver, del, del_z, dynamic_window, ...
                                      filt_size_i, filt_size_j, filt_size_k, ...
                                      filt_pow, n_th_disp)

    pad_filter_i = round(filt_size_i/(2*(1-overlap_i)));    % The filter used on the displacement will create border effects. 
    pad_filter_j = round(filt_size_j/(2*(1-overlap_j)));    % In order to avoid them, we create a rim around the matrix, with 
    pad_filter_k = round(filt_size_k/(2*(1-overlap_k)));    % size pad_filter, that is the specular image of the displacement
                                                            % field.

    % Size of the image. Keep in mind that the first index (rows) corresponds to y, the second index (columns) to x and 
    % the third (layers), to z.
    sz = size(im1);

    % We apply the PIV only to the boxes with an intensity higher than the averare intensity of the 3D stack.
    Settings.Fluo_th = mean(im2, 'all', 'omitnan');

    im1 = xExpandMatrix_3D(im1, 1, 1, 1, padding_i, padding_i, padding_j, padding_j, padding_k, padding_k, 0);
    im2 = xExpandMatrix_3D(im2, 1, 1, 1, padding_i, padding_i, padding_j, padding_j, padding_k, padding_k, 0);

    incx = round(jblocksize*(1-overlap_j));             % Define increment in X.
    incy = round(iblocksize*(1-overlap_i));             % Define increment in Y.
    incz = round(kblocksize*(1-overlap_k));             % Define increment in Z.
    % Now with more indices!!!
    
    if incx<1 || incx>jblocksize
        error('Wrong Overlap in Correlation Algorithm')
    end

    if incy<1 || incy>iblocksize
        error('Wrong Overlap in Correlation Algorithm')
    end

    if incz<1 || incz>kblocksize
        error('Wrong Overlap in Correlation Algorithm')
    end

    % x, y and z coordinates of the boxes.
    [x, y, z] = meshgrid((1:incx:sz(2)-jblocksize+1) + jblocksize/2, ...
                         (1:incy:sz(1)-iblocksize+1) + iblocksize/2, ...
                         (1:incz:sz(3)-kblocksize+1) + kblocksize/2);

    dx  = zeros(size(x));               % Displacements in x.
    dy  = zeros(size(x));               % Displacements in y.
    dz  = zeros(size(x));               % Displacements in z.
    pkh = zeros(size(x));               % Height of the xcorr peak.
%     c   = zeros(size(x));               % Cross correlation of the images.

    ws_xy = ones(size(x))*max(iblocksize, jblocksize);      % Size of each window for the longdigsp calculations.
    ws_z  = ones(size(x))*kblocksize;                       % Size of each window for the longdigsp calculations.

    dx_filt = zeros(size(x) + 2*pad_filter_j);          % Displacements in x, filtered.
    dy_filt = zeros(size(x) + 2*pad_filter_i);          % Displacements in y, filtered.
    dz_filt = zeros(size(x) + 2*pad_filter_k);          % Displacements in z, filtered.

    Convergence_Mask = zeros(size(x));                  % 1 for PIV boxes that reach convergence, 0 otherwise.
    e = zeros(size(x));
    inc_err0 = zeros(size(x));

    niterations = 1;
    exit = 0;
    conver = [0];

    % Index of the non-converged windows.
    conv_vect = num2cell(sub2ind(sz, y(:), x(:), z(:))');
    
    % Make a first PIV run with the user-privided window size. It should be fast. Use the displacement of each box to decide 
    % the window size needed to resolve the displacement of this window, i.e. the one-quarter rule: disp < window_size / 4.
    % See "Particle Image Velocimetry, a Practical Guide", by M. Raffel et al.

    % Coordinates of each point in the windows.
    x_sect = 1:jblocksize + 2*padding_j;
    y_sect = 1:iblocksize + 2*padding_i;
    z_sect = 1:kblocksize + 2*padding_k;

    time_iter0 = tic;

    conv_vect_2 = zeros(size(conv_vect), 'logical');

    % Sub-indices of the non-converged windows.
    [ky, kx, kz] = ind2sub(sz, cell2mat(conv_vect));
    kx = kx - jblocksize/2;
    ky = ky - iblocksize/2;
    kz = kz - kblocksize/2;

    % In order to accelerate the 3D PIV, we apply it only to the PIV boxes with values of intensity higher than the average 
    % intensity of the whole stack.
    parfor(ll = 1:numel(conv_vect), n_th_disp)

        % Interrogation windows, cropped to the correct size for each itaration.
        im22 = im2(y_sect + ky(ll)-1, x_sect + kx(ll)-1, z_sect + kz(ll)-1);

        conv_vect_2(ll) = any(im22(:) >= Settings.Fluo_th);

    end

    conv_vect = conv_vect(conv_vect_2);

    % Sub-indices of the non-converged windows.
    [ky, kx, kz] = ind2sub(sz, cell2mat(conv_vect));
    kx = kx - jblocksize/2;
    ky = ky - iblocksize/2;
    kz = kz - kblocksize/2;

    jj = (kx+incx-1)/incx;
    ii = (ky+incy-1)/incy;
    kk = (kz+incz-1)/incz;

    ind = sub2ind(size( x ), ii, jj, kk);

    dx_ind = zeros(size(conv_vect));
    dy_ind = zeros(size(conv_vect));
    dz_ind = zeros(size(conv_vect));

    % First PIV iteration, with the user provided box size.
    parfor(ll = 1:numel(conv_vect), n_th_disp)

        % Interrogation windows, cropped to the correct size for each iteration.
        im11 = im1(y_sect + ky(ll)-1, x_sect + kx(ll)-1, z_sect + kz(ll)-1);
        im22 = im2(y_sect + ky(ll)-1, x_sect + kx(ll)-1, z_sect + kz(ll)-1);

        % Normalize the images.
        im11 = (im11 - mean(im11, 'all'))/std_flat(im11);
        im22 = (im22 - mean(im22, 'all'))/std_flat(im22);

        % Multiply each block by a window.
        im11 = xwindow3D_mex(im11, window);
        im22 = xwindow3D_mex(im22, window);

        c = real(xcorrf3(im11, im22, 'no'));

        [xc, yc, zc] = x_cntr_3D(c, N, method);

        dx_ind(ll) = real(xc);          % Accumulated displacement across iterations.
        dy_ind(ll) = real(yc);
        dz_ind(ll) = real(zc);

    end

    dx(ind) = dx_ind;
    dy(ind) = dy_ind;
    dz(ind) = dz_ind;

    iter_time(1) = toc(time_iter0);
    time_iter0 = datestr(iter_time(1)/(60*60*24), 'HH:MM:SS')

    [n, m, p] = size(ws_xy);

    % Calculate the initial PIV box size needed, so that it is larger than four times the box displacement.
    parfor(ii=1:n, n_th_disp)
        for jj=1:m
            for kk=1:p

                % If the displacement of a box (or neighboutring boxes) is larger than one fourth of the user provided 
                % window size, it cannot be resolved with the user provided box size. In that case, enlarge the initial box size 
                % (longdisp algorithm).
                tmp = max(sqrt(...
                    dx(max(1, ii-del):min(n, ii+del), max(1, jj-del):min(m, jj+del), max(1, kk-del_z):min(p, kk+del_z)).^2 + ...
                    dy(max(1, ii-del):min(n, ii+del), max(1, jj-del):min(m, jj+del), max(1, kk-del_z):min(p, kk+del_z)).^2), ...
                    [], 'all', 'omitnan');

                ws_xy(ii, jj, kk) = ws_xy(ii, jj, kk) + ceil(max(0, 4*tmp - ws_xy(ii, jj, kk))/(2*blocksizedec))*2*blocksizedec;

                tmp = max(abs(...
                    dz(max(1, ii-del):min(n, ii+del), max(1, jj-del):min(m, jj+del), max(1, kk-del_z):min(p, kk+del_z))), ...
                    [], 'all', 'omitnan');

                ws_z(ii, jj, kk) = ws_z(ii, jj, kk) + ceil(max(0, 4*tmp - ws_z(ii, jj, kk))/(2*blocksizedec_z))*2*blocksizedec_z;

            end
        end
    end

    if dynamic_window==0
        % All the PIV boxes start with the maximum size calculated for the whole domain.
        ws_xy(:) = max(ws_xy, [], 'all');
        ws_z(:) = max(ws_z, [], 'all');
    elseif dynamic_window==2
        % All the boxes that neen enlarging will start with the maximum size. 
        % The other boxes will start with the user provided window size.
        ws_xy(ws_xy>jblocksize) = max(ws_xy(:));
        ws_z(ws_z>kblocksize) = max(ws_z(:));
    end

    ws_xy0 = ws_xy;
    ws_z0 = ws_z;

    ws_xy_max = max(ws_xy(:));
    ws_z_max = max(ws_z(:));

    pad_dec_i = max(padding_i, (ws_xy_max - iblocksize)/2);         % Padding for the increasing of the edges in the longdisp.
    pad_dec_j = max(padding_j, (ws_xy_max - jblocksize)/2);
    pad_dec_k = max(padding_k, (ws_z_max - kblocksize)/2);

    im1 = xExpandMatrix_3D(im1(1+padding_i:end-padding_i, 1+padding_j:end-padding_j, 1+padding_k:end-padding_k), 1, 1, 1, ...
                           pad_dec_i, pad_dec_i, pad_dec_j, pad_dec_j, pad_dec_k, pad_dec_k, 0);
    im2 = xExpandMatrix_3D(im2(1+padding_i:end-padding_i, 1+padding_j:end-padding_j, 1+padding_k:end-padding_k), 1, 1, 1, ...
                           pad_dec_i, pad_dec_i, pad_dec_j, pad_dec_j, pad_dec_k, pad_dec_k, 0);

    % Coordinates of each point in the largest windows.
    x_sect = 1:(jblocksize + 2*pad_dec_j);
    y_sect = 1:(iblocksize + 2*pad_dec_i);
    z_sect = 1:(kblocksize + 2*pad_dec_k);

    % Perform the iterations with box size decrease.
    while (exit~= 1)
 
        niterations
        time_iter = tic;

        % Sub-indices of the non-converged windows.
        [ky, kx, kz] = ind2sub(sz, cell2mat(conv_vect));
        kx = kx - jblocksize/2;
        ky = ky - iblocksize/2;
        kz = kz - kblocksize/2;

        jj = (kx+incx-1)/incx;
        ii = (ky+incy-1)/incy;
        kk = (kz+incz-1)/incz;

        ind = sub2ind(size(x), ii, jj, kk);
        ind_2 = sub2ind([size(x, 1) + 2*pad_filter_i, size(x, 2) + 2*pad_filter_j, size(x, 3) + 2*pad_filter_k], ...
                        ii + pad_filter_i, jj + pad_filter_j, kk + pad_filter_k);

        % Re-define variables acceptable for the parfor loop.
        conv_vect_2 = zeros(size(conv_vect));

        ws_xy_2 = ws_xy(ind);
        ws_z_2 = ws_z(ind);

        inc_err0_2 = inc_err0(ind);
        Convergence_Mask_2 = Convergence_Mask(ind);
        e_2 = e(ind);

        dx_2 = dx(ind);
        dy_2 = dy(ind);
        dz_2 = dz(ind);
        pkh_2 = pkh(ind);

        DX = dx_filt(ind_2);        % This is the subpixel amount of shifting until getting peak locking conditions
        DY = dy_filt(ind_2);
        DZ = dz_filt(ind_2);

        x_2 = x(ind);
        y_2 = y(ind);
        z_2 = z(ind);

        x_wind_dec = (x_sect(end) - ws_xy_2)/2.;        % Reduction of the cropped window.
        y_wind_dec = (y_sect(end) - ws_xy_2)/2.;
        z_wind_dec = (z_sect(end) - ws_z_2 )/2.;

        parfor(ll=1:numel(conv_vect), n_th_disp)
%         for ll=1:numel(conv_vect)
% 
%             fprintf('.');
% 
%             disp(['niterations = ', num2str(niterations), '; ii = ', num2str(ii), '; jj = ', num2str(jj)])

            % Coordinates of the window where we will perform the cross-correlation. It is at the center of x_sect, y_sect.
            x_wind = x_sect(1+x_wind_dec(ll):x_sect(end)-x_wind_dec(ll));
            y_wind = y_sect(1+y_wind_dec(ll):y_sect(end)-y_wind_dec(ll));
            z_wind = z_sect(1+z_wind_dec(ll):z_sect(end)-z_wind_dec(ll));

            x_wind_shift = x_sect(1+max(0, x_wind_dec(ll)-1):x_sect(end)-max(0, x_wind_dec(ll)-1));
            y_wind_shift = y_sect(1+max(0, y_wind_dec(ll)-1):y_sect(end)-max(0, y_wind_dec(ll)-1));
            z_wind_shift = z_sect(1+max(0, z_wind_dec(ll)-1):z_sect(end)-max(0, z_wind_dec(ll)-1));

            x_wind_shift = x_wind_shift((x_wind_shift-round( DX(ll)/2)>=1)&(x_wind_shift-round( DX(ll)/2)<=x_sect(end))& ...
                                        (x_wind_shift-round(-DX(ll)/2)>=1)&(x_wind_shift-round(-DX(ll)/2)<=x_sect(end)));
            y_wind_shift = y_wind_shift((y_wind_shift-round( DY(ll)/2)>=1)&(y_wind_shift-round( DY(ll)/2)<=y_sect(end))& ...
                                        (y_wind_shift-round(-DY(ll)/2)>=1)&(y_wind_shift-round(-DY(ll)/2)<=y_sect(end)));
            z_wind_shift = z_wind_shift((z_wind_shift-round( DZ(ll)/2)>=1)&(z_wind_shift-round( DZ(ll)/2)<=z_sect(end))& ...
                                        (z_wind_shift-round(-DZ(ll)/2)>=1)&(z_wind_shift-round(-DZ(ll)/2)<=z_sect(end)));

            % Interrogation windows, cropped to the correct size for each iteration.
            if numel(x_wind_shift) && numel(y_wind_shift) && numel(z_wind_shift)

                im11 = im1(y_sect+y_2(ll)-1-ws_xy_max/2+pad_dec_i, ...
                           x_sect+x_2(ll)-1-ws_xy_max/2+pad_dec_j, ...
                           z_sect+z_2(ll)-1-ws_z_max/2+pad_dec_k);
                im22 = im2(y_sect+y_2(ll)-1-ws_xy_max/2+pad_dec_i, ...
                           x_sect+x_2(ll)-1-ws_xy_max/2+pad_dec_j, ...
                           z_sect+z_2(ll)-1-ws_z_max/2+pad_dec_k);

                % Subpixel shift of the image.
                im11(y_wind_shift, x_wind_shift, z_wind_shift) = xsubpix_shift_3D_mex( ...
                    im11(y_wind_shift-round( DY(ll)/2), x_wind_shift-round( DX(ll)/2), z_wind_shift-round( DZ(ll)/2)), ...
                     DX(ll)/2-round( DX(ll)/2),  DY(ll)/2-round( DY(ll)/2),  DZ(ll)/2-round( DZ(ll)/2));
                im22(y_wind_shift, x_wind_shift, z_wind_shift) = xsubpix_shift_3D_mex( ...
                    im22(y_wind_shift-round(-DY(ll)/2), x_wind_shift-round(-DX(ll)/2), z_wind_shift-round(-DZ(ll)/2)), ...
                    -DX(ll)/2-round(-DX(ll)/2), -DY(ll)/2-round(-DY(ll)/2), -DZ(ll)/2-round(-DZ(ll)/2));

                % Test the convergence of the previous step. If we test the convergence of the current step, we have to
                % perform again xsubpix_shift_3D_mex(  ), which is a very expensive function. This way, reduce the calls
                % to the function by half.
                I = sum((im11(y_wind, x_wind, z_wind) - im22(y_wind, x_wind, z_wind)).^2, 'all');

                err = sqrt(I/((x_wind(end)-x_wind(1)+1)*(y_wind(end)-y_wind(1)+1)*(z_wind(end)-z_wind(1)+1)))/...
                      ((std_flat(im11(y_wind, x_wind, z_wind)) + std_flat(im22(y_wind, x_wind, z_wind)))/2);

                % Error parameters of the first iteration (it is used after the second iteration).
                if niterations==2

                    I0 = sum((im11 - im22).^2, 'all');

                    err0 = sqrt(I0/((x_sect(end)-x_sect(1)+1)*(y_sect(end)-y_sect(1)+1)*(z_sect(end)-z_sect(1)+1)))/...
                           ((std_flat(im11) + std_flat(im22))/2);

                    inc_err0_2(ll) = err0 - err;

                end

                if (abs((err - e_2(ll))/inc_err0_2(ll)) <= epsilon) && (ws_xy_2(ll)==jblocksize) && (ws_z_2(ll)==kblocksize)

                    Convergence_Mask_2(ll) = 1;
                    conv_vect_2(ll) = 1;
                    dx_2(ll) = DX(ll);
                    dy_2(ll) = DY(ll);
                    dz_2(ll) = DZ(ll);

                else

                    % Normalize the images.
                    im11(y_wind, x_wind, z_wind) = (im11(y_wind, x_wind, z_wind) - ...
                                    mean(reshape(im11(y_wind, x_wind, z_wind), [], 1)))/std_flat(im11(y_wind, x_wind, z_wind));
                    im22(y_wind, x_wind, z_wind) = (im22( y_wind, x_wind, z_wind) - ...
                                    mean(reshape(im22(y_wind, x_wind, z_wind), [], 1)))/std_flat(im22(y_wind, x_wind, z_wind));

                    % Multiply each block by a window.
                    im11(y_wind, x_wind, z_wind) = xwindow3D_mex(im11(y_wind, x_wind, z_wind), window);
                    im22(y_wind, x_wind, z_wind) = xwindow3D_mex(im22(y_wind, x_wind, z_wind), window);

                    c = real(xcorrf3(im11(y_wind, x_wind, z_wind), im22(y_wind, x_wind, z_wind), 'no'));

                    [xc, yc, zc] = x_cntr_3D(c, N, method);

                    dx_2(ll) = real(DX(ll) + xc);       % Accumulated displacement across iterations.
                    dy_2(ll) = real(DY(ll) + yc);
                    dz_2(ll) = real(DZ(ll) + zc);
                    pkh_2(ll) = max(abs(c(:)));         % Store peak height.

                end

                e_2(ll) = err;
        
            else

                if niterations==2
                    inc_err0_2(ll) = NaN;
                end

                Convergence_Mask_2(ll) = 1;
                conv_vect_2(ll) = 1;
                dx_2(ll) = NaN;
                dy_2(ll) = NaN;
                dz_2(ll) = NaN;

                e_2(ll) = NaN;

            end

            % Determine whether the window need to be further shrunk.
            ws_xy_2(ll) = max(jblocksize, ws_xy_2(ll) - 2*blocksizedec);
            ws_z_2(ll) = max(kblocksize, ws_z_2(ll) - 2*blocksizedec_z);

        end

        % Re-assign the intermediate variables used in the parfor loop.
        ws_xy(ind) = ws_xy_2;
        ws_z(ind) = ws_z_2;

        inc_err0(ind) = inc_err0_2;
        Convergence_Mask(ind) = Convergence_Mask_2;
        e(ind) = e_2;

        dx(ind) = dx_2;
        dy(ind) = dy_2;
        dz(ind) = dz_2;
        pkh(ind) = pkh_2;

        conv_vect(conv_vect_2==1) = [];

        time_iter1 = datestr(toc(time_iter)/(60*60*24), 'HH:MM:SS')

        % Add the specular image of the borders of the displacement field before filtering.
        dx_filt(1+pad_filter_i:end-pad_filter_i, 1+pad_filter_j:end-pad_filter_j, 1+pad_filter_k:end-pad_filter_k) = dx;
        dy_filt(1+pad_filter_i:end-pad_filter_i, 1+pad_filter_j:end-pad_filter_j, 1+pad_filter_k:end-pad_filter_k) = dy;
        dz_filt(1+pad_filter_i:end-pad_filter_i, 1+pad_filter_j:end-pad_filter_j, 1+pad_filter_k:end-pad_filter_k) = dz;

        dx_filt(1:pad_filter_i, :, :) = dx_filt(1+2*pad_filter_i:-1:2+pad_filter_i, :, :);
        dy_filt(1:pad_filter_i, :, :) = dy_filt(1+2*pad_filter_i:-1:2+pad_filter_i, :, :);
        dz_filt(1:pad_filter_i, :, :) = dz_filt(1+2*pad_filter_i:-1:2+pad_filter_i, :, :);

        dx_filt(end-pad_filter_i+1:end, :, :) = dx_filt(end-pad_filter_i-1:-1:end-2*pad_filter_i, :, :);
        dy_filt(end-pad_filter_i+1:end, :, :) = dy_filt(end-pad_filter_i-1:-1:end-2*pad_filter_i, :, :);
        dz_filt(end-pad_filter_i+1:end, :, :) = dz_filt(end-pad_filter_i-1:-1:end-2*pad_filter_i, :, :);

        dx_filt(:, 1:pad_filter_j, :) = dx_filt(:, 1+2*pad_filter_j:-1:2+pad_filter_j, :);
        dy_filt(:, 1:pad_filter_j, :) = dy_filt(:, 1+2*pad_filter_j:-1:2+pad_filter_j, :);
        dz_filt(:, 1:pad_filter_j, :) = dz_filt(:, 1+2*pad_filter_j:-1:2+pad_filter_j, :);

        dx_filt(:, end-pad_filter_j+1:end, :) = dx_filt(:, end-pad_filter_j-1:-1:end-2*pad_filter_j, :);
        dy_filt(:, end-pad_filter_j+1:end, :) = dy_filt(:, end-pad_filter_j-1:-1:end-2*pad_filter_j, :);
        dz_filt(:, end-pad_filter_j+1:end, :) = dz_filt(:, end-pad_filter_j-1:-1:end-2*pad_filter_j, :);

        dx_filt(:, :, 1:pad_filter_k) = dx_filt(:, :, 1+2*pad_filter_k:-1:2+pad_filter_k);
        dy_filt(:, :, 1:pad_filter_k) = dy_filt(:, :, 1+2*pad_filter_k:-1:2+pad_filter_k);
        dz_filt(:, :, 1:pad_filter_k) = dz_filt(:, :, 1+2*pad_filter_k:-1:2+pad_filter_k);

        dx_filt(:, :, end-pad_filter_k+1:end) = dx_filt(:, :, end-pad_filter_k-1:-1:end-2*pad_filter_k);
        dy_filt(:, :, end-pad_filter_k+1:end) = dy_filt(:, :, end-pad_filter_k-1:-1:end-2*pad_filter_k);
        dz_filt(:, :, end-pad_filter_k+1:end) = dz_filt(:, :, end-pad_filter_k-1:-1:end-2*pad_filter_k);

        % Filtering.
        dx_filt = inpaint_nans3(dx_filt);
        dy_filt = inpaint_nans3(dy_filt);
        dz_filt = inpaint_nans3(dz_filt);

        dx_filt = filter_3D_Disp(dx_filt, ...
                            [filt_size_i/(1 - overlap_i), filt_size_j/(1 - overlap_j), filt_size_k/(1 - overlap_k)], filt_pow);
        dy_filt = filter_3D_Disp( dy_filt, ...
                            [filt_size_i/(1 - overlap_i), filt_size_j/(1 - overlap_j), filt_size_k/(1 - overlap_k)], filt_pow);
        dz_filt = filter_3D_Disp( dz_filt, ...
                            [filt_size_i/(1 - overlap_i), filt_size_j/(1 - overlap_j), filt_size_k/(1 - overlap_k)], filt_pow);

        % Outlier removal.
        dx_filt = removeOutliers(dx_filt);
        dy_filt = removeOutliers(dy_filt);
        dz_filt = removeOutliers(dz_filt);

        dx_filt = dx_filt{1};
        dy_filt = dy_filt{1};
        dz_filt = dz_filt{1};

        % Deform images and check convergence and decide if the window needs to be shrunk or not and for which magnitude.
        fprintf('\n');

        if niterations >= 2

            conver(niterations) = sum(Convergence_Mask(:))/numel(Convergence_Mask)*100;
            disp(['At iteration number ', num2str(niterations-1), ',   ', num2str(conver(end)), ...
                  ' % of points have reached convergence'])

            exit = (niterations >= max([pad_dec_i/blocksizedec, pad_dec_j/blocksizedec, pad_dec_k/blocksizedec_z]) + N_lim) ...
                   || (conver(end)==100);

            if niterations >= N_conver+1
                exit = exit || (sum(abs(diff(conver(min(numel(conver), end-N_conver+1):end))))==0);
            end

        end

        iter_time(niterations + 1) = toc(time_iter);

        niterations = niterations + 1;

    end

    fprintf('\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%