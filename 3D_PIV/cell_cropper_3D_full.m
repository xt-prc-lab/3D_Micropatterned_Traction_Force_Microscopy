function [Settings, File] = cell_cropper_3D_full(Settings, File)

    %% CELL_CROPPER crops all images from the stack i.e. the gel fluorescence image,
    % the reference image and phase contrast image (BF). Images are cropped so that the
    % shift between the GFI and the corresponding RI is smaller than 0.5 pixels.
    % Usually a stack of RIs is acquired to account for unfocussing during the
    % time lapse experiments. Each GFI is cross-correlated with all RIs to
    % determine the RI whose focussing plane is closer to that of the GFI. Then the
    % shift between GFI and RI is calculated and both images are croped and
    % stored in the folder CROPEDDATA.

    % The parameters in the structure SETTINGS to select the ROI of the GFI
    % and RI that are crosscorrelated. The program also allows to
    % croscorrelated several ROIs of the RI and GFI to obain a more accurate
    % measurement of the shift and a better estimation of the optimal RI.

    %   Xavier Trepat 05/02/2007
    %   Raimon Sunyer 09/26/2006
    %	Based on Iva Marija Tolic-Norrelykke 'crop_onecell.m';
    %%
    % Marina Uroz: 3D.
    % Ernest Latorre and Manuel Gomez 04/21/2016: Reduce the RAM ammount used for large z-stacks and other minor modifications.
    %                                             Eliminate the need to pre-process the z-stack in ImageJ to center and crop it 
    %                                             around the more focused plane. It is done automatically here.

    disp(' ');
    disp('         *** Start image cropping ***       ');
    disp(' ');

    % Paddings used to dedrift the images.
    padding = Settings.padding;                     % Padding for the edges of the the files.
    paddingz = Settings.paddingz;                   % Padding for the z direction.

    % Remove the warnings for the poorly-formed Tiff tags. If not, many warnings will appear with the use of TIFFStack.
    w1 = warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
    w2 = warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning');

    % Pre-allocate enough memory for the largest arrays.
    [Settings.SizeY, Settings.SizeX, ~] = size(TIFFStack(File.Name(1).Trypsin));        % Size of the Trypsin stack.
    num_images3 = 0;                                                        % Thickness of the 3D reconstruction stack, if any.

    % Define the ROI settings depending on the image size.
    % Settings.Size_ROI{x,y} indicate the size of the final cropped images.
    % Check that the cropped images will not span more than the original images.
    Settings.Size_ROIx = min(Settings.Size_ROIx, Settings.SizeX);
    Settings.Size_ROIy = min(Settings.Size_ROIy, Settings.SizeY);

    % Settings.CenterRoi{X,Y} indicate the center of the final cropped images.
    if ~isfield(Settings, 'CenterRoiX')
        Settings.CenterRoiX = round(Settings.SizeX/2);
    end

    if ~isfield( Settings, 'CenterRoiY' )
        Settings.CenterRoiY = round(Settings.SizeY/2);
    end

    Settings.CenterRoiX = max(1, min( Settings.CenterRoiX, Settings.SizeX));
    Settings.CenterRoiY = max(1, min( Settings.CenterRoiY, Settings.SizeY));

    % Settings.Size_ROI{x,y}_Cropper indicate the size of the subimages used by the dedrifting algorithm.
    % Check that the de-drifted sub-images will not span more than the original image.
    Settings.Size_ROIx_Cropper = 2*min(Settings.CenterRoiX, round(Settings.Size_ROIx_Cropper/2));
    Settings.Size_ROIy_Cropper = 2*min(Settings.CenterRoiY, round(Settings.Size_ROIy_Cropper/2));

    % Settings.Roi.{i,j}{1,2} indicate the lower and upper bounds of the i and j indices of the sub-images used by the 
    % de-drifting algorithm.
    % Check that the sub-images will not span more than the original image.
    Settings.ROI.i1 = max(1, Settings.CenterRoiY-round(Settings.Size_ROIy_Cropper/2) + 1);
    Settings.ROI.i2 = min(Settings.SizeY, Settings.CenterRoiY+round(Settings.Size_ROIy_Cropper/2));
    Settings.ROI.j1 = max(1, Settings.CenterRoiX-round(Settings.Size_ROIx_Cropper/2) + 1);
    Settings.ROI.j2 = min(Settings.SizeX, Settings.CenterRoiX+round(Settings.Size_ROIx_Cropper/2));

    % Use f{i,j}{1,2} as short names for Settings.Roi.{i,j}{1,2}.
    fi1 = Settings.ROI.i1;
    fi2 = Settings.ROI.i2;
    fj1 = Settings.ROI.j1;
    fj2 = Settings.ROI.j2;

    % We calculate one timepoint each "stepbeads" times.
    stepbeads = Settings.StepBeads;

    % In order to calculate the displacements in the z-direction, we need a z-stack of "n" images around the beads surface. 
    % Typically, during the experiment, we acquire a larger stack. Later, we pre-process this stack in ImageJ by locating the 
    % central plane, choosing "n" planes around it, and deleting the remaining planes above and below them. This is a tedius 
    % and manual task. In this version of the code, we don't need to pre-process the data with ImageJ. We load the whole stack, 
    % and we find the central plane accornding to a suitable criterion. We then keep "n" planes centered around it.

    %  Load the TIFFStack object of the Trypsin image. The Tripsin image is common for the whole experiment. Load it only once.  
    tsStack = TIFFStack(File.Name(1).Trypsin);
    im1_corr = double(tsStack(:, :, :));

    [~, ~, num_images1] = size(im1_corr);       % Size of the Trypsin stack.

    % Calculate the center of the stack. The casuistic could be very diverse, depending on the geometry of the problem.
    if strcmp(Settings.exp_type, 'Wells')
        % For the 'Wells' geometry, the center of the stack will be the middle plane between the top of the gel and the bottom 
        % of the well. The top of the gel can be found by calculating the mean intensity of each plane, and the peak of this 
        % profile will be the top plane. However, there is not a clear peak for the bottom of the well using this approach. By 
        % calculating the standard deviation of the intensity of each plane, we will have two peaks loosely located near the 
        % two planes of interest. This approach is not perfect, but it works just fine.

        % Calculate the std of the intensity of each plane.
        std_inte = squeeze(std(reshape(im1_corr, Settings.SizeX*Settings.SizeY, num_images1), 0, 1, 'omitnan'));

        [pks, locs] = findpeaks(std_inte);
        [~, ind] = maxk(pks, 2);

        % Location of the central plane, i.e. the middle between the top plane of the gel and the bottom of the well.
        cent_pos = round(sum(locs(ind))/2);

    else
        % For the default case, the center of the stack is the most in-focush plane.

        % Calculate the average intensity of each plane.
        ave_inte = squeeze(mean(reshape(im1_corr, Settings.SizeX*Settings.SizeY, num_images1), 1, 'omitnan'));

        [~, cent_pos] = max(ave_inte);      % Location of the plane with maximum intensity, i.e. most in-focus plane.

    end

    if Settings.Do.Open.Fluo_1            % 3D reconstruction image.
        num_images3 = size(TIFFStack(File.Name(1).Fluo_1), 3);              % Thickness of the 3D reconstruction stack.
    elseif Settings.Do.Open.Fluo_2        % 3D reconstruction image.
        num_images3 = size(TIFFStack(File.Name(1).Fluo_2), 3);              % Thickness of the 3D reconstruction stack.
    elseif Settings.Do.Open.Fluo          % 3D reconstruction image.
        num_images3 = size(TIFFStack(File.Name(1).Fluo), 3);                % Thickness of the 3D reconstruction stack.
    end

    if Settings.Round_Fourier_Size
        P = nextpow2(Settings.num_planes);              % Change the size of the stacks to the nearest power of two.
        num_planes = 2^P;
    else
        num_planes = Settings.num_planes;
    end

    Settings.SizeZ = num_planes;

    im1 = zeros(Settings.SizeY, Settings.SizeX, num_planes, 'uint16');      % Array to store the Tripsin image.

    % We keep only "num_planes" centered around the most central plane. Here we check the boundaries to make sure we have 
    % enough planes.
    % If we don't have enough planes above or below, we try to keep "num_planes" by moving the in-focus plane from the center.
    % Anything else is zero-padded.
    if ((num_images1 - cent_pos) > (num_planes - round(num_planes/2))) && (round(num_planes/2) > cent_pos)
        % Case: After centering both stacks, the original stack overflows (or just fits) in the upper half of the new stack, 
        % and falls short in the lower part of the new stack.
        ind0 = 1;
        indf = min(num_images1, num_planes);

        ind0_d = (ind0 - cent_pos) + round(num_planes/2) - ...
                 (min(num_images1 - (num_planes - round(num_planes/2)), round(num_planes/2)) - cent_pos);
        indf_d = (indf - cent_pos) + round(num_planes/2) - ...
                 (min(num_images1 - (num_planes - round(num_planes/2)), round(num_planes/2)) - cent_pos);

        [~, temp, ext] = fileparts(File.Name(1).Trypsin);
        disp(['Not enough planes above the maximum intensity in ', [temp, ext], '.']);

    elseif ((num_images1 - cent_pos) < (num_planes - round(num_planes/2))) && (round(num_planes/2) < cent_pos)
        % Case: After centering both stacks, the original stack falls short in the upper half of the new stack, 
        % and overflows (or just fits) in the lower part of the new stack.
        ind0 = max(1, num_images1-num_planes+1);
        indf = num_images1;

        ind0_d = (ind0 - cent_pos) + round(num_planes/2) + ...
                 (min((num_planes - round(num_planes/2)) - num_images1, -round(num_planes/2)) + cent_pos);
        indf_d = (indf - cent_pos) + round(num_planes/2) + ...
                 (min((num_planes - round(num_planes/2)) - num_images1, -round(num_planes/2)) + cent_pos);

        [~, temp, ext] = fileparts(File.Name(1).Trypsin);
        disp(['Not enough planes below the maximum intensity in ', [temp, ext], '.']);

    else
        ind0 = max(1, cent_pos-round(num_planes/2)+1);
        indf = min(num_images1, ind0+num_planes-1);

        ind0_d = (ind0 - cent_pos) + round(num_planes/2);
        indf_d = (indf - cent_pos) + round(num_planes/2);

    end

    im1(:, :, ind0_d:indf_d) = im1_corr(:, :, ind0:indf);
    clear('im1_corr');

    % Define the ROI. This contains the region with the cells of interest.
    % {i,j}{1,2} indicate the lower and upper bounds of the i and j indices of the sub-images that constitute the final cropped 
    % images.
    % Check that the sub-images will not span more than the original image.
    i1 = max(1, Settings.CenterRoiY - Settings.Size_ROIy/2+1);
    i2 = min(Settings.SizeY, Settings.CenterRoiY + Settings.Size_ROIy/2);
    j1 = max(1, Settings.CenterRoiX - Settings.Size_ROIx/2+1);
    j2 = min(Settings.SizeX, Settings.CenterRoiX + Settings.Size_ROIx/2);

    % Here, we are storing the cropped tripsin image, but we are not using this cropped stack any more. Thus, 
    % we don't need to store it in a variable. If we don't store it and instead we crop it on-the-fly when 
    % saving it, we save a lot of RAM for large stacks.
    File.Name(1).CropTrypsin = [File.CroppPath, filesep, 'crop', '_', 'trypsin', '_', num2str(1), '.tif'];

    imwrite(im1(i1:i2, j1:j2, 1), File.Name(1).CropTrypsin, 'compression', 'none');

    for kk=2:num_planes
        imwrite(im1(i1:i2, j1:j2, kk), File.Name(1).CropTrypsin, 'WriteMode', 'append', 'compression', 'none');
    end

    % When the stacks are already cropped, we will skip them. However, there is some information in the "File.mat" and 
    % "Settings.mat", obtained when the images where actually cropped, which we need to access.
    if exist([File.pathname, filesep, 'File.mat'], 'file') && exist([File.pathname, filesep, 'Settings.mat'], 'file')

        File_old = load([File.pathname, filesep, 'File.mat']);

        if isfield(File_old, 'File')
            File_old = File_old.File;
        end

        save([File.pathname, filesep, 'File_old.mat'], 'File_old', '-mat');

        n_elem_old = numel(File_old.Name);

    end

    % If Settings.mat or File.mat don't exist, we will need to re-crop the stacks.
    re_calculate = ~exist([File.pathname, filesep, 'File.mat'], 'file') || ...
                   ~exist([File.pathname, filesep, 'Settings.mat'], 'file');

    % Major loop that goes through each timepoint.
    for ii=1:stepbeads:File.NFiles.Beads

        % Skip the cropping of the stacks only if: it has previously been cropped and "File.mat" and "Settings.mat" do exist.
        if ~exist([File.CroppPath, filesep, 'crop', '_', 'Beads', '_', num2str(ii), '.tif'], 'file') || re_calculate || ...
                (ii > n_elem_old)

            [~, temp, ext] = fileparts(File.Name(ii).Beads);
            disp(['Processing ', [temp, ext], ' of ', num2str(File.NFiles.Beads), ' images']);

            % Array to store the padded Beads image. Since it is larger than the original Beads image, we will store only the 
            % padded array. This way, we will save  a lot of RAM for large stacks.
            im2 = zeros(Settings.SizeY + 2*padding, Settings.SizeX + 2*padding, num_planes + 2*paddingz, 'uint16');

            %  Load the TIFFStack object.
            tsStack = TIFFStack(File.Name(ii).Beads);
            im2_corr = double(tsStack( :, :, : ));

            [~, ~, num_images2] = size(im2_corr);       % Size of the tilt-corrected Trypsin stack.

            % Center the stack around a the central plane.
            if strcmp(Settings.exp_type, 'Wells')
                % For the wells, the upper plane of the gel can be selected as the plane with highest mean intensity, but we 
                % don't have a way to detect the plane of the bottom of the well by looking at the mean intensity. By looking 
                % at the standard deviation of the intensity, the upper and bottom planes are roughly identified as the largest 
                % peaks. This criterion can be used to roughly locate the center of the well.

                % Calculate the std of the intensity of each plane.
                std_inte = squeeze(std(reshape(im2_corr, Settings.SizeX*Settings.SizeY, num_images2), 0, 1, 'omitnan'));

                [pks, locs] = findpeaks(std_inte);
                [~, ind] = maxk(pks, 2);

                cent_pos = round(sum(locs(ind))/2);         % Location of the central plane of the well.

            else
                % For other configurations, select the most in focus plane, i.e. the plane with highest mean intensity.
                % Calculate the average intensity of each plane.
                ave_inte = squeeze(mean(reshape(im2_corr, Settings.SizeY*Settings.SizeX, num_images2), 1, 'omitnan'));
                [~, cent_pos] = max(ave_inte);      % Location of the plane with maximum intensity, i.e. most in-focus plane.

            end

            % We keep only "num_planes + 2*paddingz" planes centered around the most in-focus plane. Here we check the 
            % boundaries to make sure we have enough planes.
            if ((num_images2 - cent_pos) > (num_planes - round(num_planes/2)) + paddingz) && ...
                    (round(num_planes/2)+paddingz > cent_pos)
                % Case: After centering both stacks, the original stack overflows (or just fits) in the upper half of the new 
                % stack, and falls short in the lower part of the new stack.
                ind0 = 1;
                indf = min(num_images2, num_planes+2*paddingz);

                ind0_d = (ind0 - cent_pos) + round(num_planes/2) + paddingz - ...
                         (min(num_images2 - (num_planes - round(num_planes/2)) - paddingz, ...
                              round(num_planes/2) + paddingz) - cent_pos);
                indf_d = (indf - cent_pos) + round(num_planes/2) + paddingz - ...
                         (min(num_images2 - (num_planes - round(num_planes/2)) - paddingz, ...
                              round(num_planes/2) + paddingz) - cent_pos);

                [~, temp, ext] = fileparts(File.Name(1).Trypsin);
                disp(['Not enough planes above the maximum intensity in ', [temp, ext], '.']);

            elseif ((num_images2 - cent_pos) < (num_planes - round(num_planes/2)) + paddingz) && ...
                    (round(num_planes/2)+paddingz < cent_pos)
                % Case: After centering both stacks, the original stack falls short in the upper half of the new stack, 
                % and overflows (or just fits) in the lower part of the new stack.
                ind0 = max(1, num_images2-num_planes-2*paddingz+1);
                indf = num_images2;

                ind0_d = (ind0 - cent_pos) + round(num_planes/2) + paddingz + ...
                         (min((num_planes - round(num_planes/2)) + paddingz - num_images2, ...
                              -round(num_planes/2) - paddingz) + cent_pos);
                indf_d = (indf - cent_pos) + round(num_planes/2) + paddingz + ...
                         (min((num_planes - round(num_planes/2)) + paddingz - num_images2, ...
                              -round(num_planes/2) - paddingz) + cent_pos);

                [~, temp, ext] = fileparts(File.Name(1).Trypsin);
                disp(['Not enough planes below the maximum intensity in ', [temp, ext], '.']);

            else
                ind0 = max(1, cent_pos-round(num_planes/2)+1-paddingz);
                indf = min(num_images2, ind0+num_planes-1+2*paddingz);

                ind0_d = (ind0 - cent_pos) + round(num_planes/2) + paddingz;
                indf_d = (indf - cent_pos) + round(num_planes/2) + paddingz;

            end

            % Load the z-stack.
            im2((1:Settings.SizeY)+padding, (1:Settings.SizeX)+padding, ind0_d:indf_d) = im2_corr(:, :, ind0:indf);

            % In this function, we input a cropped section of matrices im1 and im2. These sub-matrices will not be used again, 
            % so we don't need to store them, and we save a lot of RAM for large stacks.
            [~, ~, ~, shx_temp, shy_temp, shz_temp, MaxCorr] = disp_on_blocks_fast_3D(...
                double(im1(fi1:fi2, fj1:fj2, :)), ...           % reference image
                double(im2((fi1:fi2)+padding, (fj1:fj2)+padding, (1:num_planes)+paddingz)), ...         % measurement image
                fi2-fi1+1, ...                                  % resolution of the PIV grid
                fj2-fj1+1, ...                                  % resolution of the PIV grid
                0, ...                                          % overlap between blocks
                padding, ...                                    % padding for the PIV, usually 0
                paddingz, ...                                   % padding in the z-direction
                1, ...                                          % window for FFT calculation. 1 means hanning
            	1, ...                                          % method to compute max of cross-correlation. 1 means 3DCentroid
                0, ...                                          % threshold of center of mass calculation
                2, ...                                          % size of the window
                3 ...                                           % Maximum number of iteration.
                );

            shx = round(shx_temp);
            shy = round(shy_temp);
            shz = round(shz_temp);

            File.Drift(ii).shx_centroid = shx;
            File.Drift(ii).shy_centroid = shy;
            File.Drift(ii).shz_centroid = shz;
            File.Drift(ii).shx_centroid_subpix = shx_temp;
            File.Drift(ii).shy_centroid_subpix = shy_temp;
            File.Drift(ii).shz_centroid_subpix = shz_temp;
            File.Drift(ii).pkh = max(MaxCorr);

            disp(['     ... current image displacement is ', '(', num2str(File.Drift(ii).shx_centroid), ',', ...
                  num2str(File.Drift(ii).shy_centroid), ',', num2str(File.Drift(ii).shz_centroid), ') pixels... ', ...
                  ' and maximum correlation is ', num2str(File.Drift(ii).pkh)]);

            % Crop and save:
            % If correlation is smaller than this, or displacement is NaN, it means that something wrong happened to your image.
            if (File.Drift(ii).pkh > Settings.MinCorrelationTrheshold) && ...
                    ~isnan(shx_temp) && ~isnan(shy_temp) && ~isnan(shz_temp)

                File.Mask(ii).Crop = 1;             % This marks this file as good.

                % For very large z-displacements, larger than paddingz, it might be possible to load more planes and use them.
                if (floor(abs(shz_temp)) >= paddingz) && (ind0 > 1) && (indf < num_images2)

                    % We keep only "num_planes" centered around the central plane. Here we check the boundaries to make sure we 
                    % have enough planes.
                    if shz+cent_pos-round(num_planes/2)+1-paddingz < 1
                        id0 = 1;
                        idf = min(num_images2, num_planes+2*paddingz);

                    elseif shz+cent_pos+(num_planes-round(num_planes/2))+paddingz > num_images2
                        id0 = max(1, num_images2-num_planes+1-2*paddingz);
                        idf = num_images2;

                    else
                        id0 = shz+cent_pos-round(num_planes/2)+1-paddingz;
                        idf = id0+num_planes-1+2*paddingz;

                    end

                    im2(:) = 0.;

                    % The indices id0 and idf indicate the range of planes from im2_corr that we will load into im2. Ideally, 
                    % those indices will completely fill im2. ind0_delta and indf_delta indicate where the planes from im2_corr 
                    % will be placed in im2.
                    id0_delta = max(0, min(paddingz, id0-(shz+cent_pos-round(num_planes/2)+1-paddingz)));
                    idf_delta = max(0, min(paddingz, (id0+num_planes-1+2*paddingz)-idf));
                    
                    % Load the z-stack.
                    im2((1:Settings.SizeY)+padding, (1:Settings.SizeX)+padding, ...
                        ((1+id0_delta):(1 + idf-idf_delta + id0_delta-id0))) = im2_corr(:, :, id0:(idf - idf_delta));

                else
                    id0 = ind0;

                end

                clear('im2_corr');

                padval = 0 ;

                % Crop beads image.
                % Shift the matrix im2 and store it again in the same variable. This will save a considerable ammount of RAM 
                % for larger stacks.
                indi0_delta = 1 + max(0, shy);
                indj0_delta = 1 + max(0, shx);
                indk0_delta = 1 + max(0, shz-(id0-ind0));
                indif_delta = size(im2, 1) + min(0, shy);
                indjf_delta = size(im2, 2) + min(0, shx);
                indkf_delta = size(im2, 3) + min(0, shz-(id0-ind0));

                indi0_delta_2 = 1 + max(0, -shy);
                indj0_delta_2 = 1 + max(0, -shx);
                indk0_delta_2 = 1 + max(0, -(shz-(id0-ind0)));
                indif_delta_2 = size(im2, 1) + min(0, -shy);
                indjf_delta_2 = size(im2, 2) + min(0, -shx);
                indkf_delta_2 = size(im2, 3) + min(0, -(shz-(id0-ind0)));

                im2(indi0_delta_2:indif_delta_2, indj0_delta_2:indjf_delta_2, indk0_delta_2:indkf_delta_2) = ...
                    im2(indi0_delta:indif_delta, indj0_delta:indjf_delta, indk0_delta:indkf_delta);

                im2(1:(indi0_delta_2-1), :, :) = padval;
                im2((indif_delta_2+1):end, :, :) = padval;
                im2(:, 1:(indj0_delta_2-1), :) = padval;
                im2(:, (indjf_delta_2+1):end, :) = padval;
                im2(:, :, 1:(indk0_delta_2+1)) = padval;
                im2(:, :, (indkf_delta_2+1):end) = padval;

                % Here, we are storing the cropped beads image, but we are not using this cropped stack any more. Thus, 
                % we don't need to store it in a variable. If we don't store it and instead we crop it on-the-fly when 
                % saving it, we save a lot of RAM for large stacks.
                File.Name(ii).CropBeads = [File.CroppPath, filesep, 'crop', '_', 'Beads', '_', num2str(ii), '.tif'];
                imwrite(uint16(im2((i1:i2)+padding, (j1:j2)+padding, 1+paddingz)), File.Name(ii).CropBeads, ...
                        'compression', 'none');

                for kk=2:num_planes
                    imwrite(uint16(im2((i1:i2)+padding, (j1:j2)+padding, kk+paddingz)), File.Name(ii).CropBeads, ...
                            'WriteMode', 'append', 'compression', 'none');
                end

                % Crop phase contast image.
                if Settings.Do.Open.BF

                    padval = 0;

                    File.Name(ii).CropBF = [File.CroppPath, filesep, 'crop','_', 'BF', '_', num2str(ii), '.tif'];

                    %  Load the TIFFStack object.
                    tsStack = TIFFStack(File.Name(ii).BF);

                    % Dedrift.
                    imPhase = ones(size(im2, 1), size(im2, 2), size(tsStack, 3))*padval;

                    imPhase((1:Settings.SizeY)+padding, (1:Settings.SizeX)+padding, :) = tsStack(:, :, :);

                    indi0_delta = 1 + max(0, shy);
                    indj0_delta = 1 + max(0, shx);
                    indif_delta = size(imPhase, 1) + min(0, shy);
                    indjf_delta = size(imPhase, 2) + min(0, shx);

                    indi0_delta_2 = 1 + max(0, -shy);
                    indj0_delta_2 = 1 + max(0, -shx);
                    indif_delta_2 = size(imPhase, 1) + min(0, -shy);
                    indjf_delta_2 = size(imPhase, 2) + min(0, -shx);

                    imPhase(indi0_delta_2:indif_delta_2, indj0_delta_2:indjf_delta_2, :) = ...
                        imPhase(indi0_delta:indif_delta, indj0_delta:indjf_delta, :);

                    imPhase(1:(indi0_delta_2-1), :, :) = padval;
                    imPhase((indif_delta_2+1):end, :, :) = padval;
                    imPhase(:, 1:(indj0_delta_2-1), :) = padval;
                    imPhase(:, (indjf_delta_2+1):end, :) = padval;

                    % Stored the cropped BF image.
                    imwrite(uint16(imPhase((i1:i2)+padding, (j1:j2)+padding, 1)), File.Name(ii).CropBF, 'compression', 'none');

                    % Repeat for each additional plane.
                    for kk = 2:size(tsStack, 3)
                        imwrite(uint16(imPhase((i1:i2)+padding, (j1:j2)+padding, kk)), File.Name(ii).CropBF, ...
                                'WriteMode', 'append', 'compression', 'none');
                    end

                end

                % Crop the Fluo_1 image.
                if  Settings.Do.Open.Fluo_1

                    padval = 0;

                    % Load the TIFFStack object.
                    tsStack = TIFFStack(File.Name(ii).Fluo_1);

                    File.Name(ii).CropFluo_1 = [File.CroppPath, filesep, 'crop', '_', 'Fluo_1', '_', num2str(ii), '.tif'];

                    % Array to store the padded Fluo_1 image. Since it is larger than the original Fluo_1 image, we will store 
                    % only the padded array. This way, we will save  a lot of RAM for large stacks.
                    im3 = zeros(Settings.SizeY + 2*padding, Settings.SizeX + 2*padding, num_planes + 2*paddingz, 'uint16');

                    im3_corr = tsStack( :, :, : ) ;

                    % Load the z-stack.
                    % For very large z-displacements, larger than paddingz, it might be possible to load more planes and use 
                    % them.
                    if (floor(abs(shz_temp)) >= paddingz) && (ind0 > 1) && (indf < num_images2)

                        im3(:) = 0.;

                        % Load the z-stack.
                        im3((1:Settings.SizeY)+padding, (1:Settings.SizeX)+padding, ...
                            ((1+id0_delta):(1 + idf-idf_delta + id0_delta-id0))) = im3_corr(:, :, id0:(idf - idf_delta));

                    else
                        im3((1:Settings.SizeY)+padding, (1:Settings.SizeX)+padding, ind0_d:indf_d) = im3_corr(:, :, ind0:indf);

                    end

                    clear('im3_corr');

                    % Crop Fluo_1 image.
                    % Shift the matrix im3 and store it again in the same variable. This will save a considerable ammount of RAM 
                    % for larger stacks.
                    im3(indi0_delta_2:indif_delta_2, indj0_delta_2:indjf_delta_2, indk0_delta_2:indkf_delta_2) = ...
                        im3(indi0_delta:indif_delta, indj0_delta:indjf_delta, indk0_delta:indkf_delta);

                    im3(1:(indi0_delta_2-1), :, :) = padval;
                    im3((indif_delta_2+1):end, :, :) = padval;
                    im3(:, 1:(indj0_delta_2-1), :) = padval;
                    im3(:, (indjf_delta_2+1):end, :) = padval;
                    im3(:, :, 1:(indk0_delta_2+1)) = padval;
                    im3(:, :, (indkf_delta_2+1):end) = padval;

                    % Stored the cropped Fluo_1 image.
                    imwrite(uint16(im3((i1:i2)+padding, (j1:j2)+padding, 1+paddingz)), File.Name(ii).CropFluo_1, ...
                            'compression', 'none');

                    for kk=2:num_planes
                        imwrite(uint16(im3((i1:i2)+padding, (j1:j2)+padding, kk+paddingz)), File.Name(ii).CropFluo_1, ...
                                'WriteMode', 'append', 'compression', 'none');
                    end

                end

                % Crop Fluo_2 image.
                if  Settings.Do.Open.Fluo_2

                    padval = 0;

                    %  Load the TIFFStack object.
                    tsStack = TIFFStack(File.Name(ii).Fluo_2);

                    File.Name(ii).CropFluo_2 = [File.CroppPath, filesep, 'crop', '_', 'Fluo_2', '_', num2str(ii), '.tif'];

                    % Repeat for each plane of the z-stack.
                    im3 = zeros(Settings.SizeY + 2*padding, Settings.SizeX + 2*padding, num_images3 + 2*paddingz);
                    im3((1:Settings.SizeY)+padding, (1:Settings.SizeX)+padding, (1:num_images3)+paddingz) = tsStack(:, :, :);

                    indi0_delta = 1 + max(0, shy);
                    indj0_delta = 1 + max(0, shx);
                    indk0_delta = 1 + max(0, shz);
                    indif_delta = size(im3, 1) + min(0, shy);
                    indjf_delta = size(im3, 2) + min(0, shx);
                    indkf_delta = size(im3, 3) + min(0, shz);

                    indi0_delta_2 = 1 + max(0, -shy);
                    indj0_delta_2 = 1 + max(0, -shx);
                    indk0_delta_2 = 1 + max(0, -shz);
                    indif_delta_2 = size(im3, 1) + min(0, -shy);
                    indjf_delta_2 = size(im3, 2) + min(0, -shx);
                    indkf_delta_2 = size(im3, 3) + min(0, -shz);

                    im3(indi0_delta_2:indif_delta_2, indj0_delta_2:indjf_delta_2, indk0_delta_2:indkf_delta_2) = ...
                        im3(indi0_delta:indif_delta, indj0_delta:indjf_delta, indk0_delta:indkf_delta);

                    im3(1:(indi0_delta_2-1), :, :) = padval;
                    im3((indif_delta_2+1):end, :, :) = padval;
                    im3(:, 1:(indj0_delta_2-1), :) = padval;
                    im3(:, (indjf_delta_2+1):end, :) = padval;
                    im3(:, :, 1:(indk0_delta_2+1)) = padval;
                    im3(:, :, (indkf_delta_2+1):end) = padval;

                    % Stored the cropped Fluo_2 image.
                    imwrite(uint16(im3((i1:i2)+padding, (j1:j2)+padding, 1+paddingz)), File.Name(ii).CropFluo_2, ...
                            'compression', 'none');

                    for kk=2:num_images3
                        imwrite(uint16(im3((i1:i2)+padding, (j1:j2)+padding, kk+paddingz)), File.Name(ii).CropFluo_2, ...
                                'WriteMode', 'append', 'compression', 'none');
                    end

                end

                % Crop Fluo 3D reconstruction image.
                if Settings.Do.Open.Fluo

                    padval = 0;

                    %  Load the TIFFStack object.
                    tsStack = TIFFStack(File.Name(ii).Fluo);

                    File.Name(ii).CropFluo = [File.CroppPath, filesep, 'crop', '_', 'Fluo', '_', num2str(ii), '.tif'];

                    % Repeat for each plane of the z-stack.
                    im3 = zeros(Settings.SizeY + 2*padding, Settings.SizeX + 2*padding, num_images3 + 2*paddingz);
                    im3((1:Settings.SizeY)+padding, (1:Settings.SizeX)+padding, (1:num_images3)+paddingz) = tsStack(:, :, :);

                    indi0_delta = 1 + max(0, shy);
                    indj0_delta = 1 + max(0, shx);
                    indk0_delta = 1 + max(0, shz);
                    indif_delta = size(im3, 1) + min(0, shy);
                    indjf_delta = size(im3, 2) + min(0, shx);
                    indkf_delta = size(im3, 3) + min(0, shz);

                    indi0_delta_2 = 1 + max(0, -shy);
                    indj0_delta_2 = 1 + max(0, -shx);
                    indk0_delta_2 = 1 + max(0, -shz);
                    indif_delta_2 = size(im3, 1) + min(0, -shy);
                    indjf_delta_2 = size(im3, 2) + min(0, -shx);
                    indkf_delta_2 = size(im3, 3) + min(0, -shz);

                    im3(indi0_delta_2:indif_delta_2, indj0_delta_2:indjf_delta_2, indk0_delta_2:indkf_delta_2) = ...
                        im3(indi0_delta:indif_delta, indj0_delta:indjf_delta, indk0_delta:indkf_delta);

                    im3(1:(indi0_delta_2-1), :, :) = padval;
                    im3((indif_delta_2+1):end, :, :) = padval;
                    im3(:, 1:(indj0_delta_2-1), :) = padval;
                    im3(:, (indjf_delta_2+1):end, :) = padval;
                    im3(:, :, 1:(indk0_delta_2+1)) = padval;
                    im3(:, :, (indkf_delta_2+1):end) = padval;

                    % Stored the cropped Fluo image.
                    imwrite(uint16(im3((i1:i2)+padding, (j1:j2)+padding, 1+paddingz)), File.Name(ii).CropFluo, ...
                            'compression', 'none');

                    for kk=2:num_images3
                        imwrite(uint16(im3((i1:i2)+padding, (j1:j2)+padding, kk+paddingz)), File.Name(ii).CropFluo, ...
                                'WriteMode', 'append', 'compression', 'none');
                    end

                end

            else

                disp(' Correlation is too low. File discarded from further analysis! ');
                File.Mask(ii).Crop = 0;             % Mark it as a wrong file and do not save it.

            end

        end

    end

    save([File.pathname, filesep, 'Settings.mat'], 'Settings', '-mat');
    save([File.pathname, filesep, 'File.mat'], 'File', '-mat');

    disp(' ');
    disp('          *** End image cropping ***        ');
    disp(' ');

    % Restore the warnings for the poorly-formed Tiff tags.
    warning(w1);
    warning(w2);

end