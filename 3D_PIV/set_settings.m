%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                    Function used to set the general settings of the 3D PIV code and traction calculation.                    %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       Settings [struct]: structure containing the different parameters for the analysis.                                     %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       Settings [struct]: structure containing the different parameters for the analysis.                                     %
%                                                                                                                              %
%   Last Revison Date: 29/01/2024                                                                                              %
%   Based on codes created by Xavier Trepat's group.                                                                           %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Settings] = set_settings(Settings)

    % Open files settings:
    Settings.Imgfmt = 'tif';                            % Image format of input images.
%     Settings.Imgfmt = 'ome.tif';
    Settings.Do.Open.Trypsin = 1;                       % Open the Trypsin images.
    Settings.Do.Open.BF = 0;                            % Open the Bright Field images.
    Settings.Do.Open.Beads = 1;                         % Open the images of the Fluorescent Beads.
    Settings.Do.Open.Fluo = 0;                          % Open the images of the fluorescence 3D reconstruction.
    Settings.Do.Open.Fluo_1 = 1;                        % Open other fluorescence channel.
    Settings.Do.Open.Fluo_2 = 0;                        % Open another fluorescence channel.
    Settings.StepBeads = 1;                             % Step between beads timepoint images. Select how many timepoints will 
                                                        % be anlyzed. 1 means every timepoint, n means 1 every n timepoints.

    % Cropping settings:
    Settings.Do.LimitNumberOfFrames = 1;                % Maximum number of the last frame to process. Will process:
                                                        % [1:Settings.StepBeads:Settings.Do.LimitNumberOfFrames].

    Settings.MinCorrelationTrheshold = 0.2;             % For drift correction purposes. If trypsin and fluorescence 
                                                        % image correlate less that this value, then don't use this data point.
    Settings.Size_ROIx = 700;                           % Size of the ROI to be analyzed. If larger than the image size, will 
    Settings.Size_ROIy = 700;                           % use the image size.
    Settings.Size_ROIx_Cropper = 700;                   % Size of the image used by the cropper.
    Settings.Size_ROIy_Cropper = 700;
    Settings.CenterRoiX = 350;                          % Center of the ROI analyzed. If not provided, will be the center of 
    Settings.CenterRoiY = 350;                          % the original image.

    Settings.exp_type = 'Wells';                        % Type of experiment.

    Settings.num_planes = 98;                           % Number of planes kept in the cropper.
    Settings.num_planes = 126;

    Settings.Round_Fourier_Size = 0;                    % If 0, use the resolution provided. 
                                                        % Otherwise, round to the nearest power of 2.

    % Displacements settings:
    Settings.Overlap = 0.75;                            % Overlap, in xy, of the boxes for PIV analysis. 0 means no overlap.
                                                        % 0.5 means overlap of half window, 0.75 means 75%, etc.
    Settings.Overlap_z = 0.75;                          % Overlap, in z, of the boxes for PIV analysis..
    Settings.Resolution = 32;                           % PIV box size in xy, in pixels. XY resolution. In the longdisp 
                                                        % algorithm, the initial size will be larger and reduced, after each 
                                                        % iteration, until the target size is reached.
    Settings.Resolution_z = 8;                          % PIV box size in z, in planes. Z resolution. Same as above.
    Settings.padding = 0;                               % Padding in the X and Y directions for the PIV boxes.
    Settings.paddingz = 4;                              % Padding in the Z direction for the PIV boxes.
    Settings.Blocksizedec = 16;                         % In the longdisp displacement finder, the first iterations will be 
                                                        % performed with a large window size. The window will shrink at every 
                                                        % iteration. This parameter is the amount by which each side of the 
                                                        % interrogation box shrinks, in xy, at every iteration.
    Settings.Blocksizedec_z = 2;                        % Same as above, in the z direction.

    Settings.fft_window = 1;                            % Window for FFT calculation. 1 means hanning. 
                                                        % Anything else means no windowing.
    Settings.corr_method = 2;                           % Method to compute the location of the peak of the cross-correlation. 
                                                        % 1 means the centroid, 2 means Gaussian fit.
    Settings.pk_window = 1;                             % Half-size of the window, around the centroid peak, used to locate the 
                                                        % peak of the cross-correlation.
    Settings.iter_max = 100;                            % Maximum number of iterations of the displacement finder, 
                                                        % after the initial iterations decreasing the box size.
    Settings.iter_epsilon = 0.01;                       % Maximum acceptable relative error of the iteration.
    Settings.N_conver = 5;                              % Stop after N_conver iterations without improving the convergence.
    Settings.longdisp_neigh = 10;                       % For each box where longdisp is required, apply longdisp to the 
    Settings.longdisp_neigh_z = 2;                      % neighbouring "longdisp_neigh" boxes around it.
    Settings.dynamic_window = 1;                        % When determining the initial window size, we could use different 
                                                        % criteria for the different boxes. 1 if the sizes of the PIV boxes are 
                                                        % determined dynamically; 2 if all the boxes larger than 
                                                        % Settings.Resolution have the maximum initial window size; 0 if all 
                                                        % the windows have the maximum window size.
    % Tractions settings: gel and lens settings.
    Settings = gel_settings(Settings);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
