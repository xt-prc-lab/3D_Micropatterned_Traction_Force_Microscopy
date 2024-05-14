%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                                   Function used to set the gel settings of the experiment.                                   %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       Settings [struct]: structure containing the different parameters for the analysis.                                     %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       Settings [struct]: structure containing the different parameters for the analysis.                                     %
%                                                                                                                              %
%   Last Revison Date: 29/01/2024                                                                                              %
%   Created by Manuel Gomez Gonzalez and Ernest Latorre Ibars                                                                  %
%                                                                                                                              %
%   References:                                                                                                                %
%       "Refractive index and axial distance measurements in 3-D microscopy";                                                  %
%       Visser, T. D., Oud, J. L., and Brakenhoff, G. J.; SPIE MILESTONE SERIES MS; 1996; Issue 131; pp 286–289.               %
%                                                                                                                              %
%       "Tutorial: avoiding and correcting sample-induced spherical aberration artifacts in 3D fluorescence microscopy";       %
%       Diel, Erin E., Jeff W. Lichtman, and Douglas S. Richardson; Nature protocols; 2020; Volume 15; Number 9, pp 2773-2784. %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Settings] = gel_settings(Settings)

    % Filter Settings.
    Settings.filt_size = 2;                             % Size of the Predictor filter, in PIV resolution units,
    Settings.filt_size_z = 2;                           % i.e. filter_size = filt_size/(1 - overlap).
    Settings.filt_size_trac = 2;                        % Size of the Corrector filter used for the traction calculation,
                                                        % in PIV resolution units. Same as above.
    Settings.filt_pow = 0.1;                            % Power of the Predictor filter.
    Settings.filt_pow_trac = 0.25;                      % Power of the Corrector filter.

    % Gel settings. At the moment, they are not used, they are introduced in tha Abaqus model.
%     Settings.GelHeight = 100;                           % Gel thickness, in microns.
%     Settings.Young = 15000;                             % Young modulus, in Pascals.
%     Settings.Poisson = 0.48;                            % Poisson ratio.

    % Lens settings:
    Settings.PixelSizeXY = 0.114;                       % Voxel size in the x, y and z directions, in microns.
    Settings.PixelSizeZ = 0.2;                          % As measured by the microscope's objective displacement. Will be 
                                                        % corrected, below, for Refractive Index mismatch.

    Settings.NA = 1.2;                                  % Numerical aperture of the lens.
    Settings.n1 = 1.456;                                % Refractive index of the medium of the lens,
                                                        % i.e. 1 for air, 1.33 for water, 1.52 for oil (and glass), 
                                                        % 1.456 for glycerin.
    Settings.n2 = 1.35;                                 % Refractive index of the medium where the beads are located,
                                                        % i.e. ~1.35 for polyacrylamide hydrogels, ~1.36 for cells,
                                                        % ~1.34 for cell culture media, ~1.40 for soft PDMS, etc.

    % Correct the measurements in Z due to the changes in Refraction Index in the path between the objective and the imaging 
    % plane (see the references above). We can either use an empirical correction or a theoretical mode.
    Settings.Corr_idx = 0.85;                           % Empirical correction index used, instead of the theoretical one.
%     Settings.Corr_idx = sqrt((Settings.n2^2 - Settings.NA^2)/(Settings.n1^2 - Settings.NA^2));      % Theoretical correction.

    % Correct for the refraction index missmatch between objective medium and sample.
%     Settings.GelHeight = Settings.GelHeight*Settings.Corr_idx;
    Settings.PixelSizeZ = Settings.PixelSizeZ*Settings.Corr_idx;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
