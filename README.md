# Description:

This repository contains the code used in the analysis of the paper "3D micropatterned traction force microscopy: a technique to control three-dimensional cell shape while measuring cell-substrate force transmission", by Laura Faure et al.

It allows to calculate the 3D tractions exerted by a cell inside a micro-well micropatterned in a soft gel. It makes use of the following experimental data: 

  * A 3D microscopy stack image of fluorescent markers embedded on the surface of the gel in a deformed configuration (when the cell is present).
  * Another 3D microscopy stack image after the cell has been removed or relaxed (reference image).

This repository is organized in the following directories:

  * [3D_PIV](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/tree/main/3D_PIV) contains the code to calculate the 3D displacement field of the gel.
  * [Micropatterned_3D_Traction_Calculation](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/tree/main/Micropatterned_3D_Traction_Calculation) contains the code to calculate the exerted tractions from the previously calculated displacement field.
  * [Unfold_Well](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/tree/main/Unfold_Well) contains the code to project and separate the data in the different sections of the well: the upper surface of the gel, the well's walls and the well's bottom.
  * [Examples](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/tree/main/Examples) contains two example experiments to be analyzed.

# Prerequisites:

The software contained in this repository is written in Matlab and C programming languages. In order to run, it needs a valid installation of:

 * [Matlab](https://www.mathworks.com/products/matlab.html): tested with Matlab versions R2019a and R2020b.
 * A C compiler compatible with Matlab, such as [GCC](https://gcc.gnu.org/): tested with gcc versions 8.5.0 and 13.2.1.

Furthermore, in order to run the code in [Micropatterned_3D_Traction_Calculation](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/tree/main/Micropatterned_3D_Traction_Calculation), it needs a working installation of:
 * [Abaqus](https://www.3ds.com/es/productos-y-servicios/simulia/productos/abaqus/): tested with version 6.14.5.

It has only been tested under Linux (Gentoo Linux and Ubuntu 20.04), but it should work in any other operating system with minor modifications.

# Dependencies:

This code makes use of the following software dependencies, that must be downloaded and placed in the corresponding folder:

## 3D PIV:

  * [TIFFStack](https://github.com/DylanMuir/TIFFStack): fast loading of TIFF files into Matlab.

## Micropatterned 3D Traction Calculation:

  * [TIFFStack](https://github.com/DylanMuir/TIFFStack): fast loading of TIFF files into Matlab.
  * [fitcircle.m](https://www.mathworks.com/matlabcentral/fileexchange/15060-fitcircle-m): fit a set of points to a circle.
  * [inpaint_nans3](https://www.mathworks.com/matlabcentral/fileexchange/21214-inpainting-nan-elements-in-3-d): in-paints over nans in a 3-D array.
  * [removeOutliers.m](https://github.com/FranckLab/FIDVC/blob/master/Matlab%20Code/removeOutliers.m): remove outliers from PIV data.

## Unfold Well:

  * [TIFFStack](https://github.com/DylanMuir/TIFFStack): fast loading of TIFF files into Matlab.

# Instructions:

## Data Preparation:

### Microscopy image requirements:
The minimum dataset needed to run this code consists on:

  * A 3D microscopy stack image of fluorescent markers embedded on the surface of the gel in a deformed configuration (when the cell is present). It needs to span the whole depth of the well, and it must not cut the Point Spread Function of the florescent beads.
  * Another 3D microscopy stack image after the cell has been removed or relaxed (reference image).

Optionally, 

  * Any other 3D microscopy stack images of the sample can be provided and they will be aligned together with the previous two, but they will not be used for the deformation and traction calculations.

These files can be located either in the same or in different directories.

Here is an illustration of the example dataset provided [here](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/tree/main/Examples/Data/220412). It represents an orthogonal view of the data: XY slice in the upper left pannel, XZ slice in the lower left pannel and YZ slice in the upper right pannel. It is an overlay of 3 different channels: the red channel represents the beads in the reference configuration, the green channel are the beads in the deformed configuration and the blue channel represents a fluorescence channel of the cell.

<p align="center">
<img src="https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/blob/main/Images/Example_Orthogonal_Views.png" alt="Orthogonal_Views_of_Raw_Data" align="center" width=35%>
</p>

Two examples of datasets are provided [here](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/tree/main/Examples/Data).

### File image convention:

Any naming convention can be followed for the data. The only restriction is that the names of the different channels have to end with a number. This number represents the timepoint of the data aquisition. Files present in the same directory with the same name, only differing in the trailing number, will be considered as different timepoints of the same image series. Multiple timepoints will be automatically detected and can be succesively analyzed by this software. Optionally, another number series could be included before the time identifier, separated by another non-numeric character, to represent a position in the sample. This number can be used to sequentially analize multiple positions of the same sample.

As an example, the image stack provided [here](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/tree/main/Examples/Data/220412) has the following name:
```
220412_15kPa_19well_gel2_beads_f4_t0.tif
```
representing the timepoint 0 of position 4.

### Abaqus model:

In order to calculate the tractions, an Abaqus model of the well needs to be provided. Two types of Boundary Conditions need to be imposed: 

  * "encastre" at the lowest plane of the model, and
  * "user-defined" displacements at the top surface (including the well walls) of the model.

The model needs to be meshed with tetrahedral elements. An "\*.inp" Abaqus job file needs to be created and saved in a directory. Only the "\*.inp" file is needed for the traction calculation, but it is a good practice to also save the  "\*.cae" model in the same directory. Any file stored in the same directory than the  "\*.inp" job file will be copied over to the results folder. Two examples of Abaqus models can be found [here](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/tree/main/Examples/Abaqus_Models).

## 3D PIV (Displacement Calculation):

In order to calculate the 3D deformation, the [PIV_3D.m](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/blob/main/3D_PIV/PIV_3D.m) file needs to be run.

### Parameters:

Two different files contain the parameters for the PIV and the experiment, and need to be modified beforehand:

  * [set_settings.m](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/blob/main/3D_PIV/set_settings.m): contains the parameters for the PIV, such as the channels to analyze, the extension of the images, PIV box size and overlap, etc.
  * [gel_settings.m](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/blob/main/3D_PIV/gel_settings.m): contains parameters such as pixel size, filtering parameters, etc.

All of the parameters contained in those files are extensively commented.

By default, when [PIV_3D.m](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/blob/main/3D_PIV/PIV_3D.m) is run, a file selection dialog will be displayed for each file type to input (beads, trypsin and optionally other fluorescence channels). Any timepoint from the time series can be chosen, and the rest of the timepoints will be automatically detected. In order to automatize the calculation of multiple timepoints, this default behaviour can be overriden. Around line 116, multiple variable are created with the names of the files to analyze, such as File.BeadsName and File.TrypsinName. They make use of the variables BeadsName, TrypsinName, etc. defined around line 53. If the names contained in File.BeadsName, File.TrypsinName, etc. exist, the code will automatically analyze those and will not prompt to click and select the data files. This code is prepared to sequentially analyze multiple positons and timepoints without user input.

### Results:

A new directory will be created with the name indicated by File.pathname, defined around lines 53 and 94, and all the results will be placed there. In particular, for the example provided, the directory 220412_Results will be created. A directory will be created inside for each position analyzed. In particular, for the example provided, this will be the directory 220412_Results/f4. The following files and directories will be placed inside:

  * Settings.mat: structure containing the settings used in the analysis.
  * File.mat: structure containing parameters regarding the files analyzed.
  * Croppeddata: directory containing the stacks of the different channels, aligned and cropped to the desired size.
  * Displacements: directory containing the files storing the measured displacement field.

## Micropatterned 3D Traction Calculation:

In order to calculate the 3D tractions, the [Calculate_Tractions_abaqus_full3D.m](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/blob/main/Micropatterned_3D_Traction_Calculation/Calculate_Tractions_abaqus_full3D.m) file needs to be run. It will use the data stored in the "Results" folder created by the 3D PIV.

### Parameters:

Most parameters will be carried over from the PIV. However, any new parameter needed for the traction calculation will be placed in the [Calculate_Tractions_abaqus_full3D.m](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/blob/main/Micropatterned_3D_Traction_Calculation/Calculate_Tractions_abaqus_full3D.m) file. These parameters include the height of the model gel, the choice of "Constrained" or "Unconstrained" method, the application of the "Refraction Index Correction", the surface on which the displacements are imposed, the location of the Abaqus executable and model, the name of the Abaqus job file, etc. Those parameters are extensively commented in the file.

The parameters previously imposed in the PIV can be replaced here, writting them in a structure called Settings_new. An example is given in this file. Let's say the pixel size written for the PIV was not correct. It can be safely replaced, because the PIV will output its results in pixels, not in microns. Around line 63, new variables Settings_new.PixelSizeXY and Settings_new.PixelSizeZ are created with the correct values, and the code will store them in the Settings structure.

### Refraction Index Missmatch correction:

In order to apply the Refraction Index missmatch correction, the variable Settings_new.Correct_Disp_RI needs to be set to something other than zero. The empirical shape and values of the aberrations introduced by the Refraction Index missmatch needs to be coded in [correct_disp_RI_matched.m](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/blob/main/Micropatterned_3D_Traction_Calculation/correct_disp_RI_matched.m), and it will depend on the samples and microscopy conditions. For the examples provided here, two shape corrections are provided: one for the large wells, and the other for the small wells. In order to choose the appropriate correction, the variable Settings_new.Well_Size needs to be set to 15 for the small wells or 19 for the large wells.

### Constrained vs unconstrained method:

In order to calculate the tractions with the constrained method, the variable Settings_new.constrained needs to be set to 'Constrained'. Any other value, or if the variable doesn't exist will default to the unconstrained method.

### Results:
