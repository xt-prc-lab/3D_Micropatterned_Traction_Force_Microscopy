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

The model needs to be meshed with tetrahedral elements. An "*.inp" Abaqus job file needs to be created and saved in a directory. Only the "*.inp" file is needed for the traction calculation, but it is a good practice to also save the  "*.cae" model in the same directory. Any file stored in the same directory than the  "*.inp" job file will be copied over to the results folder. Two examples of Abaqus models can be found [here](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/tree/main/Examples/Abaqus_Models).
