# Description:

This repository contains the code used in the analysis of the paper "3D micropatterned traction force microscopy: a technique to control three-dimensional cell shape while measuring cell-substrate force transmission", by Laura Faure et al.

It allows to calculate the 3D tractions exerted by a cell inside a micro-well micropatterned in a soft gel. It makes use of the following experimental data: 

  * A 3D microscopy stack image of fluorescent markers embedded on the surface of the gel in a deformed configuration (when the cell is present).
  * Another image after the cell has been removed, killed or relaxed (reference image).

This repository is organized in the following directories:

  * [3D_PIV](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/tree/main/3D_PIV) contains the code to calculate the 3D displacement field of the gel.
  * [Micropatterned_3D_Traction_Calculation](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/tree/main/Micropatterned_3D_Traction_Calculation) contains the code to calculate the exerted tractions from the previously calculated displacement field.
  * [Unfold_Well](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/tree/main/Unfold_Well) contains the code to project and separate the data in the different sections of the well: the upper surface of the gel, the well's walls and the well's bottom.
  * [Examples](https://github.com/xt-prc-lab/3D_Micropatterned_Traction_Force_Microscopy/tree/main/Examples) contains an example experiment to be analyzed.

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
