# PathPlanning_GeomTools_AlignCoordinates

This repo supports regression-based methods that align one cartesian coordinate system to another.
<!-- PROJECT LOGO -->
<br>
<p align="center">
  <h2 align="center">
    PathPlanning_GeomTools_AlignCoordinates
  </h2>
  <pre align="center">
        <img src=".\Images\BigDipper_small.jpg" alt="main align coordinates picture" width="960" height="540">
        <!-- figcaption>Fig.1 - The typical progression of map generation.</figcaption -->
        <font size="-2">Photo by <a href="https://unsplash.com/@niclas_lundin?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Niclas Lundin</a> on <a href="https://unsplash.com/photos/LpLcrYtuMqs?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a>
    </font>
  </pre>
</p>

<p align="left">
    This AlignCoordinates code library hosts functions to perform regression fitting of one coordinate system to another. In the 2D case, we choose homogenous coordinates of the form below which includes translation via $(tx, ty)$, rotation by angle ($\theta$), and uniform scaling of the axes ($S$):

```Matlab
% translation via (tx, ty), rotation by angle (q), and uniform scaling of
% the axes (S):
%
% X' = [S*cos(q)  -S*sin(q)   tx ] * | x | 
% Y' = [S*sin(q)   S*cos(q)   ty ]   | y |
% 1  = [0            0        1  ]   | 1 |

```

</p>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About the Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="structure">Repo Structure</a>
     <ul>
        <li><a href="#directories">Top-Level Directories</li>
        <li><a href="#dependencies">Dependencies</li>
     </ul>
    <li><a href="#functions">Functions</li>
         <ul>
        <li><a href="#basic-support-functions">Basic Support Functions</li>
        <ul>
          <li><a href="#fcn_aligncoords_generate2dtransformmatrix">fcn_AlignCoords_generate2DTransformMatrix - fills in transform matrix given translation, rotation, and scaling.</li>
          <li><a href="#fcn_aligncoords_fillsamplepoints">fcn_AlignCoords_fillSamplePoints - fills in sample points in homogenous form, useful for testing transformations.</li>
          <li><a href="#fcn_aligncoords_regressionfitscalefactor">fcn_AlignCoords_regressionFitScaleFactor - performs regression fitting to find the best scale factor that matches one set of points to another.</li>
          <li><a href="#fcn_aligncoords_fitrotationkabsch">fcn_AlignCoords_fitRotationKabsch - performs regression fitting to find the best-fit rotation and translation that matches one set of points to another using the Kabsch algorithm.</li>
        </ul>
        <li><a href="#core-functions">Core Functions</li>
        <ul>
          <li><a href="#fcn_aligncoords_fit2dcoordinates">fcn_AlignCoords_fit2DCoordinates - performs regression fitting to find the transform that matches one 2D coordinate system to another. Uses the Kabsch and scaling algorithms.</li>
          <li><a href="#fcn_aligncoords_fitaffinexform">fcn_AlignCoords_fitAffineXform - performs regression fitting to find the affine transform that matches one point set another</li>
        </ul>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
     <ul>
     <li><a href="#general-usage">General Usage</li>
     <li><a href="#examples">Examples</li>
     </ul>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## About The Project

<!--[![Product Name Screen Shot][product-screenshot]](https://example.com)-->

Often there is a need to align one coordinate system with another, for example to align an image's coordinates with those of real-world coordinates. This library of codes shows how to regression fit points in one coordinate system that correspond to points in another coordinate system so that one can extract the scaling factor and/or rotation and/or translation, or even the generic affine transform from one coordinate system to another.

* Inputs:
  * a set of XY points in N x 2 format in one "base" coordinate system, that correspond row-by-row to points (in N x 2 format) in another "target" coordinate system.
  
* Outputs
  * The scale factor, rotation, translation, and/or transformation matrix that maps one set of points to another.

<a href="#pathplanning_geomtools_aligncoordinates">Back to top</a>

***

<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Installation

1. Make sure to run MATLAB 2020b or higher. Why? The "digitspattern" command used in the DebugTools utilities was released late 2020 and this is used heavily in the Debug routines. If debugging is shut off, then earlier MATLAB versions will likely work, and this has been tested back to 2018 releases.

2. Clone the repo

   ```sh
   git clone https://github.com/ivsg-psu/PathPlanning_GeomTools_AlignCoordinates
   ```

3. Run the main code in the root of the folder (script_demo_AlignCoordinates.m), this will download the required utilities for this code, unzip the zip files into a Utilities folder (.\Utilities), and update the MATLAB path to include the Utility locations. This install process will only occur the first time. Note: to force the install to occur again, delete the Utilities directory and clear all global variables in MATLAB (type: "clear global *").

4. Confirm it works! Run script_demo_AlignCoordinates. If the code works, the script should run without errors. This script produces numerous example images such as those in this README file.

<a href="#pathplanning_geomtools_aligncoordinates">Back to top</a>

***

<!-- STRUCTURE OF THE REPO -->
### Directories

The following are the top level directories within the repository:
<ul>
 <li>/Data folder: contains any example datasets used in the code.</li>
 <li>/Documents folder: contains reference documents used for the creation of the code.</li>
 <li>/Functions folder: The majority of the codes are found within functions in this directory. All functions as well as test scripts are provided.</li>
 <li>/Images folder: images used for this README.md file are found in this directory.</li>
 <li>/Releases folder: contains current and historic releases of the code in zip file archives.</li>
 <li>/Utilities folder: Dependencies that are utilized but not implemented in this repository are placed in the Utilities directory. These can be single files but are most often folders containing other cloned repositories.</li>
</ul>

<a href="#pathplanning_geomtools_aligncoordinates">Back to top</a>

***

### Dependencies

* [Errata_Tutorials_DebugTools](https://github.com/ivsg-psu/Errata_Tutorials_DebugTools) - The DebugTools repo is used for the initial automated folder setup, and for input checking and general debugging calls within subfunctions. The repo can be found at: <https://github.com/ivsg-psu/Errata_Tutorials_DebugTools>

The dependencies are automatically installed by running the root master script (see below).

<a href="#pathplanning_geomtools_aligncoordinates">Back to top</a>

***

<!-- FUNCTION DEFINITIONS -->
## Functions

### Basic Support Functions

#### fcn_AlignCoords_generate2DTransformMatrix

The function fcn_AlignCoords_generate2DTransformMatrix fills in transform matrix given translation, rotation, and scaling. The user must specify the order of operations via a string, with 'T' for translation, 'R' for roation, and 'S' for scaling, and the the string corresponds to the left-to-right construction. For example, the following code generates the transformation matrix, $T$, such that $T = T_{rotate} \cdot T_{translate} \cdot T_{scale}$, and the point transformation is $p_{moved} = T \cdot p_{base}$:

```MATLAB
fig_num = 1;

% Fill in a sample transform matrix
S = 2;
tx = 2;
ty = 7;
theta = -40*pi/180;

order_string = 'RTS';
T = fcn_AlignCoords_generate2DTransformMatrix( ...
    S, theta, tx, ty, order_string, fig_num);
title('Demonstration of fcn_AlignCoords_fillSamplePoints', 'Interpreter', 'none');
```

The results of this transform are shown automatically, if a figure number is given

<pre align="center">
  <img src=".\Images\fcn_AlignCoords_generate2DTransformMatrix.png" alt="fcn_AlignCoords_generate2DTransformMatrix picture" width="400" height="300">
  <figcaption>Fig.1 - The function fcn_AlignCoords_generate2DTransformMatrix fills in transform matrix given translation, rotation, and scaling.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_geomtools_aligncoordinates">Back to top</a>

***

#### fcn_AlignCoords_fillSamplePoints

The function fcn_AlignCoords_fillSamplePoints fills in sample points in homogenous form, useful for testing transformations. Optional inputs allows users to specify the test point group, and/or plot the test points.

<pre align="center">
  <img src=".\Images\fcn_AlignCoords_fillSamplePoints.png" alt="fcn_AlignCoords_fillSamplePoints picture" width="400" height="300">
  <figcaption>Fig.2 - The function fcn_AlignCoords_fillSamplePoints creates test data sets for exercising lap functions.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_geomtools_aligncoordinates">Back to top</a>

***

#### fcn_AlignCoords_regressionFitScaleFactor

The function fcn_aligncoords_plotZoneDefinition performs regression fitting to find the best scale factor that matches one set of points to another.

<pre align="center">
  <img src=".\Images\fcn_AlignCoords_regressionFitScaleFactor.png" alt="fcn_AlignCoords_regressionFitScaleFactor picture" width="400" height="300">
  <figcaption>Fig.3 - The function fcn_AlignCoords_regressionFitScaleFactor performs regression fitting to match the base point's scale factor.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_geomtools_aligncoordinates">Back to top</a>

***

#### fcn_AlignCoords_fitRotationKabsch

The function fcn_AlignCoords_fitRotationKabsch performs regression fitting to find the best-fit rotation and translation that matches one set of points to another using the Kabsch algorithm.

<pre align="center">
  <img src=".\Images\fcn_AlignCoords_fitRotationKabsch.png" alt="fcn_AlignCoords_fitRotationKabsch picture" width="400" height="300">
  <figcaption>Fig.4 - The function fcn_AlignCoords_fitRotationKabsch performs regression fitting to find the best-fit rotation and translation that matches one set of points to another.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_geomtools_aligncoordinates">Back to top</a>

***

### Core Functions

#### fcn_AlignCoords_fit2DCoordinates

The function fcn_AlignCoords_fit2DCoordinates performs regression fitting to find the transform that matches one 2D coordinate system to another. Uses the Kabsch and scaling algorithms.

<pre align="center">
  <img src=".\Images\fcn_AlignCoords_fit2DCoordinates.png" alt="fcn_AlignCoords_fit2DCoordinates picture" width="400" height="300">
  <figcaption>Fig.5 - The function fcn_AlignCoords_fit2DCoordinates performs regression fitting to find the transform that matches one 2D coordinate system to another.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_geomtools_aligncoordinates">Back to top</a>

***

#### fcn_AlignCoords_fitAffineXform

The function fcn_AlignCoords_fitAffineXform performs regression fitting to find the affine transform that matches one point set another. In general, it can have better accuracy in fitting than fitting 2D coordinates, but the resulting transform does not preserve rotation. Both algorithms give nearly identical performance.

<pre align="center">
  <img src=".\Images\fcn_AlignCoords_fitAffineXform.png" alt="fcn_AlignCoords_fitAffineXform picture" width="400" height="300">
  <figcaption>Fig.6 - The function fcn_AlignCoords_fitAffineXform performs regression fitting to find the affine transform that matches one point set another.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<a href="#pathplanning_geomtools_aligncoordinates">Back to top</a>

***

<!-- USAGE EXAMPLES -->
## Usage
<!-- Use this space to show useful examples of how a project can be used.
Additional screenshots, code examples and demos work well in this space. You may
also link to more resources. -->

### General Usage

Each of the functions has an associated test script, using the convention

```sh
script_test_fcn_fcnname
```

where fcnname is the function name as listed above.

As well, each of the functions includes a well-documented header that explains inputs and outputs. These are supported by MATLAB's help style so that one can type:

```sh
help fcn_fcnname
```

for any function to view function details.

<a href="#pathplanning_geomtools_aligncoordinates">Back to top</a>

***

### Examples

1. Run the main script to set up the workspace and demonstrate main outputs, including the figures included here:

   ```sh
   script_demo_AlignCoordinates
   ```

    This exercises the main function of this code: fcn_AlignCoords_fit2DCoordinates

2. After running the main script to define the included directories for utility functions, one can then navigate to the Functions directory and run any of the functions or scripts there as well. All functions for this library are found in the Functions sub-folder, and each has an associated test script. Run any of the various test scripts, such as:

   ```sh
   script_test_fcn_aligncoords_breakDataIntoLapIndices
   ```

<a href="#pathplanning_geomtools_aligncoordinates">Back to top</a>

***

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.

<a href="#pathplanning_geomtools_aligncoordinates">Back to top</a>

***

## Major release versions

This code is still in development (alpha testing)

<a href="#pathplanning_geomtools_aligncoordinates">Back to top</a>

***

<!-- CONTACT -->
## Contact

Sean Brennan - sbrennan@psu.edu

Project Link: [https://github.com/ivsg-psu/PathPlanning_GeomTools_AlignCoordinates](https://github.com/ivsg-psu/PathPlanning_GeomTools_AlignCoordinates)

<a href="#pathplanning_geomtools_aligncoordinates">Back to top</a>

***

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->