# Welcome to QRMITools

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2530801.svg)](https://doi.org/10.5281/zenodo.2530801) 
[![status](http://joss.theoj.org/papers/ef8bfb6c31499845d353b6a5af0d6300/status.svg)](http://joss.theoj.org/papers/ef8bfb6c31499845d353b6a5af0d6300)
[![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/dwyl/esta/issues)

------------------------------------------------------------------------

[![wolfram language](https://github.com/mfroeling/QMRITools/blob/master/docs/images/wolfram_language.png)](https://www.wolfram.com/language/)
[![wolfram workbench](https://github.com/mfroeling/QMRITools/blob/master/docs/images/wolfram_workbench.jpg)](https://www.wolfram.com/workbench/)
[![eclipse](https://github.com/mfroeling/QMRITools/blob/master/docs/images/eclipse.png)](https://www.eclipse.org/)
[![Wolfram Mathematica](https://github.com/mfroeling/QMRITools/blob/master/docs/images/wolfram_mathematica.png)](http://www.wolfram.com/mathematica/)

------------------------------------------------------------------------
 
## Content
- [Introduction](#introduction)
- [Installation](#installation)
- [Demonstrations](#demonstrations)
- [Documentation](#documentation)
- [Using the toolbox](#using-the-toolbox)
- [Functionality](#functionality)
- [Toolboxes](#toolboxes)
- [License](#license)

------------------------------------------------------------------------

## Introduction

`QMRITools` is written in Mathematica using Wolfram Workbench and Eclipse and contains a collection of tools and functions for processing quantitative MRI data. The toolbox does not provide a GUI and its primary goal is to allow for fast and batch data processing, and facilitate development and prototyping of new functions. The core of the toolbox contains various functions for data manipulation and restructuring.

The toolbox was developed mostly in the context of quantitative muscle
(Froeling et al. 2012), nerve and cardiac magnetic resonance imaging.
The library of functions grows along with the research it is used for
and started as a toolbox to analyze DWI data of muscle. Since then it
has grown to include many other features such as cardiac analysis
(tagging and T1 mapping), dixon reconstruction, EPG modeling and
fitting, j-coupling simulations and more. It currently contains over 350 custom 
functions (over 20.000 lines of code) complete with documentation and demonstrations. 

![Muscle processing](https://github.com/mfroeling/QMRITools/blob/master/docs/images/processing.PNG)

The toolbox is developed for the [Wolfram language](https://www.wolfram.com/language/) and maintained using [Wolfram workbench](https://www.wolfram.com/workbench/) for [eclipse](https://www.eclipse.org/) and runs in the latest version of [Wolfram Mathematica](http://www.wolfram.com/mathematica/).

When using the toolbox please cite one of the following references:

1. Froeling M: QMRTools: a Mathematica toolbox for quantitative MRI analysis. J Open Source Softw 2019; 4:1204. [link](https://joss.theoj.org/papers/ef8bfb6c31499845d353b6a5af0d6300)
2. Froeling M, et al.: Reproducibility of diffusion tensor imaging in human forearm muscles at 3.0 T in a clinical setting. Magn Reson Med 2010; 64:1182-1190. [link](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.22477)
3. Froeling M, et al.: Diffusion-tensor MRI reveals the complex muscle architecture of the human forearm. J Magn Reson Imaging 2012; 36:237-248. [link](https://onlinelibrary.wiley.com/doi/10.1002/jmri.23608)

------------------------------------------------------------------------

## Installation

The latest release can be found [here](https://github.com/mfroeling/QMRITools/releases). 

Install the toolbox in the Mathematica UserBaseDirectory > Applications.

`FileNameJoin[{$UserBaseDirectory, "Applications"}]`

Some functions of QMRITools call on external executables and software.
These executables need to be present in "QMRITools\Applications" and are included in the release.
If for any reason you want to use other (older/newer) versions you can replace them but functionality is not guaranteed.
For the latest version of these tools and their user license please visit their website.

- [dcm2niix](https://github.com/rordenlab/dcm2niix/)
	- dcm2niix.exe
- [Elastix](http://elastix.isi.uu.nl/)
	- elastix.exe
	- transformix.exe

All functionality is tested under Windows 10 with the latest Mathematica version.
The Mathematica code is cross platform compatible with the exception of the external tools which are compiled for each OS.
The toolbox provides compiled versions for each OS but their functionality is not guaranteed.
The Elastix version used is 4.9 with OpenCL support. Additionally Elastix needs to be compiles with the PCA metrics, all DTI related parameters and all affine related parameters.

Although cross platform compatibility is provided I have only limited options for testing so if any issues arise please let me know.  

------------------------------------------------------------------------

## Demonstrations

The release contains a zip file [DemoAndTest.zip](https://github.com/mfroeling/QMRITools/releases/download/2.0/DemoAndTest.zip) in which there is a file ``demo.nb``, a folder ``DemoData`` and a folder ``Testing``. 
To have a global overview of the functionality of the toolbox you can download this folder and run the ``demo.nb``.
By default the ``demo.nb`` looks for the folders ``DemoData`` and ``Testing`` in the same folder as the notebook.

In the first section of the demo notebook the toolbox is loaded and two tests are performed. The first test is to check of all files that are needed to run the toolbox are present. The second test runs after the toolbox is loaded and checks if all the functions and their options that are defined are correct.

------------------------------------------------------------------------

## Documentation

Documentation of all functions and their options is fully integrated in the Mathematica documentation. The toolbox always works within the latest version of Mathematica and does not support any backward compatibility.
After the toolbox is installed correctly it should show up as a package in the Mathematica add-ons. 

![QMRITools package](https://github.com/mfroeling/QMRITools/blob/master/docs/images/addons.PNG)

All code and documentation is maintained and uploaded to github using [Workbench](https://www.wolfram.com/workbench/). An online version of the full documentation can be found [here](https://mfroeling.github.io/QMRITools/htmldoc/guide/QMRITools.html).

![Guides QMRITools](https://github.com/mfroeling/QMRITools/blob/master/docs/images/Guide.PNG)

------------------------------------------------------------------------

## Using the toolbox

The toolbox can be loaded by using: `` <<QMRITools` ``
If you want to monitor the package loading you can use: `` QMRITools`$Verbose = True; <<QMRITools` ``

A list of all QMRITools packages is generated by 
	
	QMRIToolsPackages[]

A list of all DTITools functions or functions per toolbox is generated by 

	QMRIToolsFunctions[]
	QMRIToolsFunctions["toolboxname"]
	
To print the documentation of all functions use

	QMRIToolsFuncPrint[]
	QMRIToolsFuncPrint["toolboxname"]

A list off all functions and their help can be found in ``All-Functions.nb``, which is alos availible as a [pdf file](https://github.com/mfroeling/QMRITools/releases/download/2.0/All-Functions.pdf).


QMRITools contains the following toolboxes:

-   CardiacTools
-   CoilTools
-   DenoiseTools
-   DixonTools
-   ElastixTools
-   GeneralTools
-   GradientTools
-   ImportTools
-   IVIMTools
-   JcouplingTools
-   MaskingTools
-   NiftiTools
-   PhysiologyTools
-   PlottingTools
-   ProcessingTools
-   ReconstructionTools
-   RelaxometryTools
-   SimulationTools
-   SpectroTools
-   VisteTools

Under development

-   TaggingTools

------------------------------------------------------------------------

## Functionality

The toolbox contains over 350 Functions and options of processing and analyzing data.
A summary of the core functionality is listed below. 

![Overview](https://github.com/mfroeling/QMRITools/blob/master/docs/images/overview.png)

* **Diffusion Analysis**
	* Signal drift correction 
	* LLS, WLLS and iWLLS methods
	* REKINDLE outlier detection
	* IVIM fitting (fixed parameters, back-projection and Bayesian fitting)
	* Parameter fitting using histogram analysis
	* Joining and sorting of multiple series of the same volume
	* Joining multiple stacks with slice overlap into one stack

![Joining of multiple stacks with overlap into one data-set](https://github.com/mfroeling/QMRITools/blob/master/docs/images/joining.png)

* **Diffusion Gradients optimization**
	* Single and multi shell
	* Rotating and correcting Bmatrix
	* Actual b-value estimation by gradient sequence integration
	* Gradient visualization
* **Noise suppression**
	* LMMSE noise suppression
	* PCA noise suppression based on random matrix theory.
	* Anisotropic tensor smoothing using diffusion filter.

![Noise suppression](https://github.com/mfroeling/QMRITools/blob/master/docs/images/registration.gif)	
	
* **Importing and Exporting**
	* Dicom data (classing and enhanced file format)
	* Nifti data (.nii and .img .hdr, supports .gz files)
	* Compatible with ExplorDTI and Viste for fiber tractography
* **Data visualization**
	* 2D 3D and 4D viewer
	* Multiple images: Transparent overlay, difference and, checkboard overlays
	* Legend bars and image labels
	* Saving to pdf, jpg, animated gif and movie

![PlotData](https://github.com/mfroeling/QMRITools/blob/master/docs/images/visualization.PNG)
	
* **Masking**
	* Automate and threshold masking
	* Extracting parameters form masks
	* Smoothing masks
	* Smoothing muscle segmentation
* **Motion and distortion correction (Registration using elastix)**
	* Rigid, affine, b-spline and cyclic registration 
	* nD to nD registration
	* Automated series processing 
	* Slice to slice motion correction of 3D and 4D data

![PloRegister DatatData](https://github.com/mfroeling/QMRITools/blob/master/docs/images/registration.png)
	
* **Dixon Reconstruction**
	* B0 phase unwrapping
	* DIXON iDEAL reconstruction with T2star
* **Relaxometry fitting**
	* T2 fitting
	* T1rho fitting
	* Tri Exponential T2 fitting
	* EPG based T2 fitting with slice profile
* **Simulation Framework**
	* Diffuison tensor simulation and analysis
	* Bloch and EPG simulations
	* Cardiac DTI models (fiber architecture)
* **Cardiac Diffusion analysis** 
	* Breathing motion correction
	* Corrupted slice rejection
	* Local myocardial coordinate system calculation
	* helix angle and fiber architecture matrix
	* AHA 17 parameter description
	* Transmural parameter description	

**Under Construction**

* **Spectra fitting**
	* Fitting spectra using a set of basis functions
	* Phase correction and line-width optimization
	* baseline correction
* **Reconstruction Tools**
	* Basic algorithms for complex coil combination.
	* CSI data reconstruction
	* Simple Image reconstruction
* **Tagging analysis**
	* Getting displacement fields from tagging data
	* Strain, trosion and rotation from cardiac data 

------------------------------------------------------------------------

## Toolboxes

### CardiacTools
A collection of tools to analyze cardiac data. The main features are cardiac shape analysis which allows defining the hard in a local myocardial coordinate system which allows quantifying and analyzing data. When the cardiac geometry is known there are functions to analyze qMRI parameters in the AH17 model [@Cerqueira2002] or perform transmural sampling of qMRI parameters. 
Most of the functionality is demonstrated in the `demo.nb`.

![Cardiac segmentation in the AHA-17 model and estimation of the local myocardial coordinate stystem.](https://github.com/mfroeling/QMRITools/blob/master/docs/images/cardiac.png)

### CoilTools
A collection of tools to evaluate complex multi-coil data. The functions are specific for analysis of multi-coil magnitude and noise data which allows quantifying per channel SNR. Furthermore, if complex coil sensitivity maps are available it allows performing SENSE g-factor maps simulations.   
This toolbox is not demonstrated in the `demo.nb`.

### DenoiseTools
The toobox provides two algorithms that allow denoising of DWI data. The first is based on and LMMSE framework [@Aja-Fernandez2008] and the second is based on a random matrix theory and Principal component analysis framework [@Veraart2015; @Veraart2016]. Furthermore, it provides an anisotropic filter for denoising the estimated diffusion tensor which provides more reliable fiber orientation analysis [@JeeEunLee].
Most of the functionality is demonstrated in the `demo.nb`.

### DixonTools
An IDEAL based Dixon reconstruction algorithm [@Reeder2005; @Yu2008]. The method provides multi-peak fitting B0 field and T2- correction. The toolbox also provides a function for unwrapping phase data in 2D and 3D based on a best path method [@Abdul-Rahman2007; @Herraez2002]. It also contains a function that allows simulating gradient echo Dixon data.
Most of the functionality is demonstrated in the `demo.nb`.

![IDEAL based Dixon reconstruction: fitted fat fractions as a function of the imposed fat fraction, SNR and B0 field offset.](https://github.com/mfroeling/QMRITools/blob/master/docs/images/dixon.png)

### ElastixTools
A wrapper that calls the Elastix registration framework [@Klein2010; @Shamonin2013]. The toolbox determines what registration or transformations need to be performed, exports the related data to a temp folder and calls an automatically generated command line script that performs the registration. After registration is completed the data is again loaded into Mathematica.
Most of the functionality is demonstrated in the `demo.nb`.

### GeneralTools
This toolbox provides core functions used in many other functions and features. The functions comprise amongst others: data cropping, mathematical and statistical operators that ignore zero values, and data rescaling, transformation and padding.
Most of the functionality is demonstrated in the `demo.nb`.

### GradientTools
The main feature is an algorithm that uses static repulsion [@Jones1999; @Froeling2016] to generate homogeneously distributed gradient directions for DWI experiments. It also provides functions to convert bval and bvec files to bmatrix and vice versa.
Most of the functionality is demonstrated in the `demo.nb`.

![The graphical user interface of the gradient generation tool.](https://github.com/mfroeling/QMRITools/blob/master/docs/images/gradients.png)

### ImportTools
Allows importing DCM data or DCM header attributes. These functions are rarely used since the toolbox mostly uses the NIfTY data format and provides tools to convert DCM to NIfTI via [dcm2niix](https://github.com/rordenlab/dcm2niix).
This toolbox is not demonstrated in the `demo.nb`.

### IVIMTools
The toolbox includes functions to perform IVIM fitting of DWI data. There are two main functions: non linear fitting and Bayesian fitting [@Orton2014]. 
Some of the functionality is demonstrated in the `demo.nb`.

![Visualization of IVIM fitting.](https://github.com/mfroeling/QMRITools/blob/master/docs/images/ivim.png)

### JcouplingTools
A toolbox that allows simulation of NMR spectra using Hamiltonians based on methods from [FID-A](https://github.com/CIC-methods/FID-A). It allows simulating large spin systems [@Castillo2011] and was mainly implemented to investigate fat spectra in TSE [@Stokes2013].
Most of the functionality is demonstrated in the `demo.nb`.

### MaskingTools
Tools for masking and homogenization of data. It provides functions for smoothing cutting and merging masks and functions for the evaluation of data within masks.
Most of the functionality is demonstrated in the `demo.nb`.

### NiftiTools
Import and export of the NIfTI file format. Part of the code is based on previously implemented [nii-converter](https://github.com/tomdelahaije/nifti-converter). For converting DICOM data to the NIfTI file format the toolbox uses [dcm2niix](https://github.com/rordenlab/dcm2niix/releases). It also provides some specialized NIfTI import functions for specific experiments which are probably not generalizable.
Most of the functionality is demonstrated in the `demo.nb`.

### PhysiologyTools
Functions for importing and analyzing Philips physiology logging and RespirAct trace files. The functions are rarely used and not well supported.
This toolbox is not demonstrated in the `demo.nb`.

### PlottingTools
A variety of functions for visualization of various data types. The main functions are 'PlotData' and 'PlotData3D' which allow viewing 2D, 3D and 4D data.
Most of the functionality is demonstrated in the `demo.nb`.

### ProcessingTools
The toolbox comprises a variety of functions that allow data manipulation and analysis. The main functions allow joining multiple data sets into one continuous data set [@Froeling2015] or to split data of two legs into two separate data-sets. Furthermore, it contains a collection of functions for data evaluation and analysis.
Most of the functionality is demonstrated in the `demo.nb`.

### RelaxometryTools
A collection of tools to fit T2, T2*, T1rho and T1 relaxometry data. The main function of this toolbox is an extended phase graph (EPG) [@Weigel2015] method for multi-compartment T2 fitting of multi-echo spin echo data [@Marty2016]. Therefore it provides functions to simulate and evaluate EPG. 
Some of the functionality is demonstrated in the `demo.nb`.

![Demonstration of EPG based T2 fitting: the fitted water T2 relaxation as a function of B1, SNR and fat fraction.](https://github.com/mfroeling/QMRITools/blob/master/docs/images/epg-t2.png)

### SimulationTools
The main purpose of this toolbox is to simulate DTI based DWI data and contains some functions to easily perform analysis of the fit results of the simulated signals [@Froeling2013].
Some of the functionality is demonstrated in the `demo.nb`.

### TensorTools
The original toolbox where the project started. The main functions in this toolbox are to fit and evaluate the diffusion tensor model. Various fitting methods are implemented (e.g. LLS, NLS, WLLS, and iWLLS). The default method is an iterative weighted linear least squares approach [@Veraart2013]. The tensor fitting also includes outlier detections using REKINDLE [@Tax2015] and data preparation includes drift correction [@Vos2014].
Most of the functionality is demonstrated in the `demo.nb`.

![MD and FA as a function of SNR and fat fraction. Results are from simulated data using an iWLLS algorithm with outlier rejection.](https://github.com/mfroeling/QMRITools/blob/master/docs/images/dti.png)

### VisteTools
Import and export functions for tensor data which can be used in the [vIST/e](http://bmia.bmt.tue.nl/software/viste/) tractography tool.
None of the functionality is demonstrated in the `demo.nb`.

------------------------------------------------------------------------

![Full leg diffusion tensor fiber tractography](https://github.com/mfroeling/QMRITools/blob/master/docs/images/animation-small.gif)

------------------------------------------------------------------------
	
## License
https://opensource.org/licenses/BSD-3-Clause

Note that restrictions imposed by these patents (and possibly others)
exist independently of and may be in conflict with the freedoms granted
in BSD-3-Clause license, which refers to copyright of the program, not patents
for any methods that it implements. Both copyright and patent law must
be obeyed to legally use and redistribute this program and it is not the
purpose of this license to induce you to infringe any patents or other
property right claims or to contest validity of any such claims.  If you
redistribute or use the program, then this license merely protects you
from committing copyright infringement.  It does not protect you from
committing patent infringement.  So, before you do anything with this
program, make sure that you have permission to do so not merely in terms
of copyright, but also in terms of patent law.

Some code in the NiiTools packages was based on https://github.com/tomdelahaije/nifti-converter

## References