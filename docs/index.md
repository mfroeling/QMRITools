## Welcome to QRMITools

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2530801.svg)](https://doi.org/10.5281/zenodo.2530801) [![status](http://joss.theoj.org/papers/ef8bfb6c31499845d353b6a5af0d6300/status.svg)](http://joss.theoj.org/papers/ef8bfb6c31499845d353b6a5af0d6300)

---

[![wolfram language](\QMRITools\images\wolfram language.png)](https://www.wolfram.com/language/) | [![wolfram workbench](\QMRITools\images\wolfram workbench.jpg)](https://www.wolfram.com/workbench/) | [![eclipse](\QMRITools\images\eclipse.png)](https://www.eclipse.org/) | [![Wolfram Mathematica](\QMRITools\images\wolfram mathematica.png)](http://www.wolfram.com/mathematica/)

---

`QMRITools` is written in Mathematica using Wolfram Workbench and Eclipse and contains a collection of tools and functions for processing quantitative MRI data. The toolbox does not provide a GUI and its primary goal is to allow for fast and batch data processing, and facilitate development and prototyping of new functions. The core of the toolbox contains various functions for data manipulation and restructuring.

The toolbox was developed mostly in the context of quantitative muscle, nerve and cardiac magnetic resonance imaging. The library of functions grows along with the research it is used for and started as a toolbox to analyze DWI data of muscle.

The toolbox is developed for the [Wolfram language](https://www.wolfram.com/language/) and maintained using [Wolfram workbench](https://www.wolfram.com/workbench/) for [eclipse](https://www.eclipse.org/) and runs in the latest version of [Wolfram Mathematica](http://www.wolfram.com/mathematica/).

***
 
## Content
* [Latest Release](#latest-release)
* [Documentation](\QMRITools\htmldoc\guide\QMRITools.html){:target="_blank"}
* [Demonstrations](#demonstrations)
* [Functionality](#functionality)

## Latest Release

The latesest release can be found [here](https://github.com/mfroeling/QMRITools/releases){:target="_blank"}. 

Install the toolbox in the Mathematica UserBaseDirectory > Applications.

`FileNameJoin[{$UserBaseDirectory, "Applications"}]`

Some functions of QMRITools call on external executables and software.
These executables need to be present in "QMRITools\Applications" and are included in the release.
If for any reason you want to use other (older/newer) versions you can replace them but functionality is not guaranteed.
For the latest version of these tools and their user license please visit their website.

* [dcm2niix](https://github.com/rordenlab/dcm2niix/)
	* dcm2niix.exe
* [Elastix](http://elastix.isi.uu.nl/)
	* elastix.exe
	* transformix.exe

All functionality is tested under Windows 10 with the latest Mathematica version.
The Mathematica code is cross platform compatible with the exception of the external tools which are compiled for each OS.
The toolbox provides compiled versions for each OS but their functionality is not guaranteed.
The Elastix version used is 4.9 with OpenCL support. Additionally Elastix needs to be compiles with the PCA metrics, all DTI related parameters and all affine related parameters.

Although cross platform compatibility is provided I have only limited options for testing so if any issues arise please let me know.  

## Demonstrations

The release contains a zip file [DemoAndTest.zip](https://github.com/mfroeling/QMRITools/releases/download/2.0/DemoAndTest.zip) in which there is a file ``demo.nb``, a folder ``DemoData`` and a folder ``Testing``. 
To have a global overview of the functionality of the toolbox you can download this folder and run the ``demo.nb``.
By default the ``demo.nb`` looks for the folders ``DemoData`` and ``Testing`` in the same folder as the notebook.

In the first section of the demo notebook the toolbox is loaded and two tests are performed. The first test is to check of all files that are needed to run the toolbox are present. The second test runs after the toolbox is loaded and checks if all the functions and their options that are defined are correct.

## Functionality

The toolbox contains over 250 Functions and options of processing and analyzing data.
A summary of the core functionality is listed below. 

![Overview](\QMRITools\images\overview.png)

* **Diffusion Analysis**
	* Signal drift correction 
	* LLS, WLLS and iWLLS methods
	* REKINDLE outlier detection
	* IVIM fitting (fixed parameters, back-projection and Bayesian fitting)
	* Parameter fitting using histogram analysis
	* Joining and sorting of multiple series of the same volume
	* Joining multiple stacks with slice overlap into one stack
* **Diffusion Gradients optimization**
	* Single and multi shell
	* Rotating and correcting Bmatrix
	* Actual b-value estimation by gradient sequence integration
	* Gradient visualization
* **Noise suppression**
	* LMMSE noise suppression
	* PCA noise suppression based on ramom matrix theory
	* Anisotropic tensor smoothing using diffusion filter.
* **Importing and Exporting**
	* Dicom data (classing and enhanced file format)
	* Nifti data (.nii and .img .hdr, supports .gz files)
	* Compatible with ExplorDTI and Viste for fiber tractography
* **Data visualization**
	* 2D 3D and 4D viewer
	* Multiple images: Transparent overlay, difference and, checkboard overlays
	* Legend bars and image labels
	* Saving to pdf, jpg, animated gif and movie
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
* **Dixon Reconstruction**
	* B0 phase unwrapping
	* DIXON iDEAL reconstruction with T2start
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