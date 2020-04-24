# Welcome to QRMITools

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2530801.svg)](https://doi.org/10.5281/zenodo.2530801)
[![status](http://joss.theoj.org/papers/ef8bfb6c31499845d353b6a5af0d6300/status.svg)](http://joss.theoj.org/papers/ef8bfb6c31499845d353b6a5af0d6300)
[![contributions
welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/dwyl/esta/issues)

------------------------------------------------------------------------

[![wolfram
language](https://github.com/mfroeling/QMRITools/blob/master/docs/images/wolfram_language.png)](https://www.wolfram.com/language/)
[![wolfram
workbench](https://github.com/mfroeling/QMRITools/blob/master/docs/images/wolfram_workbench.jpg)](https://www.wolfram.com/workbench/)
[![eclipse](https://github.com/mfroeling/QMRITools/blob/master/docs/images/eclipse.png)](https://www.eclipse.org/)
[![Wolfram
Mathematica](https://github.com/mfroeling/QMRITools/blob/master/docs/images/wolfram_mathematica.png)](http://www.wolfram.com/mathematica/)

------------------------------------------------------------------------

## Content

-   [Introduction](#introduction)
-   [Installation](#installation)
-   [Demonstrations](#demonstrations)
-   [Documentation](#documentation)
-   [Using the toolbox](#using-the-toolbox)
-   [Functionality](#functionality)
-   [Toolboxes](#toolboxes)
-   [License](#license)

------------------------------------------------------------------------

## Introduction

`QMRITools` is written in Mathematica using Wolfram Workbench and
Eclipse and contains a collection of tools and functions for processing
quantitative MRI data. The toolbox does not provide a GUI and its
primary goal is to allow for fast and batch data processing, and
facilitate development and prototyping of new functions. The core of the
toolbox contains various functions for data manipulation and
restructuring.

The toolbox was developed mostly in the context of quantitative muscle
(Froeling et al. 2012), nerve and cardiac magnetic resonance imaging.
The library of functions grows along with the research it is used for
and started as a toolbox to analyze DWI data of muscle. Since then it
has grown to include many other features such as cardiac analysis
(tagging and T1 mapping), dixon reconstruction, EPG modeling and
fitting, j-coupling simulations and more. It currently contains over 350 custom 
functions (over 20.000 lines of code) complete with documentation and demonstrations. 

The toolbox is developed for the [Wolfram
language](https://www.wolfram.com/language/) and maintained using
[Wolfram workbench](https://www.wolfram.com/workbench/) for
[eclipse](https://www.eclipse.org/) and runs in the latest version of
[Wolfram Mathematica](http://www.wolfram.com/mathematica/).

When using the toolbox please cite one of the following references:

1.  Froeling M: QMRTools: a Mathematica toolbox for quantitative MRI
    analysis. J Open Source Softw 2019; 4:1204.
    [link](https://joss.theoj.org/papers/ef8bfb6c31499845d353b6a5af0d6300)
2.  Froeling M, et al.: Reproducibility of diffusion tensor imaging in
    human forearm muscles at 3.0 T in a clinical setting. Magn Reson Med
    2010; 64:1182-1190.
    [link](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.22477)
3.  Froeling M, et al.: Diffusion-tensor MRI reveals the complex muscle
    architecture of the human forearm. J Magn Reson Imaging 2012;
    36:237-248.
    [link](https://onlinelibrary.wiley.com/doi/10.1002/jmri.23608)

------------------------------------------------------------------------

## Installation

The latest release can be found
[here](https://github.com/mfroeling/QMRITools/releases).

Install the toolbox in the Mathematica UserBaseDirectory &gt;
Applications.

`FileNameJoin[{$UserBaseDirectory, "Applications"}]`

Some functions of QMRITools call on external executables and software.
These executables need to be present in “QMRITools” and are included in
the release. If for any reason you want to use other (older/newer)
versions you can replace them but functionality is not guaranteed. For
the latest version of these tools and their user license please visit
their website.

-   [dcm2niix](https://github.com/rordenlab/dcm2niix/)
    -   dcm2niix.exe
-   [Elastix](http://elastix.isi.uu.nl/)
    -   elastix.exe
    -   transformix.exe

All functionality is tested under Windows 10 with the latest Mathematica
version. The Mathematica code is cross platform compatible with the
exception of the external tools which are compiled for each OS. The
toolbox provides compiled versions for each OS but their functionality
is not guaranteed. The Elastix version used is 4.9 with OpenCL support.
Additionally Elastix needs to be compiles with the PCA metrics, all DTI
related parameters and all affine related parameters.

Although cross platform compatibility is provided I have only limited
options for testing so if any issues arise please let me know.

------------------------------------------------------------------------

## Demonstrations

The release contains a zip file
[DemoAndTest.zip](https://github.com/mfroeling/QMRITools/releases/download/2.0/DemoAndTest.zip)
in which there is a file `demo.nb`, a folder `DemoData` and a folder
`Testing`. To have a global overview of the functionality of the toolbox
you can download this folder and run the `demo.nb`. By default the
`demo.nb` looks for the folders `DemoData` and `Testing` in the same
folder as the notebook.

In the first section of the demo notebook the toolbox is loaded and two
tests are performed. The first test is to check of all files that are
needed to run the toolbox are present. The second test runs after the
toolbox is loaded and checks if all the functions and their options that
are defined are correct.

------------------------------------------------------------------------

## Documentation

Documentation of all functions and their options is fully integrated in
the Mathematica documentation. The toolbox always works within the
latest version of Mathematica and does not support any backward
compatibility. After the toolbox is installed correctly it should show
up as a package in the Mathematica add-ons.

![QMRITools
package](https://github.com/mfroeling/QMRITools/blob/master/docs/images/addons.PNG)

All code and documentation is maintained and uploaded to github using
[Workbench](https://www.wolfram.com/workbench/). An online version of
the full documentation can be found
[here](https://mfroeling.github.io/QMRITools/htmldoc/guide/QMRITools.html).

![Guides
QMRITools](https://github.com/mfroeling/QMRITools/blob/master/docs/images/Guide.PNG)

------------------------------------------------------------------------

## Using the toolbox

The toolbox can be loaded by using: `` <<QMRITools` `` If you want to
monitor the package loading you can use:
`` QMRITools`$Verbose = True; <<QMRITools` ``

A list of all QMRITools packages is generated by

    QMRIToolsPackages[]

A list of all DTITools functions or functions per toolbox is generated
by

    QMRIToolsFunctions[]
    QMRIToolsFunctions["toolboxname"]

To print the documentation of all functions use

    QMRIToolsFuncPrint[]
    QMRIToolsFuncPrint["toolboxname"]

A list off all functions and their help can be found in
`All-Functions.nb`, which is alos availible as a [pdf file](https://github.com/mfroeling/QMRITools/releases/download/2.0/All-Functions.pdf).

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

The toolbox contains over 250 Functions and options of processing and
analyzing data. A summary of the core functionality is listed below.

![Overview](https://github.com/mfroeling/QMRITools/blob/master/docs/images/overview.png)

-   **Diffusion Analysis**
    -   Signal drift correction
    -   LLS, WLLS and iWLLS methods
    -   REKINDLE outlier detection
    -   IVIM fitting (fixed parameters, back-projection and Bayesian
        fitting)
    -   Parameter fitting using histogram analysis
    -   Joining and sorting of multiple series of the same volume
    -   Joining multiple stacks with slice overlap into one stack
-   **Diffusion Gradients optimization**
    -   Single and multi shell
    -   Rotating and correcting Bmatrix
    -   Actual b-value estimation by gradient sequence integration
    -   Gradient visualization
-   **Noise suppression**
    -   LMMSE noise suppression
    -   PCA noise suppression based on ramom matrix theory
    -   Anisotropic tensor smoothing using diffusion filter.
-   **Importing and Exporting**
    -   Dicom data (classing and enhanced file format)
    -   Nifti data (.nii and .img .hdr, supports .gz files)
    -   Compatible with ExplorDTI and Viste for fiber tractography
-   **Data visualization**
    -   2D 3D and 4D viewer
    -   Multiple images: Transparent overlay, difference and, checkboard overlays
    -   Legend bars and image labels
    -   Saving to pdf, jpg, animated gif and movie

![PlotData](https://github.com/mfroeling/QMRITools/blob/master/docs/images/visualization.PNG)

-   **Masking**
    -   Automate and threshold masking
    -   Extracting parameters form masks
    -   Smoothing masks
    -   Smoothing muscle segmentation
-   **Motion and distortion correction (Registration using elastix)**
    -   Rigid, affine, b-spline and cyclic registration
    -   nD to nD registration
    -   Automated series processing
    -   Slice to slice motion correction of 3D and 4D data
-   **Dixon Reconstruction**
    -   B0 phase unwrapping
    -   DIXON iDEAL reconstruction with T2star
-   **Relaxometry fitting**
    -   T2 fitting
    -   T1rho fitting
    -   Tri Exponential T2 fitting
    -   EPG based T2 fitting with slice profile
-   **Reconstruction Tools**
    -   Basic algorithms for complex coil combination including RSS, Roemer, WSVD.
    -   Simple Image reconstruction (MS-2D and 3D reconstruction) 
    -   CSI data reconstruction
    -   reading Philips *.list and *.data files.
    -   Hamming filtering.
-   **Simulation Framework**
    -   Diffuison tensor simulation and analysis
    -   Bloch and EPG simulations
    -   Cardiac DTI models (fiber architecture)
-   **Spectra processing and fitting**
    -   Fitting spectra using a set of basis functions. 
    -   Phase correction and line-width optimization.
    -   Basic plotting and visualization.
    -   Automatic 0th and 1st order phase correction using Henkel matrix SVD. 
    -   baseline correction
    -   CSI data visualization. 
-   **Cardiac Diffusion analysis**
    -   Breathing motion correction
    -   Corrupted slice rejection
    -   Local myocardial coordinate system calculation
    -   helix angle and fiber architecture matrix
    -   AHA 17 parameter description
    -   Transmural parameter description

**Under Construction**

-   **Tagging analysis**
    -   Getting displacement fields from tagging data (grid and line tagging)
    -   Strain, trosion and rotation from cardiac data

------------------------------------------------------------------------

## Toolboxes

### CardiacTools

A collection of tools to analyze cardiac data. The main features are
cardiac shape analysis which allows defining the hard in a local
myocardial coordinate system which allows quantifying and analyzing
data. When the cardiac geometry is known there are functions to analyze
qMRI parameters in the AH17 model (Cerqueira et al. 2002) or perform
transmural sampling of qMRI parameters. Most of the functionality is
demonstrated in the `demo.nb`.

### CoilTools

A collection of tools to evaluate complex multi-coil data. The functions
are specific for analysis of multi-coil magnitude and noise data which
allows quantifying per channel SNR. Furthermore, if complex coil
sensitivity maps are available it allows performing SENSE g-factor maps
simulations.  
This toolbox is not demonstrated in the `demo.nb`.

### DenoiseTools

The toobox provides two algorithms that allow denoising of DWI data. The
first is based on and LMMSE framework (Aja-Fernandez et al. 2008) and
the second is based on a random matrix theory and Principal component
analysis framework (Veraart, Fieremans, and Novikov 2016; Veraart et al.
2016). Furthermore, it provides an anisotropic filter for denoising the
estimated diffusion tensor which provides more reliable fiber
orientation analysis (Lee, Chung, and Alexander 2006). Most of the
functionality is demonstrated in the `demo.nb`.

### DixonTools

An IDEAL based Dixon reconstruction algorithm (Reeder et al. 2005; Yu et
al. 2008). The method provides multi-peak fitting B0 field and T2-
correction. The toolbox also provides a function for unwrapping phase
data in 2D and 3D based on a best path method (Abdul-Rahman et al. 2007;
Herraez et al. 2002). It also contains a function that allows simulating
gradient echo Dixon data. Most of the functionality is demonstrated in
the `demo.nb`.

![IDEAL based Dixon reconstruction: fitted fat fractions as a function
of the imposed fat fraction, SNR and B0 field
offset.](https://github.com/mfroeling/QMRITools/blob/master/docs/images/dixon.png)

### ElastixTools

A wrapper that calls the Elastix registration framework (Klein et al.
2010; Shamonin 2013). The toolbox determines what registration or
transformations need to be performed, exports the related data to a temp
folder and calls an automatically generated command line script that
performs the registration. After registration is completed the data is
again loaded into Mathematica. Most of the functionality is demonstrated
in the `demo.nb`.

### GeneralTools

This toolbox provides core functions used in many other functions and
features. The functions comprise amongst others: data cropping,
mathematical and statistical operators that ignore zero values, and data
rescaling, transformation and padding. Most of the functionality is
demonstrated in the `demo.nb`.

### GradientTools

The main feature is an algorithm that uses static repulsion (Jones,
Horsfield, and Simmons 1999; Froeling et al. 2017) to generate
homogeneously distributed gradient directions for DWI experiments. It
also provides functions to convert bval and bvec files to bmatrix and
vice versa. Most of the functionality is demonstrated in the `demo.nb`.

![The graphical user interface of the gradient generation
tool.](https://github.com/mfroeling/QMRITools/blob/master/docs/images/gradients.png)

### ImportTools

Allows importing DCM data or DCM header attributes. These functions are
rarely used since the toolbox mostly uses the NIfTY data format and
provides tools to convert DCM to NIfTI via
[dcm2niix](https://github.com/rordenlab/dcm2niix). This toolbox is not
demonstrated in the `demo.nb`.

### IVIMTools

The toolbox includes functions to perform IVIM fitting of DWI data.
There are two main functions: non linear fitting and Bayesian fitting
(Orton et al. 2014). Some of the functionality is demonstrated in the
`demo.nb`.

![Visualization of IVIM
fitting.](https://github.com/mfroeling/QMRITools/blob/master/docs/images/ivim.png)

### JcouplingTools

A toolbox that allows simulation of NMR spectra using Hamiltonians based
on methods from [FID-A](https://github.com/CIC-methods/FID-A). It allows
simulating large spin systems (Castillo, Patiny, and Wist 2011) and was
mainly implemented to investigate fat spectra in TSE (Stokes et al.
2013). Most of the functionality is demonstrated in the `demo.nb`.

### MaskingTools

Tools for masking and homogenization of data. It provides functions for
smoothing cutting and merging masks and functions for the evaluation of
data within masks. Most of the functionality is demonstrated in the
`demo.nb`.

### NiftiTools

Import and export of the NIfTI file format. Part of the code is based on
previously implemented
[nii-converter](https://github.com/tomdelahaije/nifti-converter). For
converting DICOM data to the NIfTI file format the toolbox uses
[dcm2niix](https://github.com/rordenlab/dcm2niix/releases). It also
provides some specialized NIfTI import functions for specific
experiments which are probably not generalizable. Most of the
functionality is demonstrated in the `demo.nb`.

### PhysiologyTools

Functions for importing and analyzing Philips physiology logging and
RespirAct trace files. The functions are rarely used and not well
supported. This toolbox is not demonstrated in the `demo.nb`.

### PlottingTools

A variety of functions for visualization of various data types. The main
functions are ‘PlotData’ and ‘PlotData3D’ which allow viewing 2D, 3D and
4D data. Most of the functionality is demonstrated in the `demo.nb`.

### Reconstruction Tools

A variety of function for raw MRI data reconstruction. The main goal was to 
create a set of functions that allow for the reconstruction of multi coil 3D CSI 
data and and low SNR 31P imaging data. For this toolbox there is no demo.  

### ProcessingTools

The toolbox comprises a variety of functions that allow data
manipulation and analysis. The main functions allow joining multiple
data sets into one continuous data set (Froeling et al. 2015) or to
split data of two legs into two separate data-sets. Furthermore, it
contains a collection of functions for data evaluation and analysis.
Most of the functionality is demonstrated in the `demo.nb`.

### RelaxometryTools

A collection of tools to fit T2, T2\*, T1rho and T1 relaxometry data.
The main function of this toolbox is an extended phase graph (EPG)
(Weigel 2015) method for multi-compartment T2 fitting of multi-echo spin
echo data (Marty et al. 2016). Therefore it provides functions to
simulate and evaluate EPG. Some of the functionality is demonstrated in
the `demo.nb`.

![Demonstration of EPG based T2 fitting: the fitted water T2 relaxation
as a function of B1, SNR and fat
fraction.](https://github.com/mfroeling/QMRITools/blob/master/docs/images/epg-t2.png)

### SimulationTools

The main purpose of this toolbox is to simulate DTI based DWI data and
contains some functions to easily perform analysis of the fit results of
the simulated signals (Froeling et al. 2013). Some of the functionality
is demonstrated in the `demo.nb`.

### SpectroTools

The main purpose of this toolbox is to process and visualize spectra data
and allows to fit spectra using simulated basis spectra. Some of the functionality
is demonstrated in the `demo.nb`.

![Comparison of fitted and measured 31P spectra of muscle.](https://github.com/mfroeling/QMRITools/blob/master/docs/images/spectra_fit.gif)

![Resulting basis spectra of a fit of 31P spectra of muacle.](https://github.com/mfroeling/QMRITools/blob/master/docs/images/spectra_fit.jpg)

### TensorTools

The original toolbox where the project started. The main functions in
this toolbox are to fit and evaluate the diffusion tensor model. Various
fitting methods are implemented (e.g. LLS, NLS, WLLS, and iWLLS). The
default method is an iterative weighted linear least squares approach
(Veraart et al. 2013). The tensor fitting also includes outlier
detections using REKINDLE (Tax et al. 2015) and data preparation
includes drift correction (Vos et al. 2017). Most of the functionality
is demonstrated in the `demo.nb`.

![MD and FA as a function of SNR and fat fraction. Results are from
simulated data using an iWLLS algorithm with outlier
rejection.](https://github.com/mfroeling/QMRITools/blob/master/docs/images/dti.png)

### VisteTools

Import and export functions for tensor data which can be used in the
[vIST/e](http://bmia.bmt.tue.nl/software/viste/) tractography tool. None
of the functionality is demonstrated in the `demo.nb`.

------------------------------------------------------------------------

![Full leg diffusion tensor fiber
tractography](https://github.com/mfroeling/QMRITools/blob/master/docs/images/animation-small.gif)

------------------------------------------------------------------------

## License

https://opensource.org/licenses/BSD-3-Clause

Note that restrictions imposed by these patents (and possibly others)
exist independently of and may be in conflict with the freedoms granted
in BSD-3-Clause license, which refers to copyright of the program, not
patents for any methods that it implements. Both copyright and patent
law must be obeyed to legally use and redistribute this program and it
is not the purpose of this license to induce you to infringe any patents
or other property right claims or to contest validity of any such
claims. If you redistribute or use the program, then this license merely
protects you from committing copyright infringement. It does not protect
you from committing patent infringement. So, before you do anything with
this program, make sure that you have permission to do so not merely in
terms of copyright, but also in terms of patent law.

Some code in the NiiTools packages was based on
https://github.com/tomdelahaije/nifti-converter

## References

Abdul-Rahman, Hussein S., Munther A. Gdeisat, David R. Burton, Michael
J. Lalor, Francis Lilley, and Christopher J. Moore. 2007. “Fast and
robust three-dimensional best path phase unwrapping algorithm.” *Applied
Optics* 46 (26): 6623. <https://doi.org/10.1364/AO.46.006623>.

Aja-Fernandez, Santiago, Marc Niethammer, Marek Kubicki, Martha E.
Shenton, and Carl Fredrik Westin. 2008. “Restoration of DWI data using a
rician LMMSE estimator.” *IEEE Transactions on Medical Imaging* 27 (10):
1389–1403. <https://doi.org/10.1109/TMI.2008.920609>.

Castillo, Andrés M., Luc Patiny, and Julien Wist. 2011. “Fast and
accurate algorithm for the simulation of NMR spectra of large spin
systems.” *Journal of Magnetic Resonance* 209 (2). Academic Press:
123–30. <https://doi.org/10.1016/j.jmr.2010.12.008>.

Cerqueira, Manuel D., Neil J. Weissman, Vasken Dilsizian, Alice K.
Jacobs, Sanjiv Kaul, Warren K. Laskey, Dudley J. Pennell, John A.
Rumberger, Thomas Ryan, and Mario S. Verani. 2002. “Standardized
myocardial sementation and nomenclature for tomographic imaging of the
heart: A Statement for Healthcare Professionals from the Cardiac Imaging
Committee of the Council on Clinical Cardiology of the American Heart
Association.” *Circulation* 105 (4). Lippincott Williams & Wilkins:
539–42. <https://doi.org/10.1161/hc0402.102975>.

Froeling, Martijn, Aart J. Nederveen, Dennis F. R. Heijtel, Arno
Lataster, Clemens Bos, Klaas Nicolay, Mario Maas, Maarten R. Drost, and
Gustav J. Strijkers. 2012. “Diffusion-tensor MRI reveals the complex
muscle architecture of the human forearm.” *Journal of Magnetic
Resonance Imaging* 36 (1). Wiley Subscription Services, Inc., A Wiley
Company: 237–48. <https://doi.org/10.1002/jmri.23608>.

Froeling, Martijn, Aart J. Nederveen, Klaas Nicolay, and Gustav J.
Strijkers. 2013. “DTI of human skeletal muscle: The effects of diffusion
encoding parameters, signal-to-noise ratio and T2 on tensor indices and
fiber tracts.” *NMR in Biomedicine* 26 (11): 1339–52.
<https://doi.org/10.1002/nbm.2959>.

Froeling, Martijn, Jos Oudeman, G. J. Gustav J. Strijkers, Mario Maas,
M. R. Maarten R. Drost, Klaas Nicolay, and Aart J. A. J. Nederveen.
2015. “Muscle Changes Detected with Diffusion-Tensor Imaging after
Long-Distance Running.” *Radiology* 274 (2): 548–62.
<https://doi.org/10.1148/radiol.14140702>.

Froeling, Martijn, Chantal M. W. Tax, Sjoerd B. Vos, Peter R. Luijten,
and Alexander Leemans. 2017. “MASSIVE brain dataset: Multiple
acquisitions for standardization of structural imaging validation and
evaluation.” *Magnetic Resonance in Medicine* 77 (5). Milan: 1797–1809.
<https://doi.org/10.1002/mrm.26259>.

Herraez, Miguel Arevallilo, David R. Burton, Michael J. Lalor, and
Munther A. Gdeisat. 2002. “Fast two-dimensional phase-unwrapping
algorithm based on sorting by reliability following a noncontinuous
path.” *Applied Optics* 41 (35): 7437.
<https://doi.org/10.1364/AO.41.007437>.

Jones, D. K., M. A. Horsfield, and A. Simmons. 1999. “Optimal strategies
for measuring diffusion in anisotropic systems by magnetic resonance
imaging.” *Magnetic Resonance in Medicine* 42 (3). Division of Medical
Physics, University of Leicester, Leicester Royal Infirmary, Leicester
LE1 5WW, United Kingdom.: 515–25.
[https://doi.org/10.1002/(SICI)1522-2594(199909)42:3&lt;515::AID-MRM14&gt;3.0.CO;2-Q](https://doi.org/10.1002/(SICI)1522-2594(199909)42:3<515::AID-MRM14>3.0.CO;2-Q).

Klein, Stefan, Marius Staring, Keelin Murphy, Max A. Viergever, and
Josien P. W. Pluim. 2010. “Elastix: A toolbox for intensity-based
medical image registration.” *IEEE Transactions on Medical Imaging* 29
(1): 196–205. <https://doi.org/10.1109/TMI.2009.2035616>.

Lee, Jee Eun, M. K. Chung, and A. L. Alexander. 2006. “Evaluation of
Anisotropic Filters for Diffusion Tensor Imaging.” In *3<sup>rd</sup>
Ieee International Symposium on Biomedical Imaging*, 77–80. IEEE.
<https://doi.org/10.1109/ISBI.2006.1624856>.

Marty, Benjamin, Pierre Yves Baudin, Harmen Reyngoudt, Noura Azzabou,
Ericky C. A. Araujo, Pierre G. Carlier, and Paulo L. de Sousa. 2016.
“Simultaneous muscle water T2and fat fraction mapping using transverse
relaxometry with stimulated echo compensation.” *NMR in Biomedicine* 29
(4): 431–43. <https://doi.org/10.1002/nbm.3459>.

Orton, Matthew R., David J. Collins, Dow-Mu Koh, and Martin O. Leach.
2014. “Improved intravoxel incoherent motion analysis of diffusion
weighted imaging by data driven Bayesian modeling.” *Magnetic Resonance
in Medicine* 71 (1): 411–20. <https://doi.org/10.1002/mrm.24649>.

Reeder, Scott B., Angel R. Pineda, Zhifei Wen, Ann Shimakawa, Huanzhou
Yu, Jean H. Brittain, Garry E. Gold, Christopher H. Beaulieu, and
Norbert T. Pelc. 2005. “Iterative decomposition of water and fat with
echo asymmetry and least-squares estimation (IDEAL): Application with
fast spin-echo imaging.” *Magnetic Resonance in Medicine* 54 (3):
636–44. <https://doi.org/10.1002/mrm.20624>.

Shamonin, Denis. 2013. “Fast parallel image registration on CPU and GPU
for diagnostic classification of Alzheimer’s disease.” *Frontiers in
Neuroinformatics* 7 (January): 50.
<https://doi.org/10.3389/fninf.2013.00050>.

Stokes, Ashley M., Yesu Feng, Tanya Mitropoulos, and Warren S. Warren.
2013. “Enhanced refocusing of fat signals using optimized multipulse
echo sequences.” *Magnetic Resonance in Medicine* 69 (4).
Wiley-Blackwell: 1044–55. <https://doi.org/10.1002/mrm.24340>.

Tax, Chantal M.W., Willem M. Otte, Max A. Viergever, Rick M. Dijkhuizen,
and Alexander Leemans. 2015. “REKINDLE: Robust Extraction of Kurtosis
INDices with Linear Estimation.” *Magnetic Resonance in Medicine* 73
(2): 794–808. <https://doi.org/10.1002/mrm.25165>.

Veraart, Jelle, Els Fieremans, and Dmitry S. Novikov. 2016. “Diffusion
MRI noise mapping using random matrix theory.” *Magnetic Resonance in
Medicine* 76 (5): 1582–93. <https://doi.org/10.1002/mrm.26059>.

Veraart, Jelle, Dmitry S. Novikov, Daan Christiaens, Benjamin Ades-aron,
Jan Sijbers, and Els Fieremans. 2016. “Denoising of diffusion MRI using
random matrix theory.” *NeuroImage* 142 (November). Elsevier Inc.:
394–406. <https://doi.org/10.1016/j.neuroimage.2016.08.016>.

Veraart, Jelle, Jan Sijbers, Stefan Sunaert, Alexander Leemans, and Ben
Jeurissen. 2013. “Weighted linear least squares estimation of diffusion
MRI parameters: Strengths, limitations, and pitfalls.” *NeuroImage* 81
(November). Elsevier Inc.: 335–46.
<https://doi.org/10.1016/j.neuroimage.2013.05.028>.

Vos, Sjoerd B., Chantal M. W. Tax, Peter R. Luijten, Sebastien Ourselin,
Alexander Leemans, and Martijn Froeling. 2017. “The importance of
correcting for signal drift in diffusion MRI.” *Magnetic Resonance in
Medicine* 77 (1): 285–99. <https://doi.org/10.1002/mrm.26124>.

Weigel, Matthias. 2015. “Extended phase graphs: Dephasing, RF pulses,
and echoes - pure and simple.” *Journal of Magnetic Resonance Imaging*
41 (2). Wiley-Blackwell: 266–95. <https://doi.org/10.1002/jmri.24619>.

Yu, Huanzhou, Ann Shimakawa, Charles A. McKenzie, Ethan Brodsky, Jean H.
Brittain, and Scott B. Reeder. 2008. “Multiecho water-fat separation and
simultaneous R\*2 estimation with multifrequency fat spectrum modeling.”
*Magnetic Resonance in Medicine* 60 (5): 1122–34.
<https://doi.org/10.1002/mrm.21737>.
