---
title: 'QMRTools: a Mathematica toolbox for quantitative MRI analysis.'
tags:
  - Diffusion Tensor Imaging
  - Dixon Imaging
  - Extended phase graph
  - Data processing
  - Mathematica
  - Simulation
authors:
  - name: Martijn Froeling
    orcid: 0000-0003-3841-0497
    affiliation: 1
affiliations:
 - name: Department of Radiology, University Medical Center Utrecht, Utrecht, The Netherlands
   index: 1
date: 7 January 2019
bibliography: paper.bib
---

# Summary
``QMRITools`` is written in Mathematica using Wolfram Workbench and Eclipse and contains a collection of tools and functions for processing quantitative MRI data. The toolbox does not provide a GUI and its primary goal is to allow for fast and batch data processing, and facilitate development and prototyping of new functions. The core of the toolbox contains various functions for data manipulation and restructuring.

The toolbox was developed mostly in the context of quantitative muscle, nerve and cardiac magnetic resonance imaging. The library of functions grows along with the research it is used for and started as a toolbox to analyze DWI data of muscle (@Froeling2012). Although there exist many different software packages and code repositories for much of the functionality in this toolbox, I was in need of one that did all. Furthermore, most diffusion packages are optimized for brain analysis and provide limited flexibility. QMRITools works alongside other software packages (e.g. [vIST/e](http://bmia.bmt.tue.nl/software/viste/), and [ITKSnap](http://www.itksnap.org/pmwiki/pmwiki.php)) and for some functionality it calls external executables (e.g. [dcm2nii](https://www.nitrc.org/projects/dcm2nii/) and [Elastix](http://elastix.isi.uu.nl/)). The toolbox has been used in various studies (e.g. @Froeling2012, @Hooijmans2015, @Froeling2015).

![An overview of some of the functionality in the QMRITools Mathematica add-on.](overview.png)

The toolbox includes some [demo data](https://github.com/mfroeling/QMRITools/releases/download/2.0/DemoAndTest.zip) which is used in the demo file ``demo.nb``. In this notebook most of the functionality of the toolbox is demonstrated. A full list of functions and packages can be found in the [file](https://github.com/mfroeling/QMRITools/blob/master/QMRITools/All-Functions.nb) ``All-Functions.nb`` (also availible as [pdf](https://github.com/mfroeling/QMRITools/releases/download/2.0/All-Functions.pdf)). For all functions and toolboxes, help files and guides are available and the documentation is build such that is incorporated in the Wolfram documentation.  

The toolbox contains some basic functionality such as DICOM and Nifti import, and 2D, 3D and 4D data visualization. The advanced features comprise data registration (@Klein2010, @Shamonin2013), noise suppression (@Aja-Fernandez2008, @Veraart2015), diffusion drift correction (@Vos2014), gradient direction optimization (@Froeling2016), simulation framework (@Froeling2013), EPG based T2 fitting (@Marty2016 @Weigel2015) and IDEAL Dixon reconstruction (@Reeder2005 @Herraez2002). The current functional toolboxes with a short description are listed below.

## CardiacTools
A collection of tools to analyze cardiac data. The main features are cardiac shape analysis which allows defining the hard in a local myocardial coordinate system which allows quantifying and analyzing data. When the cardiac geometry is known there are functions to analyze qMRI parameters in the AH17 model (@Cerqueira2002) or perform transmural sampling of qMRI parameters. 
Most of the functionality is demonstrated in the `demo.nb`.

## CoilTools
A collection of tools to evaluate complex multi-coil data. The functions are specific for analysis of multi-coil magnitude and noise data which allows quantifying per channel SNR. Furthermore, if complex coil sensitivity maps are available it allows performing SENSE g-factor maps simulations.   
This toolbox is not demonstrated in the `demo.nb`.

## DenoiseTools
The toobox provides two algorithms that allow denoising of DWI data. The first is based on and LMMSE framework (@Aja-Fernandez2008) and the second is based on a random matrix theory and Principal component analysis framework (@Veraart2015, @Veraart2016). Furthermore, it provides an anisotropic filter for denoising the estimated diffusion tensor which provides more reliable fiber orientation analysis (@JeeEunLee).
Most of the functionality is demonstrated in the `demo.nb`.

## DixonTools
An IDEAL based Dixon reconstruction algorithm (@Reeder2005, @Yu2008). The method provides multi-peak fitting B0 field and T2* correction. The toolbox also provides a function for unwrapping phase data in 2D and 3D based on a best path method (@Abdul-Rahman2007, @Herraez2002). It also contains a function that allows simulating gradient echo Dixon data.
Most of the functionality is demonstrated in the `demo.nb`.

![IDEAL based Dixon reconstruction: fitted fat fractions as a function of the imposed fat fraction, SNR and B0 field offset.](dixon.png)

## ElastixTools
A wrapper that calls the Elastix registration framework (@Klein2010, @Shamonin2013). The toolbox determines what registration or transformations need to be performed, exports the related data to a temp folder and calls an automatically generated command line script that performs the registration. After registration is completed the data is again loaded into Mathematica.
Most of the functionality is demonstrated in the `demo.nb`.

## GeneralTools
This toolbox provides core functions used in many other functions and features. The functions comprise amongst others: data cropping, mathematical and statistical operators that ignore zero values, and data rescaling, transformation and padding.
Most of the functionality is demonstrated in the `demo.nb`.

## GradientTools
The main feature is an algorithm that uses static repulsion (@Jones1999, @Froeling2016) to generate homogeneously distributed gradient directions for DWI experiments. It also provides functions to convert bval and bvec files to bmatrix and vice versa.
Most of the functionality is demonstrated in the `demo.nb`.

![The graphical user interface of the gradient generation tool.](gradients.png)

## ImportTools
Allows importing DCM data or DCM header attributes. These functions are rarely used since the toolbox mostly uses the NIfTY data format and provides tools to convert DCM to NIfTI via [dcm2niix](https://github.com/rordenlab/dcm2niix).
This toolbox is not demonstrated in the `demo.nb`.

## IVIMTools
The toolbox includes functions to perform IVIM fitting of DWI data. There are two main functions: non linear fitting and Bayesian fitting (@Orton2014). 
Some of the functionality is demonstrated in the `demo.nb`.

![Visualization of IVIM fitting.](ivim.png)

## JcouplingTools
A toolbox that allows simulation of NMR spectra using Hamiltonians based on methods from [FID-A](https://github.com/CIC-methods/FID-A). It allows simulating large spin systems (@Castillo2011) and was mainly implemented to investigate fat spectra in TSE (@Stokes2013).
Most of the functionality is demonstrated in the `demo.nb`.

## MaskingTools
Tools for masking and homogenization of data. It provides functions for smoothing cutting and merging masks and functions for the evaluation of data within masks.
Most of the functionality is demonstrated in the `demo.nb`.

## NiftiTools
Import and export of the NIfTI file format. Part of the code is based on previously implemented [nii-converter](https://github.com/tomdelahaije/nifti-converter). For converting DICOM data to the NIfTI file format the toolbox uses [dcm2niix](https://github.com/rordenlab/dcm2niix/releases). It also provides some specialized NIfTI import functions for specific experiments which are probably not generalizable.
Most of the functionality is demonstrated in the `demo.nb`.

## PhysiologyTools
Functions for importing and analyzing Philips physiology logging and RespirAct trace files. The functions are rarely used and not well supported.
This toolbox is not demonstrated in the `demo.nb`.

## PlottingTools
A variety of functions for visualization of various data types. The main functions are 'PlotData' and 'PlotData3D' which allow viewing 2D, 3D and 4D data.
Most of the functionality is demonstrated in the `demo.nb`.

## ProcessingTools
The toolbox comprises a variety of functions that allow data manipulation and analysis. The main functions allow joining multiple data sets into one continuous data set (@froeling2015a) or to split data of two legs into two separate data-sets. Furthermore, it contains a collection of functions for data evaluation and analysis.
Most of the functionality is demonstrated in the `demo.nb`.

## RelaxometryTools
A collection of tools to fit T2, T2*, T1rho and T1 relaxometry data. The main function of this toolbox is an extended phase graph (EPG) (@Weigel2015) method for multi-compartment T2 fitting of multi-echo spin echo data (@Marty2016). Therefore it provides functions to simulate and evaluate EPG. 
Some of the functionality is demonstrated in the `demo.nb`.

![Demonstration of EPG based T2 fitting: the fitted water T2 relaxation as a function of B1, SNR and fat fraction.](epg-t2.png)

## SimulationTools
The main purpose of this toolbox is to simulate DTI based DWI data and contains some functions to easily perform analysis of the fit results of the simulated signals (@Froeling2013).
Some of the functionality is demonstrated in the `demo.nb`.

## TensorTools
The original toolbox where the project started. The main functions in this toolbox are to fit and evaluate the diffusion tensor model. Various fitting methods are implemented (e.g. LLS, NLS, WLLS, and iWLLS). The default method is an iterative weighted linear least squares approach (@Veraart2013). The tensor fitting also includes outlier detections using REKINDLE (@Tax2015) and data preparation includes drift correction (@Vos2014).
Most of the functionality is demonstrated in the `demo.nb`.

![MD and FA as a function of SNR and fat fraction. Results are from simulated data using an iWLLS algorithm with outlier rejection.](dti.png)

## VisteTools
Import and export functions for tensor data which can be used in the [vIST/e](http://bmia.bmt.tue.nl/software/viste/) tractography tool.
None of the functionality is demonstrated in the `demo.nb`.

# References
