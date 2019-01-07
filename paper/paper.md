---
title: 'DTITools: a collection of tools for processing and visualization of diffusion MRI data.'
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
``QMRITools`` is written in Mathematica using Workbench and Eclipse and contains a collection of tools and functions for processing quantitative MRI data. The toolbox does not provide a GUI and its primary goal is to allow for fast and batch data processing, and facilitate development and prototyping of new functions. The core of the toobox contains various functions for data manipulation and restructuring. The list of current functional toolboxes is:
- CardiacTools
- CoilTools
- DenoiseTools
- DixonTools
- ElastixTools
- GeneralTools
- GradientTools
- ImportTools
- IVIMTools
- JcouplingTools
- MaskingTools
- NiftiTools
- PhysiologyTools
- PlottingTools
- ProcessingTools
- RelaxometryTools
- SimulationTools
- VisteTools

The toolbox includes some [demo data](https://github.com/mfroeling/QMRITools/tree/master/testdata) which is used in the [demo file](https://github.com/mfroeling/QMRITools/blob/master/demo.nb) ``demo.nb``. In this notebook some of the functionality of the toolbox is demonstrated. A full list of functions and packages can be found in the [file](https://github.com/mfroeling/QMRITools/blob/master/QMRITools/All-Functions.nb) ``All-Functions.nb``. For all functions and toolboxes help files and guides are available and the documentation is build such that is incorporated in the Mathematica documentation.  
	 
The library of functions grows along with the research it is used for. The toolbox works along side other software packages (e.g. [vIST/e](http://bmia.bmt.tue.nl/software/viste/), and [ITKSnap](http://www.itksnap.org/pmwiki/pmwiki.php)) and for some functionality it calls external executables (e.g. [dcm2nii](https://www.nitrc.org/projects/dcm2nii/) and [Elastix](http://elastix.isi.uu.nl/)). The toolbox has been used is various studies (e.g. @Froeling2012, @Hooijmans2015, @Froeling2015).

The toolbox contains some basic functionality such as DICOM and Nifti import, and 2D,3D and 4D data visualization. 

![Overview](OverView.png)

The advanced features comprise data registration (@Klein2010, @Shamonin2013), Noise suppression (@Aja-Fernandez2008, @Veraart2015), Diffusion drift correction (@Vos2014), gradient direction optimization (@Froeling2016), simulation framework (@Froeling2013), EPG based T2 fitting (@Marty2016 @Weigel2015) and IDEAL Dixon reconstruction (@Reeder2005 @Herraez2002).

# References
