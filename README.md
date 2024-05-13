# Welcome to QRMITools

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.595302.svg)](https://zenodo.org/doi/10.5281/zenodo.595302)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.01204/status.svg)](https://doi.org/10.21105/joss.01204)
[![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/dwyl/esta/issues)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fmfroeling%2FQMRITools&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)

[![MR-Hub](https://github.com/mfroeling/QMRITools/blob/master/docs/images/MR-Hub.png)](https://ismrm.github.io/mrhub/)
[![MRSHub](https://github.com/mfroeling/QMRITools/blob/master/docs/images/MRSHub.png)](https://mrshub.org/software_analysis/#QMRITools)
[![Open and Reproducible Musculoskeletal Imaging Research](https://github.com/mfroeling/QMRITools/blob/master/docs/images/ORMIR.png)](https://ormircommunity.github.io/packages.html#other-packages)
[![OpenSourceImaging](https://github.com/mfroeling/QMRITools/blob/master/docs/images/open_source_images.png)](https://www.opensourceimaging.org/project/qmritools-mathematica-toolbox-for-quantitative-mri-data/)

[![wolfram language](https://github.com/mfroeling/QMRITools/blob/master/docs/images/wolfram_language.png)](https://www.wolfram.com/language/)
[![wolfram workbench](https://github.com/mfroeling/QMRITools/blob/master/docs/images/wolfram_workbench.png)](https://www.wolfram.com/workbench/)
[![Visual studio code](https://github.com/mfroeling/QMRITools/blob/master/docs/images/visual-studio-code.png)](https://marketplace.visualstudio.com/items?itemName=WolframResearch.wolfram)
[![eclipse](https://github.com/mfroeling/QMRITools/blob/master/docs/images/eclipse.png)](https://www.eclipse.org/)
[![Wolfram Mathematica](https://github.com/mfroeling/QMRITools/blob/master/docs/images/wolfram_mathematica.png)](http://www.wolfram.com/mathematica/)
[![Wolfram Community](https://github.com/mfroeling/QMRITools/blob/master/docs/images/community.png)](https://community.wolfram.com/groups/-/m/t/1661539)

------------------------------------------------------------------------

## Content

- [Introduction](#introduction)
- [Installation](#installation)
- [Citing](#citing)
- [Documentation](#documentation)
- [External dependencies](#external-dependencies)
- [Toolboxes](#toolboxes)
- [License](#license)

------------------------------------------------------------------------

## Introduction

`QMRITools` is written in Mathematica and contains a collection of
tools and functions for processing quantitative MRI data.
The toolbox does not provide a GUI and its
primary goal is to allow for fast and batch data processing, and
facilitate development and prototyping of new functions. The core of the
toolbox contains various functions for data manipulation and
restructuring.

For more information [visit our website](https://www.qmritools.com/)

<p align="center">
<img src="https://github.com/mfroeling/QMRITools/blob/master/docs/images/both-small.gif"
alt="bilateral whole leg diffusion tensor imaging muscle fiber tractography"
title="bilateral whole leg diffusion tensor imaging muscle fiber tractography"  
width="40%" />
</p>

## Installation

The latest release can be found [here](https://github.com/mfroeling/QMRITools/releases).
The toolbox is best installed via the Mathematica paclet system. For more information [visit the website](https://www.qmritools.com/doc/instal/).

Automatic installation:

1. Download the `QMRITools-x.x.x.paclet`.
2. Install the paclet using `PacletInstall`.

`PacletInstall["xxx\\QMRITools-x.x.x.paclet"]`  

Or alternatively you can directly install it from the latest release page

`PacletInstall["https://github.com/mfroeling/QMRITools/releases/download/x.x.x/QMRITools-x.x.x.paclet"]`

<p align="center">
<img src="https://github.com/mfroeling/QMRITools/blob/master/docs/images/processing.png"
alt="Quantitative muscle MRI processing of diffusion tensor imaging, T2 mapping and water fat chemical shift imaging."
title="Quantitative muscle MRI processing of diffusion tensor imaging, T2 mapping and water fat chemical shift imaging."
width="70%" />
</p>

## Citing

When using the toolbox please cite one of the following references:

1. Froeling M: QMRTools: a Mathematica toolbox for quantitative MRI
  analysis. J Open Source Softw 2019; 4:1204.
  [link](https://joss.theoj.org/papers/ef8bfb6c31499845d353b6a5af0d6300)
2. Froeling M, et al.: Reproducibility of diffusion tensor imaging in
  human forearm muscles at 3.0 T in a clinical setting. Magn Reson Med
  2010; 64:1182-1190.
  [link](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.22477)
3. Froeling M, et al.: Diffusion-tensor MRI reveals the complex muscle
  architecture of the human forearm. J Magn Reson Imaging 2012;
  36:237-248.
  [link](https://onlinelibrary.wiley.com/doi/10.1002/jmri.23608)
4. Schlaffke et al.: Multi‐center evaluation of stability and reproducibility of
  quantitative MRI measures in healthy calf muscles; NMR Biomed. 2019;32:e4119
  [link](https://onlinelibrary.wiley.com/doi/full/10.1002/nbm.4119)

## Media and awards

- During the 2023 ISMRM in Toronto QMRITools was awarded received the “Best Open Source Tool Award” from the Quantitative MRI study group.
- If you want to learn more about the workings of QMRITools you can watch a live discussion with the Wolfram academic outreach team about <a href="https://www.youtube.com/live/wupxxiPJkxU?si=22BV_HSSa5u7Ds3D" target="_blank">QMRITools</a> and the role of computational Wolfram technology.
- A more in depth explanation of <a href="https://www.youtube.com/live/LVUBupORthA?si=UjoNpM2szsrgB7xx" target="_blank">the paclet functionality</a> was presented to the Wolfram R&D Team.
- QMRITools is build using Wolfram language for which it was awarded the <a href="https://www.wolfram.com/events/technology-conference/innovator-award/2023/martijn-froeling" target="_blank">Wolfram Innovator Award</a> in 2023 during the Wolfram Technology conference.

<p align="center">
<img src="https://github.com/mfroeling/QMRITools/raw/master/docs/images/ToolAward.png"
alt="Best Open Source Tool Award for quantitative MRI."
title="Best Open Source Tool Award for quantitative MRI."
width="70%" />
</p>

## Documentation

An online version of the full documentation can be found [here](https://www.qmritools.com/assets/htmldoc/html/guide/qmritools).

<p align="center">
<img src="https://github.com/mfroeling/QMRITools/blob/master/docs/images/addons.PNG" alt="QMRITools package add on"  width="70%" />
</p>

## External dependencies

Some functions of QMRITools call on external executables and software.
These executables need to be present in “QMRITools” and are included in
the release. If for any reason you want to use other (older/newer)
versions you can replace them but functionality is not guaranteed. For
the latest version of these tools and their user license please visit
their website.

- [dcm2niix](https://github.com/rordenlab/dcm2niix/)
- [Elastix](https://elastix.lumc.nl/)

## Toolboxes

QMRITools contains the following [toolboxes](https://www.qmritools.com/tool/):

- CardiacTools
- CoilTools
- DenoiseTools
- DixonTools
- ElastixTools
- FasciculationTools
- GeneralTools
- GradientTools
- ImportTools
- IVIMTools
- JcouplingTools
- LoggingTools
- MaskingTools
- MuscleBidsTools
- NiftiTools
- PhysiologyTools
- PlottingTools
- ProcessingTools
- ReconstructionTools
- RelaxometryTools
- SimulationTools
- SpectroTools
- TaggingTools
- TensorTools
- TractographyTools
- VisteTools

------------------------------------------------------------------------

## License

`https://opensource.org/licenses/BSD-3-Clause`

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

Some code in the NiiTools packages was based on `https://github.com/tomdelahaije/nifti-converter`
