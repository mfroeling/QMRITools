---
title: 'DTITools: a collection of tools for processing and visualization of diffusion MRI data.'
tags:
  - Diffusion Tensor Imaging
  - Data processing
  - Mathematica
  - Simulation
authors:
  - name: Martijn Froeling
    orcid: 0000-0003-3841-0497
    affiliation: "1"
affiliations:
 - name:  affiliation: Department of Radiology, University Medical Center Utrecht, Utrecht, The Netherlands
   index: 1
date: 02 June 2018
bibliography: paper.bib
---

# Summary
DTITools is written Mathematica and contains a collection of tools and functions primaraly designed for processing diffusion MRI data. The toolbox does not provide a GUI and its primary goal is to allow for fast development and prototyping of new functions and batch data processing. The library of functions grows along with the research it is used for. The toobox works along side other software packages (e.g. [vIST/e](http://bmia.bmt.tue.nl/software/viste/) and [ExploreDTI](http://www.exploredti.com/)) and for some functionality it calls external executables (e.g. [MRIcron](https://www.nitrc.org/projects/mricron/), [dcm2nii](https://www.nitrc.org/projects/dcm2nii/) and [Elastix](http://elastix.isi.uu.nl/)). The toobox has been used is various studies (e.g. Froeling et al. 2012, Hooijmans et al., 2015, Froeling et al. 2015).

The toobox contains, amongst others, DICOM and Nifti import, 2D,3D and 4D data visualization, data registration (S. Klein et al. 2010, Shamonin et al. 2014), Noise supression (Aja-Fernandez et al. 2008, Veraart et al.2015), Diffusion drift correction (Vos et al. 2016), gradient direction optimization (Froeling et al. 2016), simulation framework (Froeling et al. 2013), Bayesian IVIM fitting (Orton et al. 2014) and more.


# References

Froeling M, Strijkers GJ, Nederveen AJ, Luijten PR. Whole heart DTI using asymmetric bipolar diffusion gradients. J. Cardiovasc. Magn. Reson. 2015;17:P15. doi: 10.1186/1532-429X-17-S1-P15.

Hooijmans MT, Damon BM, Froeling M, Versluis MJ, Burakiewicz J, Verschuuren JJGM, Niks EH, Webb AG, Kan HE. Evaluation of skeletal muscle DTI in patients with duchenne muscular dystrophy. NMR Biomed. 2015;28:1589–1597. doi: 10.1002/nbm.3427.

Froeling M, Nederveen AJ, Heijtel DFR, Lataster A, Bos C, Nicolay K, Maas M, Drost MR, Strijkers GJ. Diffusion-tensor MRI reveals the complex muscle architecture of the human forearm. J. Magn. Reson. Imaging 2012;36:237–248. doi: 10.1002/jmri.23608.

Klein S, Staring M, Murphy K, Viergever MA, Pluim JPW. Elastix: A toolbox for intensity-based medical image registration. IEEE Trans. Med. Imaging 2010;29:196–205. doi: 10.1109/TMI.2009.2035616.

Shamonin DP, Bron EE, Lelieveldt BPF, Smits M, Klein S, Staring M. Fast parallel image registration on CPU and GPU for diagnostic classification of Alzheimer’s disease. Front. Neuroinform. 2013;7:50. doi: 10.3389/fninf.2013.00050.

Aja-Fernandez S, Niethammer M, Kubicki M, Shenton ME, Westin CF. Restoration of DWI data using a rician LMMSE estimator. IEEE Trans. Med. Imaging 2008;27:1389–1403. doi: 10.1109/TMI.2008.920609.

Veraart J, Fieremans E, Novikov DS. Diffusion MRI noise mapping using random matrix theory. Magn. Reson. Med. 2015;0:n/a-n/a. doi: 10.1002/mrm.26059.

Vos SB, Tax CMW, Luijten PR, Ourselin S, Leemans A, Froeling M. The Importance of Correcting for Signal Drift in Diffusion MRI. Magn. Reson. Med. 2016;0:1–15. doi: 10.1002/mrm.26124.

Froeling M, Tax CMW, Vos SB, Luijten PR, Leemans A. “MASSIVE” Brain Dataset: Multiple Acquisitions for Standardization of Structural Imaging Validation and Evaluation. Magn. Reson. Med. 2016. doi: 10.1002/mrm.26259.

Froeling M, Nederveen AJ, Nicolay K, Strijkers GJ. DTI of human skeletal muscle: The effects of diffusion encoding parameters, signal-to-noise ratio and T2 on tensor indices and fiber tracts. NMR Biomed. 2013;26:1339–1352. doi: 10.1002/nbm.2959.

Orton MR, Collins DJ, Koh D-M, Leach MO. Improved intravoxel incoherent motion analysis of diffusion weighted imaging by data driven Bayesian modeling. Magn. Reson. Med. 2014;71:411–20. doi: 10.1002/mrm.24649.