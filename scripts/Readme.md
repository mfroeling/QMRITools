# Wolfram engine and script

- Install the free Wolfram Engine as described here:
https://support.wolfram.com/45743

- Install WolframScript as described here
https://reference.wolfram.com/language/workflow/InstallWolframScript.html

# QMRITools 

Install the latest version of QMRITools for the worlfram enging using the `Install_QMRITools.wls` WolframScript.

`wolframscript -f "path to file\Install_QMRITools.wls"`

To install an other version of the QMRITools from the gihub release [pages](https://github.com/mfroeling/QMRITools/releases) use

`wolframscript -f "path to file\Install_QMRITools.wls" https://github.com/mfroeling/QMRITools/releases/download/3.17.0/QMRITools-x.xx.x.paclet`

For the situation where you run the scrip from its own folder with the test folder also there this will be script

`wolframscript -f "Install_QMRITools.wls" https://github.com/mfroeling/QMRITools/releases/download/4.x.x/QMRITools-4.x.x.paclet`

# Segmentation script

Run the segmentation script (networks are currelty for out-phase dixon data):

`wolframscript -f "path to file\Segment_Nii.wls" "file to be segmented.nii.gz" "output file.nii"`

It does not matter what part of the leg it is and if one or two legs are in the field of view. 
For the situation where you run the scrip from its own folder with the test folder also there this will be script.

`wolframscript -f "Segment_Nii.wls" "test data\test_up.nii.gz" "test data\out_up.nii"`
`wolframscript -f "Segment_Nii.wls" "test data\test_low.nii.gz" "test data\out_low.nii"`

<p align="center">
<img src="https://github.com/mfroeling/QMRITools/blob/master/docs/images/script.png"
alt="Example script for automate leg muscle segmentation."
title="Example script for automate leg muscle segmentation."  
width="80%" />
</p>

By default the CPU is used, if you want to switch to GPU open the script in any text editor and change "CPU" for "GPU".
