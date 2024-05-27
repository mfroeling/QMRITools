# Wolfram engine and script


Install the Wolfram Engine

https://support.wolfram.com/45743

Install WolframScript

https://reference.wolfram.com/language/workflow/InstallWolframScript.html


# QMRITools 

Install QMRITools using WolframScript

The latest release link can be found on the releases (page)[https://github.com/mfroeling/QMRITools/releases]

wolframscript -f "path to file\Install_QMRITools.wls" https://github.com/mfroeling/QMRITools/releases/download/3.17.0/QMRITools-x.xx.x.paclet


# Segmentation script

Run the segmentation script:

wolframscript -f "path to file\Segment_Nii.wls" "file to be segmented.nii.gz" "output file.nii"

For the situation where you run the scrip from its own folder with the test folder also there this will be script

wolframscript -f "Segment_Nii.wls" "test data\test.nii.gz" "test data\out.nii"

<p align="center">
<img src="https://github.com/mfroeling/QMRITools/blob/master/docs/images/script.png"
alt="Example script for automate leg muscle segmentation."
title="Example script for automate leg muscle segmentation."  
width="40%" />
</p>