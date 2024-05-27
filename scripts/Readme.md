# Wolfram engine and script


Install the Wolfram Engine

https://support.wolfram.com/45743

Install WolframScript

https://reference.wolfram.com/language/workflow/InstallWolframScript.html


# QMRITools 

Install QMRITools using WolframScript

wolframscript -f "path to file\Install_QMRITools.wls" https://github.com/mfroeling/QMRITools/releases/download/3.17.0/QMRITools-x.xx.x.paclet


# Segmentation script

Run the segmentation script:

wolframscript -f "path to file\Segment_Nii.wls" "file to be segmented.nii.gz" "output file.nii"