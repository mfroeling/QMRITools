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

To see how to use the script run:

`wolframscript -f "path to file\Segment_Nii.wls" --h`

This will give the up to data help that might change over time, currently it is:

```
Usage:
  script.wls -i <input> -o <output> [-g Y|N] [-l "Legs"|"UpperLegs"|"LowerLegs"]
Flags:
  --i    Input NIfTI file (required)
  --o    Output NIfTI file (required)
  --g    Use GPU? (Y/N, default: N)
  --l    Region to segment: Legs (default), UpperLegs, LowerLegs, Shoulder
  --h, --help  Show this help message
```

Run the segmentation script (works on in or out phase data, raw echo data, water of fat reconstruction):

`wolframscript -f "path to file\Segment_Nii.wls" --i "file to be segmented.nii.gz" --o "output file.nii"`

It does not matter what part of the leg it is and if one or two legs are in the field of view. 
For the situation where you run the scrip from its own folder with the test folder also there this will be script.

`wolframscript -f "Segment_Nii.wls" --i "test data\test_up.nii.gz" --o "test data\out_up.nii"`
`wolframscript -f "Segment_Nii.wls" --i "test data\test_low.nii.gz" --o "test data\out_low.nii"`

<p align="center">
<img src="https://github.com/mfroeling/QMRITools/blob/master/docs/images/script.png"
alt="Example script for automate leg muscle segmentation."
title="Example script for automate leg muscle segmentation."  
width="80%" />
</p>
