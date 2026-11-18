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

`wolframscript -script "path to file\Segment_Nii.wls" --h`

Use `-script` and not `-f`, `wolframscript -f` does not pass arguments starting with `-` to the script.

This will give the up to data help that might change over time, currently it is:

```
Usage:
  script.wls --i <input> [--o <output>] [--g Y|N] [--l <region>]
Flags:
  --i    Input NIfTI file (required)
  --o    Output NIfTI file (default: <input>_seg.nii(.gz))
  --g    Use GPU? (Y/N, default: N)
  --l    Region to segment (default: Body), one of:
         <list of all regions, keys of $SegmentationLabels>
  --h, --help  Show this help message
```

The default region `Body` segments the full body with the full body muscle labels. The regions `Legs`, `LegsHip`, `UpperLegs`, `LowerLegs`, `Hip`, `Shoulder` and `Arm` give the older labels files, the other regions give the full body labels. The labels file is written next to the output file.

Run the segmentation script (works on in or out phase data, raw echo data, water of fat reconstruction):

`wolframscript -script "path to file\Segment_Nii.wls" --i "file to be segmented.nii.gz" --o "output file.nii"`

It does not matter what part of the leg it is and if one or two legs are in the field of view. 
For the situation where you run the scrip from its own folder with the test folder also there this will be script.

`wolframscript -script "Segment_Nii.wls" --i "test data\test_up.nii.gz" --o "test data\out_up.nii"`
`wolframscript -script "Segment_Nii.wls" --i "test data\test_low.nii.gz" --o "test data\out_low.nii"`

# Muscle-BIDS script

`MuscleBids.wls` runs the Muscle-BIDS pipeline of a study from the command line, the same steps as `MuscleBidsGUI` in
a notebook. The study folder holds the `config.json`, see the manual `BIDS-config.pdf` in the `BIDS Example` folder of
QMRITools. It needs QMRITools 5.3.1 or newer. To see the help run:

`wolframscript -script "path to file\MuscleBids.wls" --h`

Examples:

```
wolframscript -script MuscleBids.wls --dir "D:\study" --step All
wolframscript -script MuscleBids.wls --dir "D:\study" --step Process Merge --subs 001_1 002_1
wolframscript -script MuscleBids.wls --dir "D:\study" --step Tractography --method Segmentation --versioncheck
wolframscript -script MuscleBids.wls --dir "D:\study" --step All --report "D:\study\report.pdf"
wolframscript -script MuscleBids.wls --dir "D:\study" --check
wolframscript -script MuscleBids.wls --dir "D:\study" --report "D:\study\report.md"
```

| Flag | Meaning | Default |
| --- | --- | --- |
| `--dir <folder>` | study folder with the `config.json` (required) | |
| `--step <steps>` | `DcmToNii Convert Process Merge Segment Tractography Analysis` or `All`, space or comma separated. Steps always run in pipeline order | |
| `--subs <names>` | subjects to process: the dicom folder name (`001_1`) or `sub-001_ses-1`. A name without session selects session `001` | all |
| `--versioncheck` | redo data processed with an older QMRITools version | off |
| `--keepraw` | keep the raw nii of dcm2niix after Convert | removed |
| `--method <method>` | Tractography: `Full`, `Tractography` or `Segmentation` | `Full` |
| `--images <images>` | Analysis images: `All`, `Quantitative`, `Segmentation`, `Tractography` or `None` | `All` |
| `--check` | only show the problem summary of all steps, nothing is run | |
| `--report <file>` | write the problem report of all steps as `.md` or `.pdf`, after the steps or alone | |
| `--debug [tags]` | print debug output, default tag `Bids` | off |

The input is checked before anything runs: unknown flags, steps or values, a missing study folder or config, and
subject names that select no session stop the script. After the steps the script prints the number of problems per
step, the same as the Check tab of `MuscleBidsGUI` and `MuscleBidsReport`. Use `--report` to see which subjects and
which log files.

The progress lines show which subject and dataset is processed. The details are in the log files of each subject, as
in the notebook. Redirect the output to keep a record of an unattended run: `... > run.log 2>&1`.

Rerun after a crash or a stop: run the same command again. Process, Merge, Segment, Tractography and Analysis only
redo the datasets that did not finish. Convert has no check file and converts while raw files are present, so with
`--keepraw` it converts them again on every run. DcmToNii skips a session when its `raw` folder is not empty: after a
crash during DcmToNii, empty that folder to convert the session again.

Exit codes, for batch and cluster jobs:

| Code | Meaning |
| --- | --- |
| 0 | done without problems |
| 1 | invalid input, nothing was run |
| 2 | a step failed (aborted), the later steps were not run |
| 3 | done, but with messages, or with error, unreadable or missing check files in the steps that ran (outdated too with `--versioncheck`) |
