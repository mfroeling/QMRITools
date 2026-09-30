#!/usr/bin/env wolframscript
(* ::Package:: *)

dir = Directory[];
Print["\nWorking from directory "]
Print[">> ", dir];

Print["\nLoading QMRITools"];
dev = PacletFind["QMRIToolsDev"] =!= {};
t1 = First@AbsoluteTiming[load = Get[If[dev, "QMRIToolsDev`", "QMRITools`"]]];
Print[">> Loading time [s]: ", ToString[Round[t1, .1]]];

If[load === $Failed, Print["\nCould not load QMRITools"]; Quit[1]];

pack = SelectFirst[PacletFind["QMRITools"], StringStartsQ[FindFile["QMRITools`"], #["Location"], IgnoreCase->True]&];
Print[">> ", "Using QMRITools " <> pack["Version"] <> " at " <> pack["Location"] <> If[dev, " (dev)", ""]];

args = ParseCommandLine[Rest @ $ScriptCommandLine];

If[KeyExistsQ[args, "--h"] || KeyExistsQ[args, "--help"],
	Print["Usage:"];
	Print["  script.wls --i <input> [--o <output>] [--g Y|N] [--l <region>]"];
	Print["Flags:"];
	Print["  --i    Input NIfTI file (required)"];
	Print["  --o    Output NIfTI file (default: <input>_seg.nii(.gz))"];
	Print["  --g    Use GPU? (Y/N, default: N)"];
	Print["  --l    Region to segment (default: Body), one of:"];
	Print["         ", StringRiffle[Keys[$SegmentationLabels], ", "]];
	Print["  --h, --help  Show this help message"];
	Quit[]
];

infile = Lookup[args, "--i", Missing[]];
outfile = Lookup[args, "--o", Missing[]];
gpuUse = Lookup[args, "--g", "N"];
region = Lookup[args, "--l", "Body"];

gpu = Switch[gpuUse, "Y"|"y"|"True"|"true"|"Yes"|"yes", "GPU", _, "CPU"];

If[!StringQ[infile], Print["\nNo input file given."]; Quit[1]];
If[!KeyExistsQ[$SegmentationLabels, region],
	Print["\nUnknown region \"", region, "\", use one of: ", StringRiffle[Keys[$SegmentationLabels], ", "]];
	Quit[1]
];
If[!NiiFileExistQ[infile], Print["\nCould not find input file."]; Quit[1]];

If[!StringQ[outfile],
	outfile = StringReplace[infile, {".nii.gz"->"_seg.nii.gz", ".nii"->"_seg.nii"}];
	Print["\nSaving output to "];
	Print[">> ", outfile];
];

Print["\nImporting data "];
Print[">> ", infile];
{dat, vox} = ImportNii[infile];

Print["\nPerforming segmentation of region ", region];
seg = SegmentData[{dat, vox}, region, TargetDevice->gpu, Monitor->True];
If[seg === $Failed, Print["\nSegmentation failed."]; Quit[1]];

Print["\nExporting segmentation to "];
Print[">> ", outfile];

CopyFile[GetAssetLocation[$SegmentationLabels[region]],
	ConvertExtension[outfile, ".txt"],
	OverwriteTarget->True];
ExportNii[seg, vox, outfile];
