#!/usr/bin/env wolframscript
(* ::Package:: *)

dir = Directory[];
Print["\nWorking from directory "]
Echo[dir];

Print["\nLoading QMRITools"];
t1 = First@AbsoluteTiming[load = Get["QMRITools`"]];
Echo[Round[t1,.1], "Loading time [s]: "];

If[load === $Failed,
	Echo["Could not load QMRITools"];
	,
	pack = Last@PacletFind["QMRITools"];
	Echo["Using QMRITools ",pack["Location"]];

	args = ParseCommandLine[Rest @ $ScriptCommandLine];
	
	(* Help flag detection *)
	If[MemberQ[Keys[args], "-h"] || MemberQ[Keys[args], "--help"],
		Print["Usage:"];
		Print["  script.wls -i <input> -o <output> [-g Y|N] [-l \"Legs\"|\"Upperleg\"|\"Lowerlegs\"]"];
		Print["Flags:"];
		Print["  -i    Input NIfTI file (required)"];
		Print["  -o    Output NIfTI file (required)"];
		Print["  -g    Use GPU? (Y/N, default: N)"];
		Print["  -l    Region to segment: Legs (default), Upperleg, or Lowerlegs"];
		Print["  -h, --help  Show this help message"];
	];
	
	If[Keys[args] === {"-h"} || Keys[args] === {"--help"}, Quit[]];
	
	infile = Lookup[args, "-i", Missing[]];
	outfile = Lookup[args, "-o", Missing[]];
	gpuUse = Lookup[args, "-g", "N"];
	region = Lookup[args, "-l", "Legs"];
	
	gpu = Switch[gpuUse, "Y"|"y"|"True"|"true"|"Yes"|"yes", "GPU", _, "CPU"];
	If[infile===Missing[], Print["\nNo input file given."];Quit[1]];
	If[outfile===Missing[], outfile = StringReplace[infile,{".nii.gz"->"_seg.nii.gz",".nii"->"_seg.nii"}]];
	
	Print["\nLoading input data "];
	Echo[infile];
	
	in = NiiFileExistQ[infile];
	
	If[!in,
		Print["\nCould not find input file."];
		,
		Print["\nImporting and masking data"];
		{dat, vox} = ImportNii[infile];
		dat = MaskData[dat, Mask[NormalizeData[dat], 10, MaskSmoothing -> True, 
			MaskClosing->True, MaskDilation->1]];
	
		Print["\nPerforming segmentation"];
		seg = SegmentData[dat, region, TargetDevice->gpu, Monitor->True];
			
		Print["\nExporting segmentation to "]
		Echo[outfile];
		
		CopyFile[GetAssetLocation["MusclesLegLabels"], 
			ConvertExtension[outfile, ".txt"],
			OverwriteTarget->True];
		ExportNii[seg, vox, outfile];
	]
];
