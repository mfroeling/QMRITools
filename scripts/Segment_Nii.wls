#!/usr/bin/env wolframscript
(* ::Package:: *)

dir = Directory[];
Print["Working from directory ", dir];

Print["Loading QMRITools"];
load = Get["QMRITools`"];

If[load === $Failed,
	Print["Could not load QMRITools"];
	,
	pack = Last@PacletFind["QMRITools"];
	Print["  - Using QMRITools ",pack["Location"]];

	{infile, outfile} = Rest @ $ScriptCommandLine;
	
	Print["Loading input data ", infile];
	in = FileExistsQ[infile];
	
	If[!in,
		Print["Could not find input file."];
		,
		{dat, vox} = ImportNii[infile];
		
		Print["  - Performing segmentation"];
		seg = SegmentData[dat, "Legs", TargetDevice->"CPU"];
		
		Print["  - Exporting segmentation to ", outfile];
		CopyFile[GetAssetLocation["MusclesLegLabels"], 
			ConvertExtension[outfile, ".txt"],
			OverwriteTarget->True];
		ExportNii[seg, vox, outfile];
	]
];
