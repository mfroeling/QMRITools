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

	{infile, outfile} = Rest @ $ScriptCommandLine;
	
	Print["\nLoading input data "];
	Echo[infile];
	
	in = NiiFileExistQ[infile];
	
	If[!in,
		Print["\nCould not find input file."];
		,
		{dat, vox} = ImportNii[infile];
	
		Print["\nPerforming segmentation"];
		seg = SegmentData[dat, "Legs", TargetDevice->"CPU", Monitor->True];
			
		Print["\nExporting segmentation to "]
		Echo[outfile];
		
		CopyFile[GetAssetLocation["MusclesLegLabels"], 
			ConvertExtension[outfile, ".txt"],
			OverwriteTarget->True];
		ExportNii[seg, vox, outfile];
	]
];
