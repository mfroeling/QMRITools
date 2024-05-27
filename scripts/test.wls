#!/usr/bin/env wolframscript
(* ::Package:: *)

dir = Directory[];
Print["Working from directory ", dir];

Print["Loading QMRITools"];
load = Get["QMRITools`"];

SegmentDataGUI[]:=Block[{dat,vox,seg},
NotebookClose[segwindow];

diag=DialogNotebook[
status=TextCell@"";
{
TextCell["Please enter the paths for the input and output files:"],

Grid[{{
TextCell["Status: "],Dynamic[status]
},{
TextCell["Input File: "],
InputField[Dynamic[inputFile],String,FieldHint->"Enter input file path",FieldSize->{25,1}],
Button["Browse",inputFile=SystemDialogInput["FileOpen"],Method->"Queued"]
},{
TextCell["Output File: "],InputField[Dynamic[outputFile],String,FieldHint->"Enter output file path",FieldSize->{25,1}],
Button["Browse",
outputFile=SystemDialogInput["FileSave"],Method->"Queued"]
},{
TextCell["Segmentation type"],PopupMenu[Dynamic[option],{"Legs","LegBones"}]
}
},Alignment->Left],
Row[{
Button["Segment",
If[!NiiFileExistQ[inputFile],
MessageDialog["Input file could not be foud."]
,
status=TextCell@"Importing";
{dat,vox}=ImportNii[inputFile];
status=TextCell@"Segmenting Data";
seg = SegmentData[dat,"Legs", TargetDevice->"CPU"];
status=TextCell@"Exporting";
CopyFile[GetAssetLocation["MusclesLegLabels"], ConvertExtension[outputFile, ".txt"], OverwriteTarget->True];
ExportNii[seg, vox, outputFile];
status=Button["Got to "<>FileBaseName@outputFile,SystemOpen[DirectoryName@outputFile]];
]
,
Method->"Queued"
],
DefaultButton[],
CancelButton[]
}]
}
];

segwindow=CreateWindow[diag,WindowTitle->"Muscle segmentation",WindowSize->All];
];

SegmentDataGUI[]
