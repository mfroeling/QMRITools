(* ::Package:: *)

(* ::Title:: *)
(*QMRITools LoggingTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`LoggingTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`LoggingTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


QMRITools`$Log::usage = 
"QMRITools`$Log is the current log. Is a list of strings."

ResetLog::usage = 
"ResetLog[] restes the log to {}."

ShowLog::usage = 
"ShowLog[] shows the log in a popup window.
ShowLog[False] shows the log in the notebook."

ExportLog::usage = 
"ExportLog[file] exports the log as a plain text to file."

ImportLog::usage = 
"ImportLog[file] imports the log as a list of string from a plain text file."

AddToLog::usage = 
"AddToLog[list] add the list to the log at level 1. All elements of the list are converted to strings and joined with spaces.
AddToLog[list, level] add the list to the log at level.
AddToLog[list, True] add the list to the log at level 1 with a timestamp.
AddToLog[list, True, level] specifies both the level and the timestamp.
AddToLog[list, level, True] specifies both the level and the timestamp."


(* ::Subsection:: *)
(*Options*)


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*QMRITools`$Log*)


QMRITools`$Log = {};


(* ::Subsection:: *)
(*ResetLog*)


SyntaxInformation[ResetLog] = {"ArgumentsPattern" -> {}};

ResetLog[] := QMRITools`$Log = {};


(* ::Subsection:: *)
(*ShowLog*)


SyntaxInformation[ShowLog] = {"ArgumentsPattern" -> {_.}};

ShowLog[] := ShowLog[True]

ShowLog[win_] := Block[{pane},
	pane = Pane[Dynamic@Column[QMRITools`$Log], {1000, 600}, Scrollbars -> {False, True}, ImageMargins -> 20];
	If[! win, pane,
		NotebookClose[logWindow];
		logWindow = CreateWindow[DialogNotebook@pane, WindowTitle -> "Logging window", Background -> White]
	];
];


(* ::Subsection:: *)
(*Import Export*)


(* ::Subsubsection:: *)
(*ExportLog*)


SyntaxInformation[ExportLog] = {"ArgumentsPattern" -> {_}};

ExportLog[file_] := Export[file, QMRITools`$Log, "Text"];


(* ::Subsubsection:: *)
(*ImportLog*)


SyntaxInformation[ImportLog] = {"ArgumentsPattern" -> {_}};

ImportLog[file_] := If[FileExistsQ[file], Import[file, "Lines"], ResetLog[]];


(* ::Subsection:: *)
(*ShowLog*)


SyntaxInformation[AddToLog] = {"ArgumentsPattern" -> {_,_.,_.}};

AddToLog[logAdd_?ListQ] := AddToLog[logAdd, 1, False]

AddToLog[logAdd_?ListQ, date_?BooleanQ] := AddToLog[logAdd, 1, date]

AddToLog[logAdd_?ListQ, lev_?IntegerQ] := AddToLog[logAdd, lev, False]

AddToLog[logAdd_?ListQ, date_?BooleanQ, lev_?IntegerQ] := AddToLog[logAdd, lev, date]

AddToLog[logAdd_?ListQ, lev_?IntegerQ, date_?BooleanQ] := AppendTo[QMRITools`$Log, 
	StringJoin[
		If[date, DateString[{"Day", "-", "Month", "-", "YearShort", " ","Time"}] <> " / ", "                  / "],
		If[lev == 0, "", StringJoin[ConstantArray[" ", 2 lev]] <> "- "],
		StringTrim[StringJoin[StringTrim[ToString[#]] <> " " & /@ logAdd]]
	]
];

AddToLog[logAdd_, a___] := AddToLog[{logAdd}, a]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
