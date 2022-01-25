(* ::Package:: *)

(* ::Title:: *)
(*QMRITools init File*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)

(*set HistoryLength to 0 to prevent excessive memory used when working with large data*)
$HistoryLength = 0;

BeginPackage["QMRITools`"];

(*check mathematica version*)
If[$VersionNumber < 13,
	CreateDialog[Column[{Style["
	Current Mathematica version is "<>ToString[$VersionNumber]<>"
	The toolbox is tested developed in 13.0.
	You need to update! (Or I am behind)
	Some functions wont work in older versions
	", TextAlignment -> Center], DefaultButton[], ""}, 
	Alignment -> Center], WindowTitle -> "Update!"];
];

(*usage notes*)
QMRITools`$SubPackages::usage = "List of the subpackages.";
QMRITools`$Contexts::usage = "The package contexts.";
QMRITools`$ContextsFunctions::usage = "The package contexts with the list of functions.";
QMRITools`$Verbose::usage = "When set True, verbose loading is used.";
QMRITools`$InstalledVersion::usage = "The version number of the installed package.";

(*subpackages names*)
QMRITools`$SubPackages = {
	(*core packages that contain functions for other toolboxes*)
	"GeneralTools`", "MaskingTools`", "NiftiTools`", "ElastixTools`", "PlottingTools`",
	(*toolboxes for processing specific data types*)
	"DixonTools`", "IVIMTools`", "DenoiseTools`", "CardiacTools`", 
	"RelaxometryTools`", "GradientTools`", "TensorTools`",
	"JcouplingTools`","SpectroTools`", "ReconstructionTools`",
	(*general processing tools with lots of dependancys*)
	"TractographyTools`", "VisteTools`", "ProcessingTools`", 
	"SimulationTools`", "PhysiologyTools`", "CoilTools`", "TaggingTools`",
	 (*legacy import functions*)
	"ImportTools`"
};


(*define context and verbose*)
QMRITools`$Contexts = (Context[] <> # & /@ QMRITools`$SubPackages);
QMRITools`$Verbose = If[QMRITools`$Verbose===True, True, False];
QMRITools`$InstalledVersion = First[PacletFind[StringDrop[Context[],-1]]]["Version"];

(*load all the packages without error reporting such we can find the names of all the functions and options*)
Quiet[Get/@QMRITools`$Contexts];
QMRITools`$ContextsFunctions = {#, Names[# <> "*"]}& /@ QMRITools`$Contexts;

(*print the Toolbox content and version*)
If[QMRITools`$Verbose,
	Print["--------------------------------------"];
	Print["Loading ", Context[]," with version number ", QMRITools`$InstalledVersion];
	Print["--------------------------------------"];
	Print["Defined packages and functions to be loaded are: "];
	(
		Print["   - ", First@#, " with functions:"];
		Print[Last@#];
	)&/@ QMRITools`$ContextsFunctions
];


Begin["`Private`"];

End[];

EndPackage[];


(*Destroy all functions defined in the subpackages*)
If[QMRITools`$Verbose, 
	Print["--------------------------------------"];
	Print["Removing all local and global definitions of:"];
];

With[{
		global = Intersection[Names["Global`*"], "Global`" <> # & /@ Last[#]]
	},
		
	If[QMRITools`$Verbose, 
		Print["   - ", First@#];
		If[global=!={}, Print[global]]
	];

	Unprotect @@ Join[Last@#,global];
	ClearAll @@ Join[Last@#,global];
	Remove @@ global;
] &/@ QMRITools`$ContextsFunctions


(*Reload and protect all the sub packages with error reporting*)
If[QMRITools`$Verbose, 
	Print["--------------------------------------"];
	Print["Loading and protecting all definitions of:"];
];

(
	If[QMRITools`$Verbose, Print["   - ", First@#]];
	Get[First@#];
	SetAttributes[#, {Protected, ReadProtected}]& /@ Last[#]
)& /@ QMRITools`$ContextsFunctions;