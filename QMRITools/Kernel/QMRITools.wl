(* ::Package:: *)

(* ::Title:: *)
(*QMRITools init File*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(*set HistoryLength to 0 to prevent excessive memory used when working with large data*)
$HistoryLength = 0;
$ProgressReporting = False;

BeginPackage["QMRITools`"];


(*check mathematica version*)
If[$VersionNumber < 13.3, CreateDialog[Column[{Style[
"Current Mathematica version is "<>ToString[$VersionNumber]<>"
The toolbox is tested developed in 13.3+.
Some functions might not work in older versions"
	, TextAlignment -> Center], DefaultButton[], ""}, Alignment -> Center], WindowTitle -> "Update!"];
];


(*usage notes*)
QMRITools`$SubPackages::usage = "List of the subpackages in the toolbox.";
QMRITools`$Contexts::usage = "The package contexts needed for loading.";
QMRITools`$ContextsFunctions::usage = "The package contexts with the list of functions for each context.";
QMRITools`$Verbose::usage = "When set True, verbose loading is used.";
QMRITools`$InstalledVersion::usage = "The version number of the installed package.";
QMRITools`$Log::usage = "The logging file.";


(*subpackages names*)
QMRITools`$SubPackages = {
	(*core packages that contain many functions for other toolboxes*)
	"LoggingTools`", "GeneralTools`", "MaskingTools`", "NiftiTools`",
	 "ElastixTools`", "PlottingTools`", "MuscleBidsTools`",
	(*toolboxes for processing specific data types*)
	"DixonTools`", "IVIMTools`", "DenoiseTools`", "CardiacTools`",
	"RelaxometryTools`", "GradientTools`", "TensorTools`", 
	"JcouplingTools`", "SpectroTools`", "ReconstructionTools`",
	(*general processing tools with lots of dependancys*)
	"TractographyTools`", "ProcessingTools`", "FasciculationTools`",
	"SimulationTools`", "CoilTools`", "TaggingTools`", "SegmentationTools`",
	 (*legacy functions*)
	"Legacy`"
};


(*define context and verbose*)
QMRITools`$Contexts = (Context[] <> # & /@ QMRITools`$SubPackages);
QMRITools`$Verbose = If[QMRITools`$Verbose===True, True, False];
QMRITools`$InstalledVersion = First[PacletFind[StringDrop[Context[],-1]]]["Version"];


(*load all the packages without error reporting such we can find the names of all the functions and options*)
Quiet[Get/@QMRITools`$Contexts];
QMRITools`$ContextsFunctions = {#, Names[# <> "*"]}& /@ QMRITools`$Contexts;

Begin["`Private`"];

End[];

EndPackage[];


(*Load all subpackages*)


(*Echo the Toolbox content and version*)
If[QMRITools`$Verbose,
	Echo["--------------------------------------"];
	Echo["Version number "<>ToString[QMRITools`$InstalledVersion], "QMRITools"];
	Echo["--------------------------------------"];
	Echo["Defined packages and functions to be loaded are: "];
	(
		Echo["with exposed functions and options:", First@#];
		Echo[Grid[Partition[Last@#, 4, 4 , 1, ""], Alignment->Left, ItemSize->18]];
	)&/@ QMRITools`$ContextsFunctions
];


(*Destroy all functions defined in the subpackages*)
If[QMRITools`$Verbose, 
	Echo["--------------------------------------"];
	Echo["Removing all local and global definitions of:"];
];


With[{
		global = Intersection[Names["Global`*"], "Global`" <> # & /@ Last[#]]
	},
		
	If[QMRITools`$Verbose, 
		Echo["", First@#];
		If[global=!={}, Echo[global]]
	];

	Unprotect @@ Join[Last@#,global];
	ClearAll @@ Join[Last@#,global];
	Remove @@ global;
] &/@ QMRITools`$ContextsFunctions


(*Reload and protect all the sub packages with error reporting*)
If[QMRITools`$Verbose, 
	Echo["--------------------------------------"];
	Echo["Loading and protecting all definitions of:"];
];


Get["Developer`"];
(
	If[QMRITools`$Verbose, Echo["", First@#]];
	Get[First@#];
	SetAttributes[#, {Protected, ReadProtected}]& /@ Last[#]
)& /@ QMRITools`$ContextsFunctions;


Protect/@{QMRITools`$InstalledVersion, QMRITools`$SubPackages, QMRITools`$Contexts, QMRITools`$ContextsFunctions};
Unprotect/@{"QMRITools`ElastixTools`$debugElastix", "QMRITools`$Log"};