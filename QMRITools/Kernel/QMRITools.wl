(* ::Package:: *)

(* ::Title:: *)
(*QMRITools init File*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(*set HistoryLength to 0 to prevent excessive memory used when working with large data*)
$HistoryLength = 0;


(* ::Section:: *)
(*Package Definitions*)


BeginPackage["QMRITools`"];


(*check mathematica version*)
If[$VersionNumber < 13.3, CreateDialog[Column[{Style[
"Current Mathematica version is "<>ToString[$VersionNumber]<>"
The toolbox is tested developed in 13.3+.
Some functions might not work in older versions"
	, TextAlignment -> Center], DefaultButton[], ""}, Alignment -> Center], WindowTitle -> "Update!"];
];


(*usage notes*)
QMRITools`$SubPackages::usage = "List of the sub packages in the toolbox.";
QMRITools`$Contexts::usage = "The package contexts needed for loading.";
QMRITools`$ContextsFunctions::usage = "The package contexts with the list of functions for each context.";
QMRITools`$Verbose::usage = "When set True, verbose loading is used.";
QMRITools`$InstalledVersion::usage = "The version number of the installed package.";

QMRITools`ElastixTools`$lastElastixTemp::usage = "$lastElastixTemp is the last temporary directory used by Elastix.";
QMRITools`ElastixTools`$debugElastix::usage = "$debugElastix is a debug flag for Elastix functionality.";
QMRITools`SegmentationTools`$debugUnet::usage = "$debugUnet is a debug flag for Unet functionality."; 
QMRITools`MuscleBidsTools`$debugBids::usage = "$debugBids is a debug flag for Bids functionality.";
QMRITools`DenoiseTools`$debugDenoise::usage = "$debugBids is a debug flag for Denoise functionality.";

QMRITools`PlottingTools`$plotOptions::usage = "$plotOptions is a list of options for plotting.";


(* ::Section:: *)
(*Package Variables*)


(*set the toolbox name*)


(*sub packages names*)
QMRITools`$SubPackages = {
	"ScientificColorData`",
	(*core packages that contain many functions for other toolboxes*)
	"LoggingTools`", "GeneralTools`", "MaskingTools`", "NiftiTools`",
	"ElastixTools`", "PlottingTools`", "MuscleBidsTools`", "NeuralNetworkTools`", 
	(*toolboxes for processing specific data types*)
	"DixonTools`", "IVIMTools`", "DenoiseTools`", "CardiacTools`",
	"RelaxometryTools`", "GradientTools`", "TensorTools`", 
	"JcouplingTools`", "SpectroTools`", "ReconstructionTools`",
	(*general processing tools with lots of dependency's*)
	"TractographyTools`", "ProcessingTools`", "FasciculationTools`",
	"SimulationTools`", "CoilTools`", "TaggingTools`", "SegmentationTools`" ,
	"ShapeTools`"
	(*legacy functions*)
	,"Legacy`"
};


(*define context and verbose*)
QMRITools`$Contexts = (Context[] <> # & /@ QMRITools`$SubPackages);
QMRITools`$Verbose = If[QMRITools`$Verbose===True, True, False];
QMRITools`$LoadedColor = If[QMRITools`$LoadedColor===True, True, False];
QMRITools`$InstalledVersion = First[PacletFind[StringDrop[Context[],-1]]]["Version"];

QMRITools`ElastixTools`$lastElastixTemp = "";
QMRITools`ElastixTools`$debugElastix = False;
QMRITools`SegmentationTools`$debugUnet = False;
QMRITools`MuscleBidsTools`$debugBids = False;
QMRITools`DenoiseTools`$debugDenoise = False;

(*load all the packages without error reporting such we can find the names of all the functions and options*)
Quiet[Get/@QMRITools`$Contexts];
QMRITools`$ContextsFunctions = {#, Names[# <> "*"]}& /@ QMRITools`$Contexts;
tempDir = StringDrop[GetAssetLocation["ColorData"], -4];

Begin["`Private`"];

End[];

EndPackage[];


(* ::Section:: *)
(*Package Loader*)


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


(*Destroy all functions defined in the sub packages*)
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


(*getting the color functions, prevents from reloading if kernel is not restarted*)
If[!QMRITools`$LoadedColor, 
	If[QMRITools`$Verbose, 
		Echo["--------------------------------------"];
		Echo["Loading color data"];
	];
	Get["QMRITools`ScientificColorData`"];
	QMRITools`ScientificColorData`ExtractColorData[tempDir];
	QMRITools`ScientificColorData`AddScientificColors[tempDir];
	ClearAll[tempDir];
	Remove[tempDir];
	QMRITools`$LoadedColor = True;
]


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


(*Protect definitions*)
Protect/@{QMRITools`$InstalledVersion, QMRITools`$SubPackages, QMRITools`$Contexts, QMRITools`$ContextsFunctions};
Unprotect/@{
	"QMRITools`ElastixTools`$lastElastixTemp",
	"QMRITools`ElastixTools`$debugElastix", 
	"QMRITools`SegmentationTools`$debugUnet", 
	"QMRITools`MuscleBidsTools`$debugBids",
	"QMRITools`DenoiseTools`$debugDenoise",
	"QMRITools`$Log",
	"QMRITools`$LogFile"
};
