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
The toolbox is tested developed in 15.0+.
Some functions might not work in older versions"
	, TextAlignment -> Center], DefaultButton[], ""}, Alignment -> Center], WindowTitle -> "Update!"];
];


(*usage notes*)
QMRITools`$SubPackages::usage = "List of the sub packages in the toolbox.";
QMRITools`$Contexts::usage = "The package contexts needed for loading.";
QMRITools`$ContextsFunctions::usage = "The package contexts with the list of functions for each context.";
QMRITools`$Verbose::usage = "When set True, verbose loading is used.";
QMRITools`$Legacy::usage = "When set True, legacy functions are also loaded.";
QMRITools`$InstalledVersion::usage = "The version number of the installed package.";


(* ::Section:: *)
(*Package Variables*)

(*See what and how to load*)
QMRITools`$Verbose = If[QMRITools`$Verbose===True, True, False];
QMRITools`$Legacy = If[QMRITools`$Legacy===True, True, False];
QMRITools`$LoadedColor = If[QMRITools`$LoadedColor===True, True, False];
QMRITools`$Loaded = If[QMRITools`$Loaded===True, True, False];

(*sub packages names*)
QMRITools`$SubPackages = {
	"GeneralTools`", "ScientificColorData`",
	(*core packages that contain many functions for other toolboxes*)
	"LoggingTools`", "MaskingTools`", "NiftiTools`",
	"ElastixTools`", "PlottingTools`", "MuscleBidsTools`", "NeuralNetworkTools`", 
	(*toolboxes for processing specific data types*)
	"DixonTools`", "IVIMTools`", "DenoiseTools`", "CardiacTools`",
	"RelaxometryTools`", "GradientTools`", "TensorTools`", 
	"JcouplingTools`", "SpectroTools`", "ReconstructionTools`",
	(*general processing tools with lots of dependency's*)
	"TractographyTools`", "ProcessingTools`", "FasciculationTools`",
	"SimulationTools`", "CoilTools`", "TaggingTools`", "SegmentationTools`" ,
	"ShapeTools`", "AmaresTools`",
	(*legacy functions*)
	If[QMRITools`$Legacy, "Legacy`", Nothing]
};

(*define context*)
QMRITools`$Contexts = (Context[] <> # & /@ QMRITools`$SubPackages);
QMRITools`$InstalledVersion = First[PacletFind[StringDrop[Context[],-1]]]["Version"];

Begin["`Private`"];

End[];

EndPackage[];


(* ::Section:: *)
(*Package Loader*)


(*Echo the Toolbox content and version*)
If[QMRITools`$Verbose,
	Echo["--------------------------------------"];
	Echo["Version number "<>ToString[QMRITools`$InstalledVersion], "QMRITools"];
];


(*On a reload within the same kernel session (development, e.g. via QMRIToolsDev) all functions are already
loaded and protected, so they first need to be unprotected/cleared before they can be redefined.
On a first load in a fresh kernel (normal use) nothing is defined yet, so this whole step is skipped
and there is no need to load everything twice.*)
If[QMRITools`$Loaded,
	If[QMRITools`$Verbose,
		Echo["--------------------------------------"];
		Echo["Removing all local and global definitions of:"];
	];

	(*load all the packages without error reporting such we can find the names of all the functions and options*)
	Quiet[Get/@QMRITools`$Contexts];
	(
		With[{
				functions = Names[# <> "*"],
				global = Intersection[Names["Global`*"], "Global`" <> # & /@ Names[# <> "*"]]
			},

			If[QMRITools`$Verbose,
				Echo["", #];
				If[global=!={}, Echo[global]]
			];

			Unprotect @@ Join[functions, global];
			ClearAll @@ Join[functions, global];
			Remove @@ global;
		]
	) &/@ QMRITools`$Contexts;
];

(*Load and protect all the sub packages with error reporting*)
If[QMRITools`$Verbose,
	Echo["--------------------------------------"];
	Echo["Loading and protecting all definitions of:"];
];

(
	If[QMRITools`$Verbose, Echo["", #]];
	Get[#];

	(*getting the color functions, prevents from reloading if kernel is not restarted*)
	If[# == "QMRITools`ScientificColorData`" && !QMRITools`$LoadedColor,
		If[QMRITools`$Verbose,
			Echo["--------------------------------------"];
			Echo["Loading color data"];
		];

		With[{tempDir = StringDrop[QMRITools`GeneralTools`GetAssetLocation["ColorData"], -4]},
			QMRITools`ScientificColorData`ExtractColorData[tempDir];
			QMRITools`ScientificColorData`AddScientificColors[tempDir];
		];
		QMRITools`$LoadedColor = True;
	];

)& /@ QMRITools`$Contexts;

QMRITools`$Loaded = True;

QMRITools`$ContextsFunctions = {#, Names[# <> "*"]}& /@ QMRITools`$Contexts;
(SetAttributes[#, {Protected, ReadProtected}]& /@ Last[#])& /@ QMRITools`$ContextsFunctions;

If[QMRITools`$Verbose,
	Echo["--------------------------------------"];
	Echo["Loaded packages and functions: "];
	(
		Echo["with exposed functions and options:", First@#];
		Echo[Grid[Partition[Last@#, 4, 4 , 1, ""], Alignment->Left, ItemSize->18]];
	)&/@ QMRITools`$ContextsFunctions
];


(*Protect all definitions, then release those that can be edited*)
Protect/@{QMRITools`$InstalledVersion, QMRITools`$SubPackages, QMRITools`$Contexts, QMRITools`$ContextsFunctions};
Unprotect/@{
	"QMRITools`ElastixTools`$lastElastixTemp",
	"QMRITools`ElastixTools`$debugElastix", 

	"QMRITools`NeuralNetworkTools`$debugUnet", 

	"QMRITools`MuscleBidsTools`$debugBids",

	"QMRITools`DenoiseTools`$debugDenoise",
	
	"QMRITools`AmaresTools`$AmaresB1Regularization",
	"QMRITools`AmaresTools`$AmaresT1Values",

	"QMRITools`$Log",
	"QMRITools`$LogFile"
};
