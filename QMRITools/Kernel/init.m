(* ::Package:: *)

(* ::Title:: *)
(*QMRITools init File*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


$HistoryLength = 0;

(*package naem*)
QMRITools`$Package = "QMRITools`";
QMRITools`$SubPackages = {
	(*core packages that contain functions for other toolboxes*)
	"GeneralTools`", "MaskingTools`", "NiftiTools`", "ElastixTools`", "PlottingTools`",
	(*toolboxes for processing specific data types*)
	"DixonTools`", "IVIMTools`", "DenoiseTools`", "CardiacTools`", 
	"RelaxometryTools`", "GradientTools`", "TensorTools`",
	"JcouplingTools`","SpectroTools`", "ReconstructionTools`",
	(*general processing tools with lots of dependancys*)
	"TractographyTools`", "VisteTools`", "ProcessingTools`", 
	"SimulationTools`", "PhysiologyTools`", "CoilTools`", 
	 (*legacy import functions*)
	"ImportTools`"
};


(*define context and verbose*)
QMRITools`$Contexts = (QMRITools`$Package <> # & /@ QMRITools`$SubPackages);
QMRITools`$Verbose = If[QMRITools`$Verbose===True, True, False];


(*print the contexts*)
If[QMRITools`$Verbose,
	Print["--------------------------------------"];
	Print["--------------------------------------"];
	Print["All defined packages to be loaded are: "];
	Print[Column@QMRITools`$Contexts];
];


(*quietly load all the packages such we can get the function names*)
Quiet[Get/@QMRITools`$Contexts];


(*Destroy all functions defined in the subpackages*)
If[QMRITools`$Verbose, 
	Print["--------------------------------------"];
	Print["--------------------------------------"];
	Print["Removing all Local and Global definitions existing in "<>QMRITools`$Package];
	Print["--------------------------------------"];
	Print["--------------------------------------"];
];
(
	(*function names in local an global context*)
	local = Names[# <> "*"];
	global = Intersection[Names["Global`*"], "Global`" <> # & /@ local];

	If[QMRITools`$Verbose, 
		Print["Removing all existing definitions in "<>#];
		Print["- Functions defined in package: \n", local];
		Print["- Functions defined in package existing in global:\n", global];
		Print["--------------------------------------"]
	];
	
	(*clear local*)
	Unprotect @@ local;
	ClearAll @@ local;
	
	(*clear global*)
	Unprotect @@ global;
	ClearAll @@ global;
	Remove @@ global;
	
	ClearAll[local,global]
) &/@ QMRITools`$Contexts


(*reload all the sub packages with error reporting*)
If[QMRITools`$Verbose, 
	Print["--------------------------------------"];
	Print["Loading and protecting all definitions existing in "<>QMRITools`$Package];
	Print["--------------------------------------"];
	Print["--------------------------------------"];
];
(
	(*load the package*)
	If[QMRITools`$Verbose, Print["Loading all definitions in "<>#]];
	Get[#];
	
	(*get the package functions*)
	local = Names[# <> "*"];
	
	If[QMRITools`$Verbose,
		Print["Protecting all definitions in "<>#];
		Print[local];
		Print["--------------------------------------"]
	];
	
	(*protect all definitions in package*)
	SetAttributes[#,{Protected, ReadProtected}]&/@ local;
	ClearAll[local]
)&/@QMRITools`$Contexts;	

(*finish*)
If[QMRITools`$Verbose,
	Print["--------------------------------------"];
	Print["Done loading "<>QMRITools`$Package];
	Print["--------------------------------------"];
	Print["--------------------------------------"];
];

(*check mathematica version*)
If[$VersionNumber < 12,
	CreateDialog[Column[{Style["
	Current Mathematica version is "<>ToString[$VersionNumber]<>"
	The toolbox is tested developed in 12.0.
	You need to update! (Or I am behind)
	Some functions wont work in older versions
	", TextAlignment -> Center], DefaultButton[], ""}, 
	Alignment -> Center], WindowTitle -> "Update!"];
];