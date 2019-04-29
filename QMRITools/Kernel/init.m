(* ::Package:: *)

(* ::Title:: *)
(*QMRITools init File*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(*package naem*)
QMRITools`$Package = "QMRITools`";
QMRITools`$SubPackages = {
	(*core packages that contain functions for other toolboxes*)
	"GeneralTools`", "MaskingTools`", "NiftiTools`", "ElastixTools`", "PlottingTools`",
	(*toolboxes for processing specific data types*)
	"DixonTools`", "IVIMTools`", "DenoiseTools`", "CardiacTools`", 
	"RelaxometryTools`", "GradientTools`", "TensorTools`", 
	"JcouplingTools`",
	(*general processing tools with lots of dependancys*)
	"VisteTools`", "ProcessingTools`", "SimulationTools`", "PhysiologyTools`", "CoilTools`", 
	 (*legacy import functions*)
	"ImportTools`"
};


(*define context and verbose*)
QMRITools`$Contexts = (QMRITools`$Package <> # & /@ QMRITools`$SubPackages);
QMRITools`$Verbose = If[QMRITools`$Verbose===True, True, False];


(*print the contexts*)
If[QMRITools`$Verbose,
	Print["--------------------------------------"];
	Print["All defined packages to be loaded are: "];
	Print[QMRITools`$Contexts];
];


(*load all the packages without error reporting such we can find the names*)
If[QMRITools`$Verbose, Print["--------------------------------------"]];
Quiet[Get/@QMRITools`$Contexts];


(*Destroy all functions defined in the subpackages*)
(
	If[QMRITools`$Verbose, 
		Print["Removing all definitions of "<>#];
		Print["- Package functions: \n", Names[# <> "*"]];
		Print["- Package functions in global:\n", Intersection[Names["Global`*"], "Global`" <> # & /@ Names[# <> "*"]]];
	];
	
	Unprotect @@ Names[# <> "*"];
	ClearAll @@ Names[# <> "*"];
	
	Unprotect @@ Intersection[Names["Global`*"], "Global`" <> # & /@ Names[# <> "*"]];
	ClearAll @@ Intersection[Names["Global`*"], "Global`" <> # & /@ Names[# <> "*"]];
	Remove @@ Intersection[Names["Global`*"], "Global`" <> # & /@ Names[# <> "*"]];
) &/@ QMRITools`$Contexts


(*reload all the sub packages with error reporting*)
If[QMRITools`$Verbose,Print["--------------------------------------"]];
(
	If[QMRITools`$Verbose, Print["Loading all definitions of "<>#]];
	Get[#];
)&/@QMRITools`$Contexts;	


(*protect all functions*)
If[QMRITools`$Verbose,Print["--------------------------------------"]];
(
	If[QMRITools`$Verbose,
		Print["protecting all definitions of "<>#];
		Print[Names[# <> "*"]];
		Print["--------------------------------------"]
	];
	
	SetAttributes[#,{Protected, ReadProtected}]&/@ Names[# <> "*"];
)& /@ QMRITools`$Contexts;


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