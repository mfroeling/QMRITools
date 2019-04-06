(* ::Package:: *)

(* ::Title:: *)
(*QMRITools init File*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


If[$VersionNumber != 11.3,
	CreateDialog[Column[{Style["
	Current Mathematica version is "<>ToString[$VersionNumber]<>"
	The toolbox is tested developed in 11.3.
	You need to update! (Or I am behind)
	Some functions wont work in older versions
	", TextAlignment -> Center], DefaultButton[], ""}, 
	Alignment -> Center], WindowTitle -> "Update!"];
];


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
QMRITools`$Verbose = True;


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
(*
(* ::Section:: *)
(*Functions*)


(*check for latest version*)
UpdateWarning[]:=If[$VersionNumber != 11.3,
	CreateDialog[Column[{Style["
	Current Mathematica version is "<>ToString[$VersionNumber]<>"
	The toolbox is tested developed in 11.3.
	You need to update! (Or I am behind)
	Some functions wont work in older versions
	", TextAlignment -> Center], DefaultButton[], ""}, 
	Alignment -> Center], WindowTitle -> "Update!"];
];


(*Fucntions to clear, load and protect package functions*)
ClearFunctions[pack_,subpack_,print_:False]:=Module[{packageName,packageSymbols,packageSymbolsG},
	If[print,Print["--------------------------------------"]];
	Quiet[
		If[print,Print["Removing all definitions of "<>#]];
		packageName = pack <> #;
		Get[packageName];
		packageSymbols =Names[packageName <> "*"];
		packageSymbolsG="Global`"<>#&/@packageSymbols;
		
		(*remove all global and private definitions*)
		Unprotect @@ packageSymbols;
		ClearAll @@ packageSymbols;
		Remove @@ packageSymbols;
		Unprotect @@ packageSymbolsG;
		ClearAll @@ packageSymbolsG;
		Remove @@ packageSymbolsG;

		]& /@ subpack;
];


LoadPackages[pack_,subpack_,print_:False]:=Module[{},
	If[print,Print["--------------------------------------"]];
	(
		If[print, Print["Loading all definitions of "<>#]];
		Get[pack<>#];
	)&/@subpack;
]


PrintPackages[context_,print_:False]:=Module[{},
	If[print,
		Print["--------------------------------------"];
		Print["All defined packages are: "];
		Print[context];
	];
]


ProtectFunctions[pack_,subpack_,print_:False]:=Module[{},
	If[print,Print["--------------------------------------"]];
	(
		If[print,Print["protecting all definitions of "<>#]];
		SetAttributes[#,{Protected, ReadProtected}]&/@ Names[pack <> # <> "*"];
		If[print,Print[Names[pack <> # <> "*"]]];
		If[print,Print["--------------------------------------"]];
	)& /@ subpack;
]


(* ::Section:: *)
(*Settings*)


(*Change Default settings*)
(*prevents the excessive use of memory with large data sets*)
$HistoryLength = 0;


(*add all mathematica packages to the context path*)
package = "QMRITools`";

subPackages = {
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


(*define all the toolbox contexts*)
System`$QMRIToolsContextPaths = (package <> # & /@ subPackages);
$ContextPath = Union[$ContextPath, System`$QMRIToolsContextPaths];


(*state if verbose is true to monitor initialization*)
QMRITools`verbose = False;


(* ::Section:: *)
(*Initialize all packages*)


(*print all the packages*)
PrintPackages[System`$QMRIToolsContextPaths, QMRITools`verbose];

(*check mathematica version*)
UpdateWarning[];

(*clear all definitions from the subPacakges*)
ClearFunctions[package,subPackages,QMRITools`verbose];

(*load all packages*)
LoadPackages[package,subPackages,QMRITools`verbose];

(*Protect functions*)
ProtectFunctions[package,subPackages,QMRITools`verbose];

*)