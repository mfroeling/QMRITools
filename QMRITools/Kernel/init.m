(* ::Package:: *)

(* ::Title:: *)
(*QMRITools init File*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


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


If[QMRITools`verbose,
Print["--------------------------------------"];
Print[System`$QMRIToolsContextPaths];
];


(*check mathematica version*)
UpdateWarning[];
(*clear all definitions from the subPacakges*)
ClearFunctions[package,subPackages,QMRITools`verbose];
(*load all packages*)
LoadPackages[package,subPackages,QMRITools`verbose];
(*Protect functions*)
ProtectFunctions[package,subPackages,QMRITools`verbose];
