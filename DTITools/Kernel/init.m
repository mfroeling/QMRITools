(* ::Package:: *)

(* ::Title:: *)
(*DTITools init File*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Functions*)


UpdateWarning[]:=If[$VersionNumber != 11.3,
	CreateDialog[Column[{Style["
	Current Mathematica version is "<>ToString[$VersionNumber]<>"
	The toolbox is tested developed in 11.3.
	You need to update! (Or I am behind)
	Some functions wont work in older versions
	", TextAlignment -> Center], DefaultButton[], ""}, 
	Alignment -> Center], WindowTitle -> "Update!"];
];


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


LoadPackages[pack_,subpack_,print_:False,run_:1]:=Module[{},
	If[print,Print["--------------------------------------"]];
	(
		If[print,If[run==1,
			Print["Loading all definitions of "<>#<>", run 1."],
			Print["Loading all definitions of "<>#<>", run 2."]
		];];
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
$HistoryLength = 0; ParallelEvaluate[$HistoryLength = 0];


(*add all mathematica packages to the context path*)
package = "DTITools`";

subPackages = {
	"CardiacTools`", "DenoiseTools`", "DixonTools`",
	"ElastixTools`", "ExportTools`", "GeneralTools`", 
	"GradientTools`", "ImportTools`", "IVIMTools`", 
	"ManipulationTools`", "MaskingTools`", "NiftiTools`", 
	"PhysiologyTools`", "PlottingTools`", "ProcessingTools`", 
	"RelaxometryTools`", "SimulationTools`"};


(*define all the toolbox contexts*)
(*System`$DTIToolsContextPaths::usage = "$DTIToolsContextPaths lists all the diffusion packages"*)
System`$DTIToolsContextPaths = (package <> # & /@ subPackages);
$ContextPath = Union[$ContextPath, System`$DTIToolsContextPaths]


(*
Needs["CCompilerDriver`"]
System`$DTIToolsCompiler = If[Length[CCompilers[]] > 0, "WVM", "WVM"];
*)
System`$DTIToolsCompiler = "WVM";


(*state if verbose is true to monitor initialization*)
DTITools`verbose = False;


(* ::Section:: *)
(*Initialize all packages*)


If[DTITools`verbose,
Print["--------------------------------------"];
Print[System`$DTIToolsContextPaths];
];

(*check mathematica version*)
UpdateWarning[];
(*clear all definitions from the subPacakges*)
ClearFunctions[package,subPackages,DTITools`verbose];
(*load all packages*)
LoadPackages[package,subPackages,DTITools`verbose,1];
(*needs to be done twice else things dont work, don't understand why*)
LoadPackages[package,subPackages,DTITools`verbose,2];
(*Protect functions*)
ProtectFunctions[package,subPackages,DTITools`verbose];
