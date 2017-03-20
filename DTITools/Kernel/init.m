(* ::DTITools:: *)

(* Mathematica Init File *)

(*initialization functions*)
ClearFunctions[pack_,subpack_,print_:False]:=Module[{packageName,packageSymbols,packageSymbolsG},
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

UpdateWarning[]:=If[$VersionNumber != 11.,
 CreateDialog[Column[{Style["
      Current Mathematica version is 11.0
      You need to update!
      Some functions wont work in older versions
      ", TextAlignment -> Center], DefaultButton[], ""}, 
    Alignment -> Center], WindowTitle -> "Update!"];
 ];

LoadPackages[pack_,subpack_,print_:False]:=(
	If[print,Print["Loading all definitions of "<>#]];
	Get[package<>#];
)&/@subPackages;


(*Change Default settings*)
$HistoryLength = 1;

(*add all mathematica packages to the context path*)
package= "DTITools`";

subPackages = {"CardiacTools`", "DenoiseTools`", "ElastixTools`", "ExportTools`", 
   "GeneralTools`", "GradientTools`", "ImportTools`", "IVIMTools`", 
   "ManipulationTools`", "MaskingTools`", "NiftiTools`", 
   "PhysiologyTools`", "PlottingTools`", "ProcessingTools`", "RelaxometryTools`", "SimulationTools`"};

(*define all the toolbox contexts*)
System`$DTIToolsContextPaths::usage = "$DTIToolsContextPaths lists all the diffusion packages"
System`$DTIToolsContextPaths = (package <> # & /@ subPackages);

Needs["CCompilerDriver`"]
System`$DTIToolsCompiler =If[Length[CCompilers[]] > 0, "C", "WVM"];

$ContextPath = Union[$ContextPath, System`$DTIToolsContextPaths]

(*state if verbose is true to monitor initialization*)
DTITools`verbose = False;

(*clear all definitions from the subPacakges*)
ClearFunctions[package,subPackages,DTITools`verbose];

(*check mathematica version*)
UpdateWarning[];

(*load all packages*)
LoadPackages[package,subPackages,DTITools`verbose];
LoadPackages[package,subPackages,DTITools`verbose];