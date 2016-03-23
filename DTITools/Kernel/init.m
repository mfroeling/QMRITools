(* ::DTITools:: *)

(* Mathematica Init File *)

(*initialization functions*)
ClearFunctions[pack_,subpack_]:=Module[{packageName,packageSymbols,packageSymbolsG},
	Quiet[
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

UpdateWarning[]:=If[$VersionNumber != 10.4,
 CreateDialog[Column[{Style["
      Current Mathematica version is 10.4
      You need to update!
      Some functions wont work in 10.3 or older versions
      ", TextAlignment -> Center], DefaultButton[], ""}, 
    Alignment -> Center], WindowTitle -> "Update!"];
 ];

LoadPackages[print_]:=(	If[print,Print[#]];	Get[package<>#];)&/@subPackages;


(*Change Default settings*)
$HistoryLength = 1;

(*add all mathematica packages to the context path*)
package= "DTITools`";

subPackages = {"CardiacTools`", "DenoiseTools`", "ElastixTools`", "ExportTools`", 
   "GeneralTools`", "GradientTools`", "ImportTools`", "IVIMTools`", 
   "ManipulationTools`", "MaskingTools`", "NiftiTools`", 
   "PhysiologyTools`", "PlottingTools`", "ProcessingTools`", 
   "RegistrationTools`", "SimulationTools`"};

(*define all the toolbox contexts*)
System`$DTIToolsContextPaths::usage = "$DTIToolsContextPaths lists all the diffusion packages"
System`$DTIToolsContextPaths = (package <> # & /@ subPackages);

(*clear all definitions from the subPacakges*)
ClearFunctions[package,subPackages];

(*check mathematica version*)
UpdateWarning[];

(*load all packages*)
LoadPackages[False];
LoadPackages[False];