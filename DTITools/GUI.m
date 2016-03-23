(* ::Package:: *)

(* ::Title:: *)
(*DTITools GradientTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling*)
(*m.froeling@gmail.com*)


BeginPackage["DTITools`GUI`", { 
	"DTITools`ExportTools`", 
	"DTITools`GradientTools`", 
	"DTITools`ImportTools`", 
	"DTITools`ManipulationTools`", 
	"DTITools`MaskingTools`", 
	"DTITools`PlottingTools`", 
	"DTITools`ProcessingTools`", 
	"DTITools`RegistrationTools`", 
	"DTITools`SimulationTools`"
	}]
(* Exported symbols added here with SymbolName::usage *)  

Unprotect @@ Names["DTITools`GUI`*"];
ClearAll @@ Names["DTITools`GUI`*"];



(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)

(*
StartGUI::usage = "StartGUI[] starts the graphical user interface for the DTITools package"
*)

(* ::Subsection:: *)
(*Options*)


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsubsection::Closed:: *)
(*StartGUI*)


StartGUI[]:=Module[{},
	CreatePalette[DynamicModule[{},
		"test"
		]
	]
]




(* ::Section:: *)
(*End Package*)


End[] 

SetAttributes[#,{Protected, ReadProtected}]& /@ Names["DTITools`GUI`*"];

EndPackage[]
