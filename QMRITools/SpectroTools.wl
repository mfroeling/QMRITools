(* ::Package:: *)

(* ::Title:: *)
(*QMRITools SpectroTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`SpectroTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`SpectroTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


ReadjMRUI::usage = 
"ReadJMRUI[file] read a jMRUI spectrum file. 
Output is the {time, spec, {begintime, samplingInterval}}."




(* ::Subsection:: *)
(*Options*)


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*Import data*)


(* ::Subsubsection::Closed:: *)
(*ReadJMRUI file*)


SyntaxInformation[ReadjMRUI]={"ArgumentsPattern"->{_}}

ReadjMRUI[file_]:=Block[{imp,data,head,series,pts,spec,time},
imp=Import[file,"Data"];
data=Select[imp,AllTrue[#,NumericQ]&&#=!={}&];
head=Flatten[Select[imp,!AllTrue[#,NumericQ]&&#=!={}&]];
head=(StringTrim/@StringSplit[#<>" ",":"])&/@Select[head,StringContainsQ[#,":"]&];
head[[2;;,2]]=(ToExpression[StringReplace[#,"E"->" 10^" ]]&/@head[[2;;,2]])/.Null->0;
head=Thread[head[[All,1]]->head[[All,2]]];
series="DatasetsInFile"/.head;
pts="PointsInDataset"/.head;
data=Partition[data,pts];
spec=Reverse/@(data[[All,All,3]]+data[[All,All,4]]I);
time=(data[[All,All,1]]+data[[All,All,2]]I);
{time,spec,N@{"BeginTime","SamplingInterval"}/1000/.head}
]





(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
