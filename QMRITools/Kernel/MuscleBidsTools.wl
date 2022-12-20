(* ::Package:: *)

(* ::Title:: *)
(*QMRITools MuscleBidsTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`MuscleBidsTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`MuscleBidsTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


ImportJSON::usage = 
"ImportJSON[file] impors a json file as rawJSON."

GetJSONPosition::usage = 
"GetJSONPosition[{json..}, {{key, value}..}] gets the position from a list of JSON association lists where keys have the given value.
GetJSONPosition[{json..}, {{key, value}..}, sortkey] same but finaly sorts the positions for the value of the sortkey."

MergeJSON::usage = 
"MergeJSON[{json..}] merges a list of JSON association lists where duplicate keys with same values are removed and duplicate keys with different values are merges."

AddToJson::usage = 
"AddToJson[json, <|key->value..|>] adds new keys and values to the JSON list where duplicte keys are eitehr removed or joined.
AddToJson[json, \"QMRITools\"] adds the QMRITools software version to the json.
"



(* ::Subsection:: *)
(*Options*)


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*JSON*)


(* ::Subsubsection:: *)
(*ImportJSON*)


ImportJSON[file_]:=Import[file,"RawJSON"]


(* ::Subsubsection:: *)
(*GetJSONPosition*)


GetJSONPosition[json_,selection_]:=GetJSONPosition[json,selection,""]

GetJSONPosition[json_,selection_,sort_]:=Block[{seli,self,list,key,val,inds,pos},
	seli = ToLowerCase[Last[Flatten[{#1/.#3}]]]===ToLowerCase[#2]&;
	self = (list=#1;key=#2[[1]];val=#2[[2]];Select[list,seli[key,val,json[[#]]]&])&;
	inds = Range[Length[json]];
	pos = Fold[self,inds,selection];
	If[sort==="", pos, pos[[Ordering[sort/.json[[pos]]]]]]
]


(* ::Subsubsection:: *)
(*MergeJSON*)


MergeJSON[json:{_?AssociationQ..}]:=Block[{keys},
keys=DeleteDuplicates[Flatten[Keys/@json]];
Association[If[#[[2]]==={},Nothing,#]&/@Thread[
keys->(If[Length[#]===1,First@#,#]&/@(
(DeleteDuplicates/@Transpose[(#/@keys)&/@json])/.Missing[___]->Nothing))]
]]



(* ::Subsubsection:: *)
(*AddToJson*)


AddToJson[json_,add_]:=MergeJSON[{json,
	Switch[add,
		"QMRITools",<|"ConversionSoftware"->"QMRITools.com","ConversionSoftwareVersion"->QMRITools`$InstalledVersion|>,
		_,{json,add}]}
]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
