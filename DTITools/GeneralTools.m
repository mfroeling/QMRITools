(* ::Package:: *)

(* ::Title:: *)
(*DTITools GeneralTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`GeneralTools`"];
$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];

Unprotect @@ Names["DTITools`GeneralTools`*"];
ClearAll @@ Names["DTITools`GeneralTools`*"];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


FileSelect::usage = 
"FileSelect[action] creates a systemdialog wicht returs file/foldername action can be \"FileOpen\", \"FileSave\" or \"Directory\".
FileSelect[action, {type}] same but allows the definition of filetypes for \"FileOpen\" and \"FileSave\" e.g. \"jpg\" or \"pdf\"."

DTItoolFunctions::usage = 
"DTItoolFunctions[] give list of all the DTItool functions."

DTItoolPackages::usage =
"DTItoolPackages[] give list of all the DTItool pacakges."

DTItoolFuncPrint::usage = 
"DTItoolFuncPrint[] gives a list of all the DTItool functions with their usage infomation."

MemoryUsage::usage = 
"MemoryUsage[] gives a table of which definitions use up memory.
MemoryUsage[n] gives a table of which definitions use up memory, where n is the amout of definitions to show."

ClearTemporaryVariables::usage = 
"ClearTemporaryVariables[] Clear temporary variables."

NumberTableForm::usage = 
"NumberTableForm[data] makes a right aligned table of the numbers with 3 decimal percision.
NumberTableForm[data, n] makes a right aligned table of the numbers with n decimal percision.";

CompilebleFunctions::usage = 
"CompilebleFunctions[] generates a list of all compilable functions."

Transpose2C::usage = "Transpose2C[data,trans] transposes a 2D array"
Transpose3C::usage = "Transpose3C[data,trans] transposes a 3D array"
Transpose4C::usage = "Transpose4C[data,trans] transposes a 4D array"
Transpose5C::usage = "Transpose5C[data,trans] transposes a 5D array"


(* ::Subsection:: *)
(*General Options*)


TableMethod::usage = 
"TableMethod is an option for NumberTableForm. It specifies which number form to uses. Values can be NumberForm, ScientificForm or EngineeringForm"


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(* File Select*)


Options[FileSelect]={WindowTitle -> Automatic};

SyntaxInformation[FileSelect] = {"ArgumentsPattern" -> {_, _.,_., OptionsPattern[]}};

FileSelect[action_, opts:OptionsPattern[]] := FileSelect[action, {""}, "*" ,opts]

FileSelect[action_, type:{_String ..}, opts:OptionsPattern[]] := FileSelect[action, type, "*", opts]

FileSelect[action_String, type : {_String ..}, name_String, opts:OptionsPattern[]] := Module[{input},
  If[!Element[action, {"FileOpen", "FileSave", "Directory"}], Return[]];
  input = If[(action == "FileOpen" || action == "FileSave"),
  	  	SystemDialogInput[action, {Directory[], {name ->{"*"<>type}}},opts],
  	  	SystemDialogInput["Directory", Directory[],opts]
    ];
  If[input === $Canceled, 
  	Print["Canceled!"], 
  	If[action == "Directory",
  		StringDrop[input,-1],
  		input
  		]
  	]
]

(* ::Subsection::Closed:: *)
(*DTItoolPackages*)


DTItoolPackages[] := 
 Last[StringSplit[#, "`"]] & /@ 
  DeleteDuplicates[
   StringReplace[#, "Private`" -> ""] & /@ Contexts["DTITools`*`"]]


(* ::Subsection::Closed:: *)
(*DTItoolFunctions*)


SyntaxInformation[DTItoolFunctions] = {"ArgumentsPattern" -> {_.,_.}};

DTItoolFunctions[]:=DTItoolFunctions[""];

DTItoolFunctions[toolb_String]:= Module[{packages},
	
  If[toolb==="",
  packages = 
   DeleteDuplicates[
    StringReplace[#, "Private`" -> ""] & /@ Contexts["DTITools`*`"]];
    ToExpression[Sort[Flatten[Names[# <> "*"] & /@ packages]]]
    ,
    Names["DTITools`" <> toolb <> "`*"]
  ]
  ]

DTItoolFunctions[p_Integer]:=Partition[DTItoolFunctions[], p, p, 1, ""] // Transpose // TableForm

DTItoolFunctions[toolb_String,p_Integer]:=Partition[DTItoolFunctions[toolb], p, p, 1, ""] // Transpose // TableForm


(* ::Subsection::Closed:: *)
(*DTItoolFuncPrint*)


SyntaxInformation[DTItoolFuncPrint] = {"ArgumentsPattern" -> {}};

DTItoolFuncPrint[]:=Module[{functions,packs},
	packs = Flatten[StringCases[$ContextPath, "DTITools`" ~~ x__ -> x]];
	functions = "DTITools`" <> # -> Names["DTITools`" <> # <> "*"] & /@ packs;
 (
     Print[Style[#[[1]], Bold, Large, Black]];
     (
        Print[Style[#, Bold, Medium, Black]];
        Information[#]
        ) & /@ #[[2]]
     ) & /@ functions;
]


(* ::Subsection::Closed:: *)
(*MemoryUsage*)


SyntaxInformation[MemoryUsage] = {"ArgumentsPattern" -> {_.}}

MemoryUsage[n__:100]:=With[{
  listing = myByteCount /@ Names[]
  },
 Labeled[
  Grid[
   Reverse@Take[Sort[listing], -n], Frame -> True, Alignment -> Left
   ],
  Column[
   {
    Style[
     "ByteCount for symbols without attributes Protected and \
ReadProtected in all contexts", 16, FontFamily -> "Times"]
    ,
    Style[
     Row@{"Total: ", Total[listing[[All, 1]]], " bytes for ", 
       Length[listing], " symbols"}, Bold]
    }, Center, 1.5
   ]
  , Top]
 ]
 
 
 myByteCount[symbolName_String] := Replace[
   ToExpression[symbolName, InputForm, Hold],
   Hold[x__] :> If[MemberQ[Attributes[x], Protected | ReadProtected],
     Sequence @@ {},
     {ByteCount[
       Through[{OwnValues, DownValues, UpValues, SubValues, 
          DefaultValues, FormatValues, NValues}[Unevaluated@x, 
         Sort -> False]]
       ],
      symbolName}
     ]
   ];


(* ::Subsection::Closed:: *)
(*MemoryUsage*)


SyntaxInformation[ClearTemporaryVariables] = {"ArgumentsPattern" -> {_.}}

ClearTemporaryVariables[] := Block[{names, attr},
  names = Names["DTITools`*`Private`*"];
  attr = Attributes /@ names;
  MapThread[If[#1 === {Temporary}, Clear[#2]] &, {attr, names}];
  ]


(* ::Subsection::Closed:: *)
(*NumberTableForm*)


Options[NumberTableForm] = Join[{TableMethod -> NumberForm}, Options[TableForm]];

SyntaxInformation[
   NumberTableForm] = {"ArgumentsPattern" -> {_, _., 
     OptionsPattern[]}};

NumberTableForm[dat_, opts : OptionsPattern[]] := 
  NumberTableForm[dat, 3, opts];

NumberTableForm[dat_, depth_, opts : OptionsPattern[]] := 
  Block[{opt, met},
   met = OptionValueBox[Method];
   met = If[
     MemberQ[{NumberForm, ScientificForm, EngineeringForm}, met], met,
      NumberForm];
   TableForm[Map[met[#, {20*depth, depth}] &, dat, {ArrayDepth[dat]}],
     FilterRules[{opts}, Options[TableForm]], TableAlignments -> Right]
   ];


(* ::Subsection::Closed:: *)
(*Compilble functions*)

SyntaxInformation[CompilebleFunctions] = {"ArgumentsPattern" -> {}};

CompilebleFunctions[]:=(Partition[Compile`CompilerFunctions[] // Sort, 50, 50, 1, 1] // Transpose) /. 1 -> {} // TableForm


(* ::Subsection::Closed:: *)
(*TransposeC*)


Transpose2C = Compile[{{dat, _Real, 2}, {trans, _Integer, 1}}, Transpose[dat, trans]];
Transpose3C = Compile[{{dat, _Real, 3}, {trans, _Integer, 1}}, Transpose[dat, trans]];
Transpose4C = Compile[{{dat, _Real, 4}, {trans, _Integer, 1}}, Transpose[dat, trans]];
Transpose5C = Compile[{{dat, _Real, 5}, {trans, _Integer, 1}}, Transpose[dat, trans]];


(* ::Section:: *)
(*End Package*)


End[]

SetAttributes[#,{Protected, ReadProtected}]&/@Names["DTITools`GeneralTools`*"];

EndPackage[]
