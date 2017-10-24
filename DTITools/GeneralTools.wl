(* ::Package:: *)

(* ::Title:: *)
(*DTITools GeneralTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`GeneralTools`", {"Developer`"}];
$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];

Unprotect @@ Names["DTITools`GeneralTools`*"];
ClearAll @@ Names["DTITools`GeneralTools`*"];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
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

MeanNoZero::usage = 
"MeanNoZero[data] calculates the mean of the data ignoring the zeros."

MeanStd::usage = 
"MeanStd[data] calculates the mean and standard deviation and reports it as a string."

MeanRange::usage = 
"MeanRange[Range] calculates the medain (50%) and standard deviation (14% and 86%) range and reports it as a string."

PadToDimensions::usage = 
"PadToDimensions[data, dim] pads the data to dimensions dim." 

SumOfSquares::usage = 
"SumOfSquares[{data1, data2, .... datan}] calculates the sum of squares of the datasets.
Output is the SoS and the weights, or just the SoS."

DevideNoZero::usage = 
"DevideNoZero[a, b] devides a/b but when b=0 the result is 0. a can be a number or vector."

LogNoZero::usage = 
"LogNoZero[val] return the log of the val which can be anny dimonsion array. if val=0 the output is 0."


(* ::Subsection::Closed:: *)
(*General Options*)


TableMethod::usage = 
"TableMethod is an option for NumberTableForm. It specifies which number form to uses. Values can be NumberForm, ScientificForm or EngineeringForm"

PadValue::usage = 
"PadValue is an option for PadToDimensions. It specifies the value of the padding."

OutputWeights::usage = 
"OutputWeights is an option for SumOfSqares. If True it also output the SoS weights."


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection:: *)
(*NoZeroFunctions*)


(* ::Subsubsection::Closed:: *)
(*DevideNoZero*)


SyntaxInformation[DevideNoZero] = {"ArgumentsPattern" -> {_,_}};

DevideNoZero[sig_,tot_]:=DevideNoZeroi[sig,tot]

DevideNoZeroi = Compile[{{sig, _Real, 1}, {tot, _Real, 0}}, If[tot == 0., sig 0., sig/tot], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> True];
DevideNoZeroi = Compile[{{sig, _Real, 0}, {tot, _Real, 0}}, If[tot == 0., 0., sig/tot], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> True];


(* ::Subsubsection::Closed:: *)
(*LogNoZero*)


SyntaxInformation[LogNoZero] = {"ArgumentsPattern" -> {_}};

LogNoZero[val_] := LogNoZeroi[val]

LogNoZeroi = Compile[{{val, _Real, 0}},If[val == 0., 0., Log[val]],RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> True]


(* ::Subsubsection::Closed:: *)
(*MeanNoZero*)


SyntaxInformation[MeanNoZero] = {"ArgumentsPattern" -> {_, _.}};

MeanNoZero[datai_] := Block[{data},
  data = N@Chop@TransData[datai, "l"];
  N@Chop@Map[Mean[DeleteCases[#, 0.] /. {} -> {0.}] &, data, {ArrayDepth[data] - 1}]
  ]


(* ::Subsection::Closed:: *)
(*SumOfSquares*)


Options[SumOfSquares] = {OutputWeights -> True}

SyntaxInformation[SumOfSquares] = {"ArgumentsPattern" -> {_,OptionsPattern[]}};

SumOfSquares[data_, OptionsPattern[]] := Block[{sos, weights, dataf},
  dataf = TransData[data, "l"];
  sos = SumOfSquaresi[dataf];
  
  If[OptionValue[OutputWeights],
   weights = DevideNoZero[dataf, sos];
   weights = Sqrt[DevideNoZero[weights, Mean[TransData[weights, "r"]]]];
   {sos, TransData[weights,"r"]},
   sos
   ]
  ]

SumOfSquaresi = Compile[{{sig, _Real, 1}}, Sqrt[Total[sig^2]], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> True];


(* ::Subsection::Closed:: *)
(*PadToDimensions*)


Options[PadToDimensions]={PadValue->0.}

SyntaxInformation[PadToDimensions] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

PadToDimensions[data_, dim_, OptionsPattern[]] := Block[{diffDim, padval, pad},
  padval = OptionValue[PadValue];
  diffDim = dim - Dimensions[data];
  pad = Transpose@{Floor[diffDim/2], Ceiling[diffDim/2]};
  ArrayPad[data, pad, padval]
  ]


(* ::Subsection::Closed:: *)
(*File Select*)


Options[FileSelect]={WindowTitle -> Automatic};

SyntaxInformation[FileSelect] = {"ArgumentsPattern" -> {_, _.,_., OptionsPattern[]}};

FileSelect[action_, opts:OptionsPattern[]] := FileSelect[action, {""}, "*" ,opts]

FileSelect[action_, type:{_String ..}, opts:OptionsPattern[]] := FileSelect[action, type, "*", opts]

FileSelect[action_String, type : {_String ..}, name_String, opts:OptionsPattern[]] := Module[{input},
  If[!Element[action, {"FileOpen", "FileSave", "Directory"}], Return[]];
  input = If[(action == "FileOpen" || action == "FileSave"),
  	  	SystemDialogInput[action, {Directory[], {name ->type}},opts],
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


(* ::Subsection:: *)
(*Package Functions*)


(* ::Subsubsection::Closed:: *)
(*DTItoolPackages*)


SyntaxInformation[DTItoolPackages] = {"ArgumentsPattern" -> {}};

DTItoolPackages[] := 
 TableForm[Last[StringSplit[#, "`"]] & /@ 
  DeleteDuplicates[
   StringReplace[#, "Private`" -> ""] & /@ Contexts["DTITools`*`"]]]


(* ::Subsubsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
(*DTItoolFuncPrint*)


SyntaxInformation[DTItoolFuncPrint] = {"ArgumentsPattern" -> {_.}};

DTItoolFuncPrint[]:=DTItoolFuncPrint[""]

DTItoolFuncPrint[toolb_String]:=Module[{functions,packs},
	packs = If[toolb==="",
		Sort[Flatten[StringCases[$ContextPath, "DTITools`" ~~ x__ -> x]]],
		{toolb}		
	];
	functions = "DTITools`" <> # -> Names["DTITools`" <> # <> "*"] & /@ packs;
 (
     Print[Style[#[[1]], Bold, Large, Black]];
     (
        Print[Style[#, Bold, Medium, Black]];
        Information[#]
        ) & /@ #[[2]]
     ) & /@ functions;
]


(* ::Subsubsection::Closed:: *)
(*Compilable functions*)


SyntaxInformation[CompilebleFunctions] = {"ArgumentsPattern" -> {}};

CompilebleFunctions[]:=(Partition[Compile`CompilerFunctions[] // Sort, 50, 50, 1, 1] // Transpose) /. 1 -> {} // TableForm


(* ::Subsection:: *)
(*Memory functions*)


(* ::Subsubsection::Closed:: *)
(*MemoryUsage*)


SyntaxInformation[MemoryUsage] = {"ArgumentsPattern" -> {_.}}

MemoryUsage[n__:100]:=With[{
  listing = myByteCount /@ Names[]
  },
 Labeled[
  Grid[
   Reverse@Take[Sort[listing], -n], Frame -> True, Alignment -> Center
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
     {Round[ByteCount[
       Through[{OwnValues, DownValues, UpValues, SubValues, 
          DefaultValues, FormatValues, NValues}[Unevaluated@x, 
         Sort -> False]]
       ]/1000000.,.01],
      symbolName}
     ]
   ];


(* ::Subsubsection::Closed:: *)
(*MemoryUsage*)


SyntaxInformation[ClearTemporaryVariables] = {"ArgumentsPattern" -> {_.}}

ClearTemporaryVariables[] := Block[{names, attr},
  names = Names["DTITools`*`Private`*"];
  attr = Attributes /@ names;
  MapThread[If[#1 === {Temporary}, Clear[#2]] &, {attr, names}];
  ]


(* ::Subsection:: *)
(*Number Functions*)


(* ::Subsubsection::Closed:: *)
(*NumberTableForm*)


Options[NumberTableForm] = Join[{TableMethod -> NumberForm}, Options[TableForm]];

SyntaxInformation[NumberTableForm] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

NumberTableForm[dat_, opts : OptionsPattern[]] := NumberTableForm[dat, 3, opts];

NumberTableForm[dat_, depth_, opts : OptionsPattern[]] := 
  Block[{opt, met},
   met = OptionValueBox[Method];
   met = If[
     MemberQ[{NumberForm, ScientificForm, EngineeringForm}, met], met,
      NumberForm];
   TableForm[Map[met[#, {20*depth, depth}] &, dat, {ArrayDepth[dat]}],
     FilterRules[{opts}, Options[TableForm]], TableAlignments -> Right]
   ];


(* ::Subsubsection::Closed:: *)
(*MeanStd*)


MeanStd[inp_] := Block[{dat}, 
	dat = inp /. {Mean[{}] -> Nothing, 0. -> Nothing};
	Quiet@Row[{NumberForm[Round[Mean[dat], .001], {3, 2}], NumberForm[Round[StandardDeviation[dat], .001], {3, 2}]}, "\[PlusMinus]"]
  ]


(* ::Subsubsection::Closed:: *)
(*MeanRange*)


MeanRange[inp_] := Block[{q1, q2, q3},
  {q1, q2, q3} = Quantile[inp /. {Mean[{}] -> Nothing, 0. -> Nothing}, {.14, .5, .86}];
  Quiet@Row[{NumberForm[Round[q2, .001], {3, 2}], "  (", NumberForm[Round[q1, .001], {3, 2}], " - ", NumberForm[Round[q3, .001], {3, 2}], ")"}]
  ]
  
MeanRange[inp_,quant_] := Block[{q1, q2, q3},
  {q1, q2, q3} = Quantile[inp /. {Mean[{}] -> Nothing, 0. -> Nothing}, {quant[[1]],.5,quant[[2]]}];
  Quiet@Row[{NumberForm[Round[q2, .001], {3, 2}], "  (", NumberForm[Round[q1, .001], {3, 2}], " - ", NumberForm[Round[q3, .001], {3, 2}], ")"}]
  ]


(* ::Section:: *)
(*End Package*)


End[]

If[DTITools`verbose,Print[Names["DTITools`GeneralTools`*"]]];
SetAttributes[#,{Protected, ReadProtected}]&/@Names["DTITools`GeneralTools`*"];

EndPackage[]
