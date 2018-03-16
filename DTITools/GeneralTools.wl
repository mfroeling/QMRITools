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

SetupDataStructure::usage = 
"SetupDataStructure[dcmFolder] makes nii folders and generates nii files for a directory of dmc data where the data is structured per subject."

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

MedianNoZero::usage = 
"MedianNoZero[data] calculates the Median of the data ignoring the zeros."

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

ExpNoZero::usage = 
"ExpNoZero[val] return the Exp of the val which can be anny dimonsion array. if val=0 the output is 0."

RMSNoZero::usage = 
"RMSNoZero[vec] return the RMS error of the vec which can be anny dimonsion array. if vec={0...} the output is 0. Zeros are ignored"

MADNoZero::usage = 
"MADNoZero[vec] return the MAD error of the vec which can be anny dimonsion array. if vec={0...} the output is 0. Zeros are ignored"


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
(*ExpNoZero*)


SyntaxInformation[ExpNoZero] = {"ArgumentsPattern" -> {_}};

ExpNoZero[val_] := ExpNoZeroi[val]

ExpNoZeroi = Compile[{{val, _Real, 0}},If[val == 0., 0., Exp[val]],RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> True]


(* ::Subsubsection::Closed:: *)
(*RMSNoZero*)


SyntaxInformation[RMSNoZero] = {"ArgumentsPattern" -> {_}};

RMSNoZero[vec_] := RMSNoZeroi[If[ArrayDepth[vec] > 1, TransData[vec, "l"], vec]]

RMSNoZeroi = Compile[{{vec, _Real, 1}}, If[AllTrue[vec, # === 0. &], 0., RootMeanSquare[DeleteCases[vec, 0.]]], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*MADNoZero*)


SyntaxInformation[MADNoZero] = {"ArgumentsPattern" -> {_}};

MADNoZero[vec_] := MADNoZeroi[If[ArrayDepth[vec] > 1, TransData[vec, "l"], vec]]

MADNoZeroi = Compile[{{vec, _Real, 1}}, If[AllTrue[vec, # === 0. &], 0., MedianDeviation[DeleteCases[vec, 0.]]], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*MeanNoZero*)


SyntaxInformation[MeanNoZero] = {"ArgumentsPattern" -> {_, _.}};

MeanNoZero[datai_] := Block[{data},
  data = N@Chop@TransData[datai, "l"];
  N@Chop@Map[Mean[DeleteCases[#, 0.] /. {} -> {0.}] &, data, {ArrayDepth[data] - 1}]
  ]


(* ::Subsubsection::Closed:: *)
(*MedianNoZero*)


SyntaxInformation[MedianNoZero] = {"ArgumentsPattern" -> {_, _.}};

MedianNoZero[datai_] := Block[{data},
  data = N@Chop@TransData[datai, "l"];
  N@Chop@Map[Median[DeleteCases[#, 0.] /. {} -> {0.}] &, data, {ArrayDepth[data] - 1}]
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
   (*weights = Sqrt[DevideNoZero[weights, Mean[TransData[weights, "r"]]]];*)
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
(*SetupDataStructure*)

SetupDataStructure[dcmFolder_] := 
 Module[{folderdcm, foldernii, folderout, folders,fol, niiFolder, outFolder},
  folderdcm = Directory[] <> "\\" <> # & /@ Select[FileNames["*", "dcm"], DirectoryQ];
  foldernii = StringReplace[#, "dcm" -> "nii"] & /@ folderdcm;
  folderout = StringReplace[#, "dcm" -> "out"] & /@ folderdcm;
  folders = Transpose[{folderdcm, foldernii, folderout}];
  
  fol = Last@FileNameSplit[dcmFolder];
  niiFolder = StringReplace[dcmFolder, fol -> "nii"];
  outFolder = StringReplace[dcmFolder, fol -> "out"];
  If[! DirectoryQ[niiFolder], CreateDirectory[niiFolder]];
  If[! DirectoryQ[outFolder], CreateDirectory[outFolder]];
  
  (*create nii files*)
  If[! DirectoryQ[#[[2]]], CreateDirectory[#[[2]]]; DcmToNii[#[[1 ;; 2]]]] & /@ folders;
  
  folders
]


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

MemoryUsage[size_: 1] := Module[{},
  NotebookClose[memwindow];
  memwindow = 
   CreateWindow[DialogNotebook[{CancelButton["Close", DialogReturn[]],
      Manipulate[
       If[listing === {},
        "Noting found",
        TableForm[Join[#, {
             Row[Dimensions[ToExpression["Global`" <> #[[1]]]], "x"],
             Head[ToExpression["Global`" <> #[[1]]]],
             ClearButton["Global`" <> #[[1]]]
             }] & /@ listing, 
         TableHeadings -> {None, {"Name", "Size (MB)", "Dimensions", "Head"}}]
        ]
       ,
       {{msize, size, "minimum size [MB]"}, {1, 10, 100, 1000}},
       Button["Update", listing = MakeListing[msize]],
       {listing, ControlType -> None},
       ContentSize -> {500, 600},
       Paneled -> False,
       AppearanceElements -> None,
       Initialization :> {
         MakeListing[mb_] := Reverse@SortBy[Select[myByteCount[Names["Global`*"]], #[[2]] > mb &], Last],
         ClearButton[name_] := Button["Clear",       
           Replace[ToExpression[name, InputForm, Hold], Hold[x__] :> Clear[Unevaluated[x]]];
           listing = MakeListing[msize]
         ],
         listing = MakeListing[size];
         }
       ]}, WindowTitle -> "Plot data window", Background -> White]
    ];
  ]


SetAttributes[myByteCount, Listable ]

myByteCount[symbolName_String] := Replace[
	ToExpression[symbolName, InputForm, Hold], Hold[x__] :> If[MemberQ[Attributes[x], Protected | ReadProtected],
     Sequence @@ {},
     (*output size in MB and name*)
     {StringDelete[symbolName, "Global`"], 
      Round[ByteCount[Through[{OwnValues, DownValues, UpValues, SubValues, DefaultValues, FormatValues, NValues}[Unevaluated@x, Sort -> False]]]/1000000., .01]}
     ]
   ];

(*
MemoryUsage2[n__: 100] := 
 With[{listing = myByteCount /@ Names[]}, 
  Labeled[Grid[Reverse@Take[Sort[listing], -n], Frame -> True, 
    Alignment -> Center], 
   Column[{Style[
      "ByteCount for symbols without attributes Protected and \
ReadProtected in all contexts", 16, FontFamily -> "Times"], 
     Style[Row@{"Total: ", Total[listing[[All, 1]]], " bytes for ", 
        Length[listing], " symbols"}, Bold]}, Center, 1.5], Top]]

myByteCount[symbolName_String] := 
  Replace[ToExpression[symbolName, InputForm, Hold], 
   Hold[x__] :> 
    If[MemberQ[Attributes[x], Protected | ReadProtected], 
     Sequence @@ {}, {Round[
       ByteCount[
         Through[{OwnValues, DownValues, UpValues, SubValues, 
            DefaultValues, FormatValues, NValues}[Unevaluated@x, 
           Sort -> False]]]/1000000., .01], symbolName}]];
*)


(* ::Subsubsection::Closed:: *)
(*MemoryUsage*)


SyntaxInformation[ClearTemporaryVariables] = {"ArgumentsPattern" -> {_.}}

ClearTemporaryVariables[] := Block[{names, attr},
  names = Names["DTITools`*`Private`*"];
  attr = Attributes /@ names;
  MapThread[If[#1 === {Temporary}, ClearAll[#2]] &, {attr, names}];
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

EndPackage[]
