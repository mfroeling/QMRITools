(* ::Package:: *)

(* ::Title:: *)
(*DTITools NiftiTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`NiftiTools`"];
$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];

Unprotect @@ Names["DTITools`NiftiTools`*"];
ClearAll @@ Names["DTITools`NiftiTools`*"];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


DcmToNii::usage = 
"DcmToNii[\"action\"] converts dicom to nii. \"action\" can be \"folder\" or \"file\", which converts all files in a folder to nii or just one file."

ImportNii::usage = 
"ImportNii[] promts to select the nii file to import.
ImportNii[\"file\"] imports the nii file. 
The option Method ->\"type\" can be \"data\" which exports {data, vox}, \"header\" which exports {hdr} or \"all\" which exports {data, vox, hdr}."

ImportNiiDiff::usage = 
"ImportNiiDiff[] will promt for the *.nii, *.bvec and *.bval file to import.
ImportNiiDiff[*.nii,*.bvec,*.bval] will import the given files.
The output will be {data,grad,bvec,vox}."

ImportBvalvec::usage =
"ImportBvalvec[] will promt to select the *.bval and *.bvec files.
ImportBvalvec[\"file.bvec\",\"file.bval\"] imports the given *.bval and *.bvec files." 

ImportBmat::usage =
"ImportBmat[] will promt to select the *.txt file containing the bmatrix.
ImportBmat[\"file.txt\"] imports the given *.txt file containing the bmatrix." 

ExportNii::usage = 
"ExportNii[data, vox] exports the nii file and will promt for a file name.
ExportNii[data, vox, \"file\"] exports the nii file to the location \"file\"."

ExportBval::usage = 
"ExportBval[bvals] exports the diffusion bvalues to exploreDTI format.
ExportBval[bvals, \"file\"] exports the diffusion bvalues to \"file\" in the exploreDTI format."

ExportBvec::usage = 
"ExportBvec[grad] exports the diffusion gradients to exploreDTI format.
ExportBvec[grad, \"file\"] exports the diffusion gradients to \"file\" in the exploreDTI format."

ExportBmat::usage = 
"ExportBmat[bmat] exports the diffusion bmatrix to exploreDTI format.
ExportBmat[bmat, \"file\"] exports the diffusion bmatrix to \"file\" in the exploreDTI format."

OpenMRIcron::usage = 
"OpenMRIcron[] promts to select the nii file to open in MRIcron.
ImOpenMRIcron[\"file\"] opens the nii file in MRIcron."


(* ::Subsection::Closed:: *)
(*Options*)


NumberType::usage = "NumberType of Nii file can be \"Integer\" of \"Real\"."

ImportResult::usage = "ImportResult is an option for OpenMRIcron and can be True or False";

NumberOfResults::usage = "NumberOfResults is an option for OpenMRIcron and should be an integer";

FlipBvec::usage = "FlipBvec is an option for ImportBvalvec"

RotateGradients::usage="RotateGradients is an option for ImportNiiDiff"


(* ::Subsection::Closed:: *)
(*Error Messages*)


DcmToNii::notfount = "dcm2nii.exe not found in $UserBaseDirectory or $BaseDirectory please install DTItools in correct directory.";

DcmToNii::type = "Input should be \"file\" or \"folder\".";

ImportNii::wht = "should be \"data\", \"header\" or \"all\".";

ExportNii::type = "NumberType should be \"Integer\" of \"Real\".";

OpenMRIcron::notfount = 
"mricron.exe not found in $UserBaseDirectory or $BaseDirectory please install DTItools in correct directory.";

OpenMRIcron::fil = "the file `1` does not exist.";

ImportNii::notfount = "the file `1` does not exist.";

(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection::Closed:: *)
(* DcmToNii*)


dcm2nii = $UserBaseDirectory <>"\\Applications\\DTITools\\Applications\\dcm2nii.exe";
dcm2nii = If[!FileExistsQ[dcm2nii],$BaseDirectory <>"\\Applications\\DTITools\\Applications\\dcm2nii.exe",dcm2nii];
dcm2nii = If[!FileExistsQ[dcm2nii], Message[DcmToNii::notfount],dcm2nii];

SyntaxInformation[DcmToNii] = {"ArgumentsPattern" -> {_,_.}};

DcmToNii[action_] := DcmToNii[action,""]; 

DcmToNii[action_,fstr_] := Module[{act,filfolin,folout,add,title},
	
	Print["Using Chris Rorden's dcm2nii.exe
http://www.mccauslandcenter.sc.edu/mricro/mricron/dcm2nii.html"];
		
	{act,title}=Switch[action,
		"folder",{"Directory","Select direcotry containig the dcm files"},
		"file",{"FileOpen","Select the dcm file to convert"},
		_,Return[Message[DcmToNii::type]]
		];
	If[
		Length[fstr]==2,
		filfolin=If[fstr[[1]]=="",FileSelect[act,WindowTitle->title],fstr[[1]]];
		If[filfolin==Null||folout==Null,Return[]];
		folout=If[fstr[[2]]=="",FileSelect[act,WindowTitle->"Select directory to put nii files in"],fstr[[2]]];
		,
		filfolin=If[fstr=="",FileSelect[act,WindowTitle->title],fstr];
		If[filfolin==Null||folout==Null,Return[]];
		folout=FileSelect["Directory",WindowTitle->"Select directory to put nii files in"];
		];
	If[filfolin==Null||folout==Null,Return[]];
	
	add=If[action == "file", "-v N "," "];
	
	Monitor[
		Quiet[Get[("!"<>dcm2nii<>" -a Y -b dcm2nii.ini -c Y -d N -e Y -f Y -g N -i Y -n Y -p Y "<>add<>"-x Y -o \""<>folout<>"\" \""<>filfolin<>"\"")]];
		,
		ProgressIndicator[Dynamic[Clock[Infinity]], Indeterminate]
	];
]


(* ::Subsection::Closed:: *)
(*ImportNii*)


Options[ImportNii]={Method->"data"};

SyntaxInformation[ImportNii] = {"ArgumentsPattern" -> {_.,OptionsPattern[]}};

ImportNii[opts:OptionsPattern[]]:=ImportNii["",opts];

ImportNii[fil_String:"",OptionsPattern[]] := Module[{strm, hdr, file, precision, adim, ddim, dim, dataf, data, vox,what,rotmat},
	
	what=OptionValue[Method];
	
	If[!MemberQ[{"data","header","all"},what],Return[Message[ImportNii::wht]]];
	
	file=If[fil=="",FileSelect["FileOpen",{"*.nii"},WindowTitle->"Select the nii file to import"],fil];	
	If[file == Null || file === $Canceled, Return[]];
	
	If[!FileExistsQ[file],Return[Message[ImportNii::notfount,file]]];
	
	strm = OpenRead[file, BinaryFormat -> True];
	hdr = "header" -> {
     "headerKey" -> {
       "size" -> BinaryReadList[strm, "Integer32", 1],
       "dataType" -> FromCharacterCode[DeleteCases[BinaryReadList[strm, "UnsignedInteger8", 10],0]](*directchar*),
       "dbName" ->FromCharacterCode[DeleteCases[BinaryReadList[strm, "UnsignedInteger8", 18],0]](*directchar*),
       "extents" -> BinaryReadList[strm, "Integer32", 1],
       "sessionError" -> BinaryReadList[strm, "Integer16", 1],
       "regular" ->FromCharacterCode[DeleteCases[BinaryReadList[strm, "UnsignedInteger8", 1],0]](*directchar*),
       "dimInfo" -> BinaryReadList[strm, "UnsignedInteger8", 1]
       }
     ,
     "imageDimensions" -> {
       "dim" -> BinaryReadList[strm, "Integer16", 8],
       "intentP1" -> BinaryReadList[strm, "Real32", 1],
       "intentP2" -> BinaryReadList[strm, "Real32", 1],
       "intentP3" -> BinaryReadList[strm, "Real32", 1],
       "intentCode" -> BinaryReadList[strm, "Integer16", 1],
       "dataType" -> BinaryReadList[strm, "Integer16", 1],
       "bitPix" -> BinaryReadList[strm, "Integer16", 1],
       "sliceStart" -> BinaryReadList[strm, "Integer16", 1],
       "pixDim" -> BinaryReadList[strm, "Real32", 8],
       "voxOffset" -> BinaryReadList[strm, "Real32", 1],
       "scaleSlope" -> BinaryReadList[strm, "Real32", 1],
       "scaleInteger" -> BinaryReadList[strm, "Real32", 1],
       "sliceEnd" -> BinaryReadList[strm, "Integer16", 1],
       "sliceCode" -> BinaryReadList[strm, "UnsignedInteger8", 1],
       "xyztUnits" -> BinaryReadList[strm, "UnsignedInteger8", 1],
       "calMax" -> BinaryReadList[strm, "Real32", 1],
       "calMin" -> BinaryReadList[strm, "Real32", 1],
       "sliecDuration" -> BinaryReadList[strm, "Real32", 1],
       "tOffset" -> BinaryReadList[strm, "Real32", 1],
       "glMax" -> BinaryReadList[strm, "Integer32", 1],
       "glMin" -> BinaryReadList[strm, "Integer32", 1]
       }
     ,
     "dataHistory" -> {
       "descrip" ->FromCharacterCode[DeleteCases[BinaryReadList[strm, "UnsignedInteger8", 80],0]](*directchar*),
       "auxFile" ->FromCharacterCode[DeleteCases[BinaryReadList[strm, "UnsignedInteger8", 24],0]](*directchar*),
       "qformCode" -> BinaryReadList[strm, "Integer16", 1],
       "sformCode" -> BinaryReadList[strm, "Integer16", 1],
       "quaternB" -> BinaryReadList[strm, "Real32", 1],
       "quaternC" -> BinaryReadList[strm, "Real32", 1],
       "quaternD" -> BinaryReadList[strm, "Real32", 1],
       "qOffsetx" -> BinaryReadList[strm, "Real32", 1],
       "qOffsety" -> BinaryReadList[strm, "Real32", 1],
       "qOffsetz" -> BinaryReadList[strm, "Real32", 1],
       "sRowx" -> BinaryReadList[strm, "Real32", 4],
       "sRowy" -> BinaryReadList[strm, "Real32", 4],
       "sRowz" -> BinaryReadList[strm, "Real32", 4],
       "intentName" -> FromCharacterCode[DeleteCases[BinaryReadList[strm, "UnsignedInteger8", 16],0]](*directchar*),
       "magic" ->FromCharacterCode[DeleteCases[BinaryReadList[strm, "UnsignedInteger8", 4],0]](*directchar*)
       }
     };
  
  precision = Switch [
    ("dataType" /. ("imageDimensions" /. ("header" /. hdr)))[[1]],
    1, "Bit",
    2, "Integer8",
    4, "Integer16",
    8, "Integer32",
    16, "Real32",
    32, "Real32",
    64, "Real64",
    128, "Unsignedinteger8",
    256, "Integer8",
    511, "Real32",
    512, "Integer16"(*"UnsignedInteger16"*),
    768, "Unsignedinteger32",
    1024, "Integer64",
    1280, "Unsignedinteger64",
    1792, "Real64",
    _, Message["This datatype is not supported"];
    ];
  
  adim = ("dim" /. ("imageDimensions" /. ("header" /. hdr)));
  ddim = adim[[1]];
  dim = adim[[2 ;; ddim + 1]];
  BinaryReadList[strm, "Real32", 1];(*???*)
  dataf = BinaryReadList[strm, precision];
  Close[strm];
  
  (*(("scaleSlope" /. ("imageDimensions" /. ("header" /. hdr)))[[1]]) * *)
  
  data = Fold[Partition, dataf, Drop[dim, {-1}]];
  If[ddim == 4, data = Transpose[data,{2,1,3,4}]];
  
  data = Map[Reverse[#] &, data, {ddim - 2}];
  
    
  If[Positive[("sRowx" /. ("dataHistory" /. ("header"/. hdr)))[[1]]],
  	data = Map[Reverse[#] &, data, {ddim - 1}]
  	];
  (*data = Map[Reverse[#] &, data, {ddim - 1}];*)
  (*If[Negative[("sRowy" /. ("dataHistory" /. ("header"/. hdr)))[[2]]],data = Map[Reverse[#] &, data, {ddim - 3}];];  
  If[Negative[("sRowz" /. ("dataHistory" /. ("header"/. hdr)))[[3]]],data = Reverse[data];];*)
  

  vox = If[ddim==2,
    	Reverse[("pixDim" /. ("imageDimensions" /. ("header" /. hdr)))[[2 ;; 3]]]
    	,
    	Reverse[("pixDim" /. ("imageDimensions" /. ("header" /. hdr)))[[2 ;; 4]]]
    ];
  
  Switch[what,
  	"data",{data,vox},
  	"header",hdr,
  	"all",
  		  rotmat = ({1, 1, -1} {"sRowx", "sRowy", "sRowz"} /. ("dataHistory" /. ("header" /. hdr)))[[All,1 ;; 3]]/ConstantArray[Reverse[vox], 3];
  		  rotmat = DiagonalMatrix[{1, -1, 1}].ConstantArray[Diagonal[Sign[rotmat]], 3] rotmat;
  		  {data, vox, hdr, rotmat}
  ]
]


(* ::Subsection::Closed:: *)
(*ImportBvalvec*)

Options[ImportBvalvec]={FlipBvec->False};

SyntaxInformation[ImportBvalvec] = {"ArgumentsPattern" -> {_.,_.,OptionsPattern[]}};

ImportBvalvec[file_?StringQ,opts:OptionsPattern[]] := Module[{valf, vecf},
  valf = file<>".bval";
  vecf = file<>".bvec";
  ImportBvalvec[valf,vecf,opts]
  ]

ImportBvalvec[valf_,vecf_,OptionsPattern[]] := Module[{grads},
	{
	Flatten[ LineToList[valf] ]//N,
	grads=Round[LineToList[vecf],0.0001];
	grads=If[OptionValue[FlipBvec],
		{1, -1, 1}#&/@grads,
		{1,-1,1}RotateLeft[#]&/@grads
		];
	If[Negative[#[[3]]],-#,#]&/@grads
	}
]

ImportBvalvec[___,opts:OptionsPattern[]] := Module[{valf, vecf},
  valf = FileSelect["FileOpen", {"*.bval"}, WindowTitle -> "Select *.bval"];
  vecf = FileSelect["FileOpen", {"*.bvec"}, WindowTitle -> "Select *.bvec"];
  ImportBvalvec[valf,vecf,opts]
  ]

LineToList[file_] := Module[{grads,tmp},
  grads = Import[file, "Lines"];
  grads = Transpose[(
  	tmp = StringReplace[#, {" " -> ",", "\t" -> ",", "E" -> "*10^"}];
	tmp = If[StringTake[tmp, -1] == ",", StringDrop[tmp, -1], tmp];
	tmp = ToExpression["{" <> tmp <> "}"]
  	) & /@ grads]
  ]


(* ::Subsection::Closed:: *)
(*ImportBmat*)


SyntaxInformation[ImportBmat] = {"ArgumentsPattern" -> {_.,_.}};

ImportBmat[] := ImportBmat[FileSelect["FileOpen", {"*.txt"}, WindowTitle -> "Select *.txt"]]

ImportBmat[fil_String] := Module[{bmati},
  If[fil == Null, Return[]];
  bmati = DeleteCases[#, ""] & /@ Import[fil, "Data"];
  Append[{-1, -1, -1, -1, -1, -1} #, 1] & /@ bmati[[All, {4, 1, 6, 2, 5, 3}]]
  ]


(* ::Subsection::Closed:: *)
(*ImportNiiDiff*)

Options[ImportNiiDiff]={RotateGradients->True,FlipBvec->True}

SyntaxInformation[ImportNiiDiff]= {"ArgumentsPattern" -> {_.,_.,_.,OptionsPattern[]}};

ImportNiiDiff[OptionsPattern[]]:=Module[{data,grad,bvec,vox,hdr,mat},
	{data,vox,hdr,mat}=ImportNii[Method -> "all"];
	{bvec, grad}=ImportBvalvec[FlipBvec->OptionValue[FlipBvec]];
	{data,Round[If[OptionValue[RotateGradients],grad.mat, grad],.0001],bvec,vox}
]

ImportNiiDiff[file_String,OptionsPattern[]]:=Module[{data,grad,bvec,vox,hdr,mat},
	{data,vox,hdr,mat}=ImportNii[file,Method -> "all"];
	{bvec, grad}=ImportBvalvec[StringDrop[file,-4],FlipBvec->OptionValue[FlipBvec]];
	{data,Round[If[OptionValue[RotateGradients],grad.mat, grad],.0001],bvec,vox}
]

ImportNiiDiff[fnii_String,fvec_String,fval_String,OptionsPattern[]]:=Module[{data,grad,bvec,vox,hdr,mat},
	{data,vox,hdr,mat}=ImportNii[fnii,Method -> "all"];
	{bvec, grad} = ImportBvalvec[fval, fvec,FlipBvec->OptionValue[FlipBvec]];
	{data,Round[If[OptionValue[RotateGradients],grad.mat, grad],.0001],bvec,vox}
]


(* ::Subsection::Closed:: *)
(*ExportNii*)

ArrangeData[data_]:=ArrangeData[data,ArrayDepth[data]]

ArrangeData[data_,depth_]:= Flatten[
	If[depth == 4,
		Transpose[Reverse[Reverse[data, depth], depth - 1]],
		Reverse[Reverse[data, depth], depth - 1]
		]
	]

SyntaxInformation[ExportNii] = {"ArgumentsPattern" -> {_,_,_., OptionsPattern[]}};

Options[ExportNii]={NumberType->"Integer"}

ExportNii[dato_, vox_,opts:OptionsPattern[]] := ExportNii[dato, vox, "",opts]

ExportNii[dato_, vox_, fil_,OptionsPattern[]] := Block[{datao, type, dim, depth, dimo, voxo, srowx, srowy, srowz, ftype, 
   strm, fileo, typeo},
  
  fileo = If[fil == "", FileSelect["FileSave",{"*.nii"},"nifti",WindowTitle->"Select the destination file"], fil];
  If[fileo == Null, Return[]];
  
  ftype=OptionValue[NumberType];
  
  dim = Dimensions[dato];
  depth = ArrayDepth[dato];
  
  Switch[
  	ftype,
  	"Integer", 
  	datao = Round[ArrangeData[dato,depth]];
  	type=If[Max[datao]>32768 || Min[datao]<-32768,"Integer32","Integer16"];
  	,
  	"Real",
  	datao = N[ArrangeData[dato,depth]]; 
  	type = "Real32";
];
    
  typeo = Switch [
    type,
    "Bit", {1, 1},
    "Unsigned8", {2, 8},
    "Integer16", {4, 16},
    "Integer32", {8, 32},
    "Real32", {16, 32},
    "Real32", {32, 32},
    "Real64", {64, 64},
    "Unsignedinteger8", {128, 8},
    "Integer8", {256, 8},
    "Real32", {511, 32},
    "Unsignedinteger16", {512, 16},
    "Unsignedinteger32", {768, 32},
    "Integer64", {1024, 64},
    "Unsignedinteger64", {1280, 64},
    "Real64", {1792, 64},
    _, Message["This datatype is not supported"]
    ];

(*
  Print["Reverse stuff"];Pause[5];
  (*datao = Map[Reverse[#] &, datao, {depth - 1}];
  datao = Map[Reverse[#] &, datao, {depth - 2}];*)
  datao = Reverse[Reverse[datao, depth], depth - 1];
    
  Print["Transpose stuff"];Pause[5];
  datao = If[depth == 4, Transpose[datao], datao];
  (*datao = Map[Reverse[#] &, datao, {depth - 2}];*)
  (*datao = Map[Reverse[#] &, datao, {depth - 3}];*)
   
  Print["flatten stuff"];Pause[5];
  datao = Flatten[datao];
*) 
 
  dimo = PadRight[Flatten[{depth, If[Length[dim] == 4, dim[[{4, 3, 1, 2}]], Reverse[dim]]}], 8, 1];
  voxo = PadRight[Flatten[{1., Reverse[vox]}], 8, 0];
  
  srowx = {vox[[3]], 0., 0., -N[vox[[3]] dim[[depth]]/2]};
  srowy = {0., vox[[2]], 0., -N[vox[[2]] dim[[depth - 1]]/2]};
  srowz = {0., 0., vox[[1]], -N[vox[[1]] dim[[1]]/2]};
  
  strm = OpenWrite[fileo, BinaryFormat -> True];
  (*size"*)BinaryWrite[strm, 348, "Integer32"];
  (*dataType*)BinaryWrite[strm, ConstantArray[0, {10}],"UnsignedInteger8"](*directchar*);
  (*dbName*)BinaryWrite[strm, ConstantArray[0, {18}],"UnsignedInteger8"](*directchar*);
  (*extents*)BinaryWrite[strm, 0, "Integer32"];
  (*sessionError*)BinaryWrite[strm, 0, "Integer16"];
  (*regular*)BinaryWrite[strm, ToCharacterCode["r"],"UnsignedInteger8"](*directchar*);
  (*dimInfo*)BinaryWrite[strm, 0, "UnsignedInteger8"];
  
  (*dim*)BinaryWrite[strm, dimo, "Integer16"];
  (*intentP1*)BinaryWrite[strm, 0., "Real32"];
  (*intentP2*)BinaryWrite[strm, 0., "Real32"];
  (*intentP3*)BinaryWrite[strm, 0., "Real32"];
  (*intentCode*)BinaryWrite[strm, 0, "Integer16"];
  (*dataType*)BinaryWrite[strm, typeo[[1]], "Integer16"];
  (*bitPix*)BinaryWrite[strm, typeo[[2]], "Integer16"];
  (*sliceStart*)BinaryWrite[strm, 0, "Integer16"];
  (*pixDim*)BinaryWrite[strm, voxo, "Real32"];
  (*voxOffset*)BinaryWrite[strm, 352., "Real32"];
  (*scaleSlope*)BinaryWrite[strm, 1., "Real32"];
  (*scaleInteger*)BinaryWrite[strm, 0., "Real32"];
  (*sliceEnd*)BinaryWrite[strm, 0, "Integer16"];
  (*sliceCode*)BinaryWrite[strm, 0, "UnsignedInteger8"];
  (*xyztUnits*)BinaryWrite[strm, 2, "UnsignedInteger8"];(*0-unknown 1-meters 2-millimeter 3-micrometers*)
  (*calMax*)BinaryWrite[strm, 0., "Real32"];
  (*calMin*)BinaryWrite[strm, 0., "Real32"];
  (*sliecDuration*)BinaryWrite[strm, 0., "Real32"];
  (*tOffset*)BinaryWrite[strm, 0., "Real32"];
  (*glMax*)BinaryWrite[strm, 0, "Integer32"];
  (*glMin*)BinaryWrite[strm, 0, "Integer32"];
  
  (*descrip*)BinaryWrite[strm,PadRight[ToCharacterCode["Created with DTItools"], 80, 0],"UnsignedInteger8"](*directchar*);
  (*auxFile*)BinaryWrite[strm, PadRight[ToCharacterCode["none"], 24, 0],"UnsignedInteger8"](*directchar*);
  (*qformCode*)BinaryWrite[strm, 2, "Integer16"];
  (*sformCode*)BinaryWrite[strm, 1, "Integer16"];
  (*quaternB*)BinaryWrite[strm, 0., "Real32"];
  (*quaternC*)BinaryWrite[strm, 0., "Real32"];
  (*quaternD*)BinaryWrite[strm, 0., "Real32"];
  (*qOffsetx*)BinaryWrite[strm, srowx[[4]], "Real32"];
  (*qOffsety*)BinaryWrite[strm, srowy[[4]], "Real32"];
  (*qOffsetz*)BinaryWrite[strm, srowz[[4]], "Real32"];
  (*sRowx*)BinaryWrite[strm, srowx, "Real32"];
  (*sRowy*)BinaryWrite[strm, srowy, "Real32"];
  (*sRowz*)BinaryWrite[strm, srowz, "Real32"];
  (*intentName*)BinaryWrite[strm, ConstantArray[0, 16],"UnsignedInteger8"](*directchar*);
  (*magic*)BinaryWrite[strm, {110, 43, 49, 0}, "UnsignedInteger8"](*directchar*);
  
  (*???*)BinaryWrite[strm, 0., "Real32"];
  
  (*data*)BinaryWrite[strm, datao, type];
  Close[strm];
  fileo
]

(* ::Subsection::Closed:: *)
(*ExportBval*)


SyntaxInformation[ExportBval] = {"ArgumentsPattern" -> {_,_.}};

ExportBval[bv_] := ExportBval[bv, ""]

ExportBval[bv_,fil_String] := Module[{file,bve},
  
  file = If[fil == "", FileSelect["FileSave", {"*.bval"}, "bval file", WindowTitle -> "Select the destination file"], fil];
  If[file === Null, Return[]];
  file = If[StringTake[file, -5] == ".bval", file, file <> ".bval"];
  
  bve=StringJoin[ToString[#] <> " " & /@ bv];
  Export[file, bve, "Text"]
  ]


(* ::Subsection::Closed:: *)
(*ExportBvec*)


SyntaxInformation[ExportBvec] = {"ArgumentsPattern" -> {_,_.}};

ExportBvec[grad_] := ExportBvec[grad, ""]

ExportBvec[grad_, fil_String] := Module[{file,grade},
  
  file = If[fil == "", FileSelect["FileSave", {"*.bvec"}, "bvec file", WindowTitle -> "Select the destination file"], fil];
  If[file === Null, Return[]];
  file = If[StringTake[file, -5] == ".bvec", file, file <> ".bvec"];
  
  grade=StringJoin[(ToString[#] <> " " & /@ #)] & /@ Transpose[Round[{1,1,-1}RotateRight[#]&/@grad, .0001]];
  Export[file, grade, "Text"]
  ]


(* ::Subsection::Closed:: *)
(*ExportBmat*)


SyntaxInformation[ExportBmat] = {"ArgumentsPattern" -> {_,_.}};

ExportBmat[bmat_] := ExportBmat[bmat, ""]

ExportBmat[bmat_, fil_String] := Module[{bmate, file},
  
  file = If[fil == "", FileSelect["FileSave", {"*.txt"}, "txt file", WindowTitle -> "Select the destination file"], fil];
  If[file == Null, Return[]];
  file = If[StringTake[file, -4] == ".txt", file, file <> ".txt"];
  bmate = If[Length[bmat[[1]]]==7,
  	StringReplace[ToString[-{1,1,1,1,1,1}#&/@Round[bmat[[All, {2,4,6,1,5,3}]], 0.0001]], {"{{" -> "","}}" -> "", "}, {" -> "\n", ", " -> "\t\t"}],
  	StringReplace[ToString[{1,1,1,1,1,1}#&/@Round[bmat[[All, {2,4,6,1,5,3}]], 0.0001]], {"{{" -> "","}}" -> "", "}, {" -> "\n", ", " -> "\t\t"}]
  ];
  
  Export[file, bmate]
  ]


(* ::Subsection::Closed:: *)
(*OpenMRIcron*)


FindMRIcron[] := Module[{mricron},
  mricron = $UserBaseDirectory <> 
    "\\Applications\\DTITools\\Applications\\mricron.exe";
  mricron = 
   If[! FileExistsQ[mricron], $BaseDirectory <> 
     "\\Applications\\DTITools\\Applications\\mricron.exe", mricron];
  If[! FileExistsQ[mricron], Return[Message[OpenMRIcron::notfount]], 
   mricron]
  ]

Clear[OpenMRIcron]
Options[OpenMRIcron] = {ImportResult -> True, NumberOfResults -> 1};
SyntaxInformation[
   OpenMRIcron] = {"ArgumentsPattern" -> {_., OptionsPattern[]}};
OpenMRIcron[opts : OptionsPattern[]] := OpenMRIcron["", opts];
OpenMRIcron[filei_, OptionsPattern[]] := 
 Block[{file, mricron, command, num},
  mricron = FindMRIcron[];
  If[! (mricron === Null),
   file = If[filei === "", FileSelect["FileOpen"], filei] /. 
      Null -> "Cancel";
   
   If[FileExistsQ[file],
    command = "!START /MAX " <> mricron <> " " <> file;
    Get[command];
    num = OptionValue[NumberOfResults];
    num = If[NumericQ[num], num, 1];
    If[OptionValue[ImportResult], 
     If[num == 1, ImportNii[], Transpose@Table[ImportNii[], {i, num}]]]
    ,
    Return[Message[OpenMRIcron::fil, file]];
    ]
   ]
  ]


(* ::Section:: *)
(*End Package*)


End[]

SetAttributes[#,{Protected, ReadProtected}]&/@Names["DTITools`NiftiTools`*"];

EndPackage[]
