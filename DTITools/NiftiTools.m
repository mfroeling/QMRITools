(* ::Package:: *)

(* ::Title:: *)
(*DTITools NiftiTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`NiftiTools`", {"Developer`"}];
$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];

Unprotect @@ Names["DTITools`NiftiTools`*"];
ClearAll @@ Names["DTITools`NiftiTools`*"];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


DcmToNii::usage =
"DcmToNii[] converts a dicom folder to nii. 
DcmToNii[\"action\"] converts dicom to nii. \"action\" can be \"folder\" or \"file\", which converts all files in a folder to nii or just one file.
DcmToNii[{\"input\",\"ouput\"}] converts the \"input\" dicom folder to nii files which are place in the \"output\" folder."

ImportNii::usage = 
"ImportNii[] promts to select the nii file to import.
ImportNii[\"file\"] imports the nii file. 

The default output is {data, vox}, however using NiiMethod various outputs can be given."

ImportNiiDiff::usage = 
"ImportNiiDiff[] will promt for the *.nii, *.bvec and *.bval file to import.
ImportNiiDiff[*.nii,*.bvec,*.bval] will import the given files.
The output will be {data,grad,bvec,vox}."

ImportBvalvec::usage =
"ImportBvalvec[] will promt to select the *.bval and *.bvec files.
ImportBvalvec[\"file.bvec\",\"file.bval\"] imports the given *.bval and *.bvec files." 

ImportBval::usage = 
"ImportBval[] will promt to select the *.bval file.
ImportBval[\"file.bval\"] imports the given *.bval file." 

ImportBvec::usage = 
"ImportBvec[] will promt to select the *.bvec file.
ImportBvec[\"file.bvec\"] imports the given *.bvec file." 

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

ExtractNiiFiles::usage =
"ExtractNiiFiles[] extracts all nii.gz files to .nii files in the selected folder.
ExtractNiiFiles[folder] extracts all nii.gz files to .nii files in folder."

CompressNiiFiles::usage =
"CompressNiiFiles[] compresses all nii files to .nii.gz files in the selected folder.
ECompressNiiFilesfolder] compresses all nii files to .nii.gz files in folder."


(* ::Subsection:: *)
(*Options*)


NiiMethod::usage = "NiiMethod is an option for ImportNIi. Values can be \"data\", \"dataTR\", \"header\", \"scaling\", \"headerMat\", \"rotation\", \"all\"."

NiiScaling::usage = "NiiScaling is an option for ImportNii. It scales the nii values with scale slope and offset for quantitative data."

CompressNii::usage = "CompressNii is an option for DcmToNii and ExportNii. If set True .nii.gz files will be created."

NumberType::usage = "NumberType is an option of Export Nii. The number type of Nii file can be \"Integer\", \"Real\", \"Complex\", or \"Automatic\"."

RotateGradients::usage = "RotateGradients is an option for ImportNiiDiff."

FlipBvec::usage = "FlipBvec is an option for ImportBvalvec."

ImportResult::usage = "ImportResult is an option for OpenMRIcron and can be True or False."

NumberOfResults::usage = "NumberOfResults is an option for OpenMRIcron and should be an integer."


(* ::Subsection:: *)
(*Error Messages*)


DcmToNii::notfount = "dcm2nii.exe not found in $UserBaseDirectory or $BaseDirectory please install DTItools in correct directory."

DcmToNii::type = "Input should be \"file\" or \"folder\"."

ImportNii::wht = "should be \"data\", \"dataTR\", \"header\", \"scaling\", \"headerMat\", \"rotation\", \"all\"."

ImportNii::notfount = "the file `1` does not exist."

Import::hdr = "The file `1` has an invalid header.";

ExportNii::type = "NumberType should be \"Integer\", \"Real\", \"Complex\", or \"Automatic\"."

OpenMRIcron::notfount = "mricron.exe not found in $UserBaseDirectory or $BaseDirectory please install DTItools in correct directory."

OpenMRIcron::fil = "the file `1` does not exist."


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection::Closed:: *)
(*General Definitions*)


unitsNii = {0 -> "DimensionlessUnit", 1 -> "Meters", 
   2 -> "Millimeters", 3 -> "Micrometers", 8 -> "Seconds", 
   16 -> "Milliseconds", 24 -> "Microseconds", 32 -> "Hertz", 
   40 -> "PartsPerMillion", 48 -> "Radians/Seconds"};

directionsNii = {0 -> Undefined, 1 -> "x", 2 -> "y", 3 -> "z"};

coordinateNii = {0 -> "None", 1 -> "Scanner Posistion", 
   2 -> "Coregistration", 3 -> "Normalized Tal", 
   4 -> "Normalized MNI152ach" 5 -> "Normalized MNI152"};

dataTypeNii = {
   1 -> "Bit",
   2 -> "Integer8",
   4 -> "Integer16",
   8 -> "Integer32",
   16 -> "Real32",
   32 -> "Real32",
   64 -> "Real64",
   128 -> "Unsignedinteger8",
   256 -> "Integer8",
   511 -> "Real32",
   512 -> "Integer16"(*"UnsignedInteger16"*),
   768 -> "Unsignedinteger32",
   1024 -> "Integer64",
   1280 -> "Unsignedinteger64",
   1792 -> "Real64",
   0 -> Undefined
   };

JoinCharacters = StringTrim@StringJoin@StringCases[Cases[#, _?StringQ], character_ /; MemberQ[CharacterRange[" ", "~"], character]]&;


(* ::Subsection::Closed:: *)
(*DcmToNii*)


dcm2nii = $UserBaseDirectory <>"\\Applications\\DTITools\\Applications";
dcm2nii = If[!FileExistsQ[dcm2nii<>"\\dcm2niix.exe"],$BaseDirectory <>"\\Applications\\DTITools\\Applications",dcm2nii];
dcm2nii = If[!FileExistsQ[dcm2nii<>"\\dcm2niix.exe"], Message[DcmToNii::notfount],dcm2nii];

Options[DcmToNii]={CompressNii->True}

SyntaxInformation[DcmToNii] = {"ArgumentsPattern" -> {_.,_.,OptionsPattern[]}};

DcmToNii[opt:OptionsPattern[]]:=DcmToNii["folder","",opt];

DcmToNii[action_?StringQ,opt:OptionsPattern[]] := DcmToNii[action,"",opt]; 

DcmToNii[fstr_?ListQ,opt:OptionsPattern[]] := DcmToNii["folder",fstr,opt];

DcmToNii[action_,fstr_,OptionsPattern[]] := Module[{act,filfolin,folout,add,title,log,command,compress},
	
	Print["Using Chris Rorden's dcm2niix.exe (https://github.com/rordenlab/dcm2niix)"];
	
	(*select the correct message if no file or folder is given*)		
	{act,title}=Switch[action,
		"folder",{"Directory","Select direcotry containig the dcm files"},
		"file",{"FileOpen","Select the dcm file to convert"},
		_,Return[Message[DcmToNii::type]]
		];
	
	(*generate a popup to select the file or folder*)
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
	
	(*create the cmd window command to run dcm2niix.exe*)
	log=" > \"" <> folout <> "\\output.txt";
	add=If[action == "file", "-v N "," "];
	compress=If[OptionValue[CompressNii],"y","n"];
	command="cd " <> dcm2nii <>"\n dcm2niix.exe  -f %f_%i_%m_%n_%p_%q_%s_%t -z "<>compress<>" -o \""<>folout<>"\" \"" <> filfolin <> "\"" <> log<>"\" \n exit \n";
	
	(*perform teh conversion*)
	Monitor[RunProcess[$SystemShell,"StandardOutput",command],ProgressIndicator[Dynamic[Clock[Infinity]], Indeterminate]];
]


(* ::Subsection::Closed:: *)
(*ImportNii*)


(* ::Subsubsection::Closed:: *)
(*Register Import*)


ImportExport`RegisterImport[
  	"Nii",
  {
   "Header" :> ImportNiiHeader,
   "Data" -> ImportNiiInfo,
   "VoxelSize" -> ImportNiiInfo,
   "TR" -> ImportNiiInfo,
   "RotationMatrix" -> ImportNiiInfo,
   "Units" -> ImportNiiInfo,
   "Scaling" -> ImportNiiInfo,
   ImportNiiDefault
   },
  {},
  "AvailableElements" -> {"Data", "Header", "VoxelSize", "TR", 
    "RotationMatrix", "Units", "Scaling"},
  "OriginalChannel" -> True,
  "Options" -> {
    NiiScaling
    }
  ];
  

(* ::Subsubsection::Closed:: *)
(*ConvertNiiExtention*)  


ConvertNiiExtention[file_, channel_, imghdr_] := Switch[imghdr,
  "hdr",
   If[
   StringMatchQ[FileExtension[file], "img", IgnoreCase -> True],
   First[System`ConvertersDump`Decode[
     StringReplace[channel, 
      RegularExpression["img(?!.*img)"] -> "hdr", 
      IgnoreCase -> True], {Automatic}]],
   file
   ],
  "img",
  If[
   StringMatchQ[FileExtension[file], "hdr", IgnoreCase -> True],
   First[System`ConvertersDump`Decode[
     StringReplace[channel, 
      RegularExpression["hdr(?!.*hdr)"] -> "img", 
      IgnoreCase -> True], {Automatic}]],
   file
   ]
  ]

  

(* ::Subsubsection::Closed:: *)
(*GetNiiInformation*) 


GetNiiInformation[hdr_] := 
 Module[{type, size, dim, ddim, offSet, vox, voxU, TR, TRU, slope, intercept, rotmat},
  
  dim = "dim" /. hdr;
  ddim = dim[[1]];
  dim = dim[[2 ;; ddim + 1]];
  size = Times @@ dim;
  
  type = "dataType" /. hdr;
  offSet = Round["voxOffset" /. hdr];
  
  vox = Reverse[("pixDim" /. hdr)[[2 ;; Clip[ddim + 1, {3, 4}]]]];
  TR = ("pixDim" /. hdr)[[5]];
  {slope, intercept} = {"scaleSlope", "scaleInteger"} /. hdr;
  
  rotmat = ({1, 1, -1} {"sRowx", "sRowy", "sRowz"} /. hdr)[[All, 1 ;; 3]]/ConstantArray[Reverse[vox], 3];
  rotmat = DiagonalMatrix[{1, -1, 1}].ConstantArray[Diagonal[Sign[Sign[rotmat] + 0.0000001]], 3] rotmat;
  
  {voxU, TRU} = "xyztUnits" /. hdr;
  
  {{type, size, dim, offSet},{vox, voxU, TR, TRU, slope, intercept, rotmat, ddim}}
]


(* ::Subsubsection::Closed:: *)
(*ImportNiiDefault*)  


Options[ImportNiiDefault] = {"Channel" -> Null, 
   "ExtensionParsing" -> False, NiiScaling -> False, NiiInfo -> False};

ImportNiiDefault[file_, opts : OptionsPattern[]] := 
 Module[{hdr, data, byteOrder, dataInfo, info, scaling},
  hdr = ImportNiiHeader[file, opts];
  byteOrder = "ByteOrder" /. hdr;
  hdr = "Header" /. hdr;
  
  If[hdr === $Failed, Return[$Failed, Module]];
  
  (*get all the data and info*)
  {dataInfo, info} = GetNiiInformation[hdr];
  data = ImportNiiData[file, dataInfo, byteOrder, opts];
  
  (*flip dimensions*)
  If[info[[8]] === 4, data = Transpose[data, {2, 1, 3, 4}]];
  data = Map[Reverse[#] &, data, {info[[8]] - 2}];
  If[Positive[("sRowx" /. hdr)[[1]]], 
   data = Map[Reverse[#] &, data, {info[[8]] - 1}]];
  
  (*scale data*)
  scaling = info[[6]] + info[[5]] # &;
  If[OptionValue[NiiScaling], data = scaling[data]];
  
  If[OptionValue[NiiInfo],
   {
    "Data" -> data,
    "VoxelSize" -> info[[1]],
    "TR" -> info[[3]],
    "RotationMatrix" -> info[[7]],
    "Units" -> info[[{2, 4}]],
    "Scaling" -> info[[5 ;; 6]]
    }
   ,
   {data, info[[1]]}
   ]
  ]


(* ::Subsubsection::Closed:: *)
(*ImportNiiInfo*)  


ImportNiiInfo[file_, opts : OptionsPattern[ImportNiiDefault]] := ImportNiiDefault[file, NiiInfo -> True, opts]


(* ::Subsubsection::Closed:: *)
(*ImportNiiData*)  


ImportNiiData[file_, {type_, size_, dim_, off_}, byteorder_, 
  OptionsPattern[ImportNiiDefault]] := Module[{datafile, strm, data},
  
  datafile = ConvertNiiExtention[file, OptionValue["Channel"], "img"];
  
  strm = OpenRead[datafile, BinaryFormat -> True];
  SetStreamPosition[strm, off];
  data = BinaryReadList[strm, type, size, ByteOrdering -> byteorder];
  Close[strm];
  
  data = ArrayReshape[data, Reverse@dim];
  
  data
  ]


(* ::Subsubsection::Closed:: *)
(*ImportNiiData*)  


ImportNiiHeader[file_, OptionsPattern[ImportNiiDefault]] := 
 Module[{strm, verNii,hdrfile, byteorder, hdrValues, hdr, val},
  
  hdrfile = ConvertNiiExtention[file, OptionValue["Channel"], "hdr"];
  byteorder = $ByteOrdering;
  
  strm = OpenRead[hdrfile, BinaryFormat -> True];
  
  (*determine version*)
  verNii = Switch[Quiet[BinaryRead[strm, "Integer32", ByteOrdering -> byteorder]], 
     348, 1, 540, 2, _, Undefined
     ];
  (*check for bite ordering*)
  If[verNii === Undefined,
   byteorder *= -1;
   SetStreamPosition[strm, 0];
   verNii = Switch[Quiet[BinaryRead[strm, "Integer32", ByteOrdering -> byteorder]], 
   	348, 1, 540, 2, _, Undefined];
   ];
 
	(*reset stream and read header in version is defined*)
	 SetStreamPosition[strm, 0];
	 
	 hdrValues = Switch[verNii,
	    1,
	    {
	     {"size", "Integer32", 1},
	     {"data_Type", "Character8", 10},
	     {"dbName", "Character8", 18},
	     {"extents", "Integer32", 1},
	     {"sessionError", "Integer16", 1},
	     {"regular", "Character8", 1},
	     {"dimInfo", "UnsignedInteger8", 1},
	     {"dim", "Integer16", 8},
	     {"intentP1", "Real32", 1},
	     {"intentP2", "Real32", 1},
	     {"intentP3", "Real32", 1},
	     {"intentCode", "Integer16", 1},
	     {"dataType", "Integer16", 1},
	     {"bitPix", "Integer16", 1},
	     {"sliceStart", "Integer16", 1},
	     {"pixDim", "Real32", 8},
	     {"voxOffset", "Real32", 1},
	     {"scaleSlope", "Real32", 1},
	     {"scaleInteger", "Real32", 1},
	     {"sliceEnd", "Integer16", 1},
	     {"sliceCode", "UnsignedInteger8", 1},
	     {"xyztUnits", "UnsignedInteger8", 1},
	     {"calMax", "Real32", 1},
	     {"calMin", "Real32", 1},
	     {"sliecDuration", "Real32", 1},
	     {"tOffset", "Real32", 1},
	     {"glMax", "Integer32", 1},
	     {"glMin", "Integer32", 1},
	     {"descrip", "Character8", 80},
	     {"auxFile", "Character8", 24},
	     {"qformCode", "Integer16", 1},
	     {"sformCode", "Integer16", 1},
	     {"quaternB", "Real32", 1},
	     {"quaternC", "Real32", 1},
	     {"quaternD", "Real32", 1},
	     {"qOffsetX", "Real32", 1},
	     {"qOffsetY", "Real32", 1},
	     {"qOffsetZ", "Real32", 1},
	     {"sRowx", "Real32", 4},
	     {"sRowy", "Real32", 4},
	     {"sRowz", "Real32", 4},
	     {"intentName", "Character8", 16},
	     {"magic", "Character8", 4},
	      {"ECode", "Byte", 4}
	     }
	    ,
	    2,
	    {
	     {"size", "Integer32", 1},
	     {"magic", "Character8", 8},
	     {"dataType", "Integer16", 1},
	     {"bitPix", "Integer16", 1},
	     {"dim", "Integer64", 8},
	     {"intentP1", "Real64", 1},
	     {"intentP2", "Real64", 1},
	     {"intentP3", "Real64", 1},
	     {"pixDim", "Real64", 8},
	     {"voxOffset", "Integer64", 1},
	     {"scaleSlope", "Real64", 1},
	     {"scaleInteger", "Real64", 1},
	     {"calMax", "Real64", 1},
	     {"calMin", "Real64", 1},
	     {"sliecDuration", 1},
	     {"tOffset", "Real64", 1},
	     {"sliceStart", "Integer64", 1},
	     {"sliceEnd", "Integer64", 1},
	     {"descrip", "Character8", 80},
	     {"auxFile", "Character8", 24},
	     {"qformCode", "Integer32", 1},
	     {"sformCode", "Integer32", 1},
	     {"quaternB", "Real64", 1},
	     {"quaternC", "Real64", 1},
	     {"quaternD", "Real64", 1},
	     {"qOffsetX", "Real64", 1},
	     {"qOffsetY", "Real64", 1},
	     {"qOffsetZ", "Real64", 1},
	     {"sRowX", "Real64", 4},
	     {"sRowY", "Real64", 4},
	     {"sRowZ", "Real64", 4},
	     {"sliceCode", "Integer32", 1},
	     {"xyztUnits", "Integer32", 1},
	     {"intentCode", "Integer32", 1},
	     {"intentName", "Character8", 16},
	     {"dimInfo", "UnsignedInteger8", 1},
	     {"unusedString", "Character8", 15},
	     {"ECode", "Byte", 4}
	     }
	    ,
	    _, Message[Import::hdr, file]; Return[$Failed, Module]
	    ];
	 
	 (*read the header and convert all values to readable text and numbers*)
	 hdr = (val = BinaryReadList[strm, #[[2]], #[[3]], ByteOrdering -> byteorder];
	      #[[1]] -> If[#[[3]] === 1, First@val, If[#[[2]] === "Character8", JoinCharacters@val, val]]
	      ) & /@ hdrValues;
	 
	 (*check if header is valid*)
	 If[hdr[[-1, -1]] === EndOfFile, Message[Import::hdr, file]; Return[$Failed, Module]];
	 
	 (*replace values with readable text*)
	 hdr = hdr /. {("dimInfo" -> ___) -> ("dimInfo" -> (FromDigits[Reverse[#], 2] & /@ Partition[Reverse[IntegerDigits["dimInfo" /. hdr, 2, 8]], 2] /. directionsNii))};
	 hdr = hdr /. {("xyztUnits" -> ___) -> ("xyztUnits" -> (FromDigits[Reverse[#], 2] & /@ Partition[Reverse[IntegerDigits["xyztUnits" /. hdr, 2, 8]][[;; 6]], 3] {1, 8} /. unitsNii))};
	 hdr = hdr /. {("dataType" -> ___) -> ("dataType" -> (("dataType" /. hdr) /. dataTypeNii))};
	 hdr = hdr /. {("qformCode" -> ___) -> ("qformCode" -> (("qformCode" /. hdr) /. coordinateNii))};
	 hdr = hdr /. {("sformCode" -> ___) -> ("sformCode" -> (("sformCode" /. hdr) /. coordinateNii))};
	 
	 Close[strm];
	 
	 {"Header" -> hdr, "ByteOrder" -> byteorder}
 ]


(* ::Subsubsection::Closed:: *)
(*ImportNii*)  


Options[ImportNii] = {NiiMethod -> "default", NiiScaling -> True};

SyntaxInformation[ImportNii] = {"ArgumentsPattern" -> {_., OptionsPattern[]}};

ImportNii[opts : OptionsPattern[]] := ImportNii["", opts];

ImportNii[fil_String: "", OptionsPattern[]] := Module[
	{file,what, out},
	
	what = OptionValue[NiiMethod];
  
  If[! MemberQ[{"default", "data", "dataTR", "header", "scaling", 
      "headerMat", "rotation", "all"}, what], 
   Return[Message[ImportNii::wht]]];
  
  (*select a file if none was given*)
  file = If[fil == "", FileSelect["FileOpen", {"*.nii", "*.hdr", "*.nii.gz"}, "nifti files ", WindowTitle -> "Select the nii file to import"], fil];
  If[file == Null || file === $Canceled, Return[]];
  
  Switch[what,
   "data", Import[fil, {"nii", "Data"}],
   "header", Import[fil, {"nii", {"Data", "VoxelSize", "Header"}}],
   "headerMat", Import[fil, {"nii", {"Data", "VoxelSize", "Header", "RotationMatrix"}}],
   "dataTR", Import[fil, {"nii", {"Data", "VoxelSize", "TR"}}], 
   "rotation", Import[fil, {"nii", {"Data", "VoxelSize", "RotationMatrix"}}],
   "scaling", Import[fil, {"nii", {"Data", "VoxelSize", "Scaling"}}],
   "all", out = Import[fil, {"nii", {"Data", "VoxelSize", "Header", "TR", "Scaling", "RotationMatrix", "Units"}}]; 
   	{out[[1]], out[[2]], out[[3 ;;]]}, 
   _, Import[fil, "nii"]]
  
  ]


(* ::Subsection::Closed:: *)
(*ImportBval*)


SyntaxInformation[ImportBval] = {"ArgumentsPattern" -> {_.}};

ImportBval[]:=ImportBval[FileSelect["FileOpen", {".bval"}, WindowTitle -> "Select *.bval"]]

ImportBval[file_]:=Flatten[LineToList[file] ]//N


(* ::Subsection::Closed:: *)
(*ImportBval*)


Options[ImportBvec]={FlipBvec->False};

SyntaxInformation[ImportBvec] = {"ArgumentsPattern" -> {_.,OptionsPattern[]}};

ImportBvec[file_?StringQ, OptionsPattern[]]:=Module[{grads},
	grads=Round[LineToList[file],0.0001];
	grads=If[OptionValue[FlipBvec],
		{1, -1, 1}#&/@grads,
		{1,-1,1}RotateLeft[#]&/@grads
		];
	If[Negative[#[[3]]],-#,#]&/@grads
]

ImportBvec[opts:OptionsPattern[]]:=ImportBvec[FileSelect["FileOpen", {".bvec"}, WindowTitle -> "Select *.bvec"], opts]


(* ::Subsection::Closed:: *)
(*ImportBvalvec*)


Options[ImportBvalvec]={FlipBvec->False};

SyntaxInformation[ImportBvalvec] = {"ArgumentsPattern" -> {_.,_.,OptionsPattern[]}};

ImportBvalvec[file_?StringQ,opts:OptionsPattern[]] := Module[{valf, vecf},
  valf = FileBaseName[file]<>".bval";
  vecf = FileBaseName[file]<>".bvec";
  ImportBvalvec[valf,vecf,opts]
  ]

ImportBvalvec[valf_,vecf_,opts:OptionsPattern[]] := {ImportBval[valf],ImportBvec[vecf,opts]}

ImportBvalvec[___,opts:OptionsPattern[]] := Module[{valf, vecf},
  valf = FileSelect["FileOpen", {"*.bval"}, WindowTitle -> "Select *.bval"];
  vecf = FileSelect["FileOpen", {"*.bvec"}, WindowTitle -> "Select *.bvec"];
  ImportBvalvec[valf,vecf,opts]
  ]

LineToList[file_] := Module[{grads,tmp},
  grads = Import[file, "Lines"];
  grads = Transpose[(
  	tmp = StringReplace[#, {" " -> ",", "\t" -> ",", "E" -> "*10^", "e" -> "*10^"}];
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


Options[ImportNiiDiff]={RotateGradients->False,FlipBvec->True}

SyntaxInformation[ImportNiiDiff]= {"ArgumentsPattern" -> {_.,_.,_.,OptionsPattern[]}};

ImportNiiDiff[OptionsPattern[]]:=Module[{data,grad,bvec,vox,hdr,mat},
	{data,vox,hdr,mat}=ImportNii[NiiMethod -> "headerMat"];
	{bvec, grad}=ImportBvalvec[FlipBvec->OptionValue[FlipBvec]];
	{data,Round[If[OptionValue[RotateGradients],grad.Inverse[mat], grad],.0001],bvec,vox}
]

ImportNiiDiff[file_String,OptionsPattern[]]:=Module[{data,grad,bvec,vox,hdr,mat},
	{data,vox,hdr,mat}=ImportNii[file,NiiMethod -> "headerMat"];
	{bvec, grad}=ImportBvalvec[StringDrop[file,-4],FlipBvec->OptionValue[FlipBvec]];
	{data,Round[If[OptionValue[RotateGradients],grad.Inverse[mat], grad],.0001],bvec,vox}
]

ImportNiiDiff[fnii_String,fvec_String,fval_String,OptionsPattern[]]:=Module[{data,grad,bvec,vox,hdr,mat},
	{data,vox,hdr,mat}=ImportNii[fnii,NiiMethod -> "headerMat"];
	{bvec, grad} = ImportBvalvec[fval, fvec,FlipBvec->OptionValue[FlipBvec]];
	{data,Round[If[OptionValue[RotateGradients],grad.Inverse[mat], grad],.0001],bvec,vox}
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

ExportNii[dato_, voxi_, opts:OptionsPattern[]] := ExportNii[dato, voxi, "" ,opts]

ExportNii[dato_, voxi_, fil_,OptionsPattern[]] := Block[{datao, type, dim, depth, dimo, voxo, srowx, srowy, srowz, 
	ftype, echo, strm, fileo, typeo, vox},
  
  fileo = If[fil == "", FileSelect["FileSave",{"*.nii"},"nifti",WindowTitle->"Select the destination file"], fil];
  If[fileo == Null, Return[]];
  
  ftype=OptionValue[NumberType];
  
  dim = Dimensions[dato];
  depth = ArrayDepth[dato];
  
  If[ListQ[voxi[[1]] && NumberQ[voxi[[2]]]],
 		vox = voxi[[1]]; echo = voxi[[2]];
 		,
 		vox = voxi;echo=0;];
  
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
  voxo[[5]]=echo;
  
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
  (*xyztUnits*)BinaryWrite[strm, If[echo!=0,10,2], "UnsignedInteger8"];(*0-unknown 1-meters 2-millimeter 3-micrometers 10-mm/sec*)
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


(* ::Subsection:: *)
(*OpenMRIcron*)


(* ::Subsubsection::Closed:: *)
(*OpenMRIcron*)


Options[OpenMRIcron] = {ImportResult -> True, NumberOfResults -> 1};

SyntaxInformation[OpenMRIcron] = {"ArgumentsPattern" -> {_., OptionsPattern[]}};

OpenMRIcron[opts : OptionsPattern[]] := OpenMRIcron["", opts];
OpenMRIcron[filei_, OptionsPattern[]] := 
 Block[{file, mricron, command, num},
  mricron = FindMRIcron[];
  If[! (mricron === Null),
   file = If[filei === "", FileSelect["FileOpen"], filei] /. Null -> "Cancel";
   
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


(* ::Subsubsection::Closed:: *)
(*FindMRIcron*)


FindMRIcron[] := Module[{mricron},
  mricron = $UserBaseDirectory <> 
    "\\Applications\\DTITools\\Applications\\mricron.exe";
  mricron = 
   If[! FileExistsQ[mricron], $BaseDirectory <> 
     "\\Applications\\DTITools\\Applications\\mricron.exe", mricron];
  If[! FileExistsQ[mricron], Return[Message[OpenMRIcron::notfount]], 
   mricron]
  ]


(* ::Subsection::Closed:: *)
(*ExtractNiiFiles*)


ExtractNiiFiles[] := ExtractNiiFiles[FileSelect["Directory", WindowTitle -> "Select direcotry containig the nii files"]]
ExtractNiiFiles[folder_] := Module[{files},
	Quiet[
  files = FileNames["*.nii.gz", folder];
  DeleteFile[StringDrop[#, -3]] & /@ files;
  ExtractArchive /@ files;
  DeleteFile /@ files;
  ]
  ]


(* ::Subsection::Closed:: *)
(*CompressNiiFiles*)

CompressNiiFiles[] := CompressNiiFiles[FileSelect["Directory", WindowTitle -> "Select direcotry containig the nii files"]]
CompressNiiFiles[folder_] := Module[{files},
  Quiet[
   files = FileNames["*.nii", folder];
   
   CreateArchive[#, # <> ".gz"] & /@ files;
   DeleteFile /@ files;
   ]
  ]


(* ::Section:: *)
(*End Package*)


End[]

If[DTITools`verbose,Print[Names["DTITools`NiftiTools`*"]]];
SetAttributes[#,{Protected, ReadProtected}]&/@Names["DTITools`NiftiTools`*"];

EndPackage[]
