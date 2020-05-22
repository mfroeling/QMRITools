(* ::Package:: *)

(* ::Title:: *)
(*QMRITools NiftiTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`NiftiTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`NiftiTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


DcmToNii::usage =
"DcmToNii[] converts a dicom folder to nii, you will be prometed for the location of the folders. 
DcmToNii[{\"input\",\"ouput\"}] converts the \"input\" dicom folder to nii files which are place in the \"output\" folder.
For this function to work the dcm2niix.exe file should be present in the QMRITools aplication folder."


ImportNii::usage = 
"ImportNii[] promts to select the nii file to import.
ImportNii[\"file\"] imports the nii file. 
The default output is {data, vox}, however using NiiMethod various outputs can be given.
The Nii import is also suported using the native Import function from Mathematica."


ImportNiiDiff::usage = 
"ImportNiiDiff[] will promt for the *.nii, *.bvec and *.bval file to import.
ImportNiiDiff[*.nii] will import the *.nii file and automatically also imports the *.bvec and *.bval is they have the same name.
ImportNiiDiff[*.nii,*.bvec,*.bval] will import the given files.
The output will be {data,grad,bvec,vox}."

ImportNiiDix::usage = 
"ImportNiiDix[\"file\"] imports the dixon nii file which should contain all possible outputs given by the scanner and corrects them accordingly."

ImportNiiT2::usage = 
"ImportNiiT2[\"file\"] imports the T2 file which should contain the echos and the T2map calculated by the scanner and corrects them accordingly."

ImportNiiT1::usage = 
"ImportNiiT1[\"file\"] imports the T1 file which should contain the echos and the T1map calculated by the scanner and corrects them accordingly."

ImportExploreDTItens::usage = 
"ImportExploreDTItens[\"file\"] imports the *.nii export for the tensor from explore DTI."


ImportBvalvec::usage =
"ImportBvalvec[] will promt to select the *.bval and *.bvec files.
ImportBvalvec[file] if file is either a *.bval or *.bvec it will automatically import the *.bval and *.bvec files.
ImportBvalvec[*.bvec,*.bval] imports the given *.bval and *.bvec files." 

ImportBval::usage = 
"ImportBval[] will promt to select the *.bval file.
ImportBval[*.bval] imports the given *.bval file." 

ImportBvec::usage = 
"ImportBvec[] will promt to select the *.bvec file.
ImportBvec[*.bvec] imports the given *.bvec file." 

ImportBmat::usage =
"ImportBmat[] will promt to select the *.txt file containing the bmatrix.
ImportBmat[*.txt] imports the given *.txt file containing the bmatrix." 


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


ExtractNiiFiles::usage =
"ExtractNiiFiles[] promts for a folder. It then extracts all nii.gz files to .nii files in the selected folder.
ExtractNiiFiles[folder] extracts all nii.gz files to .nii files in folder."

CompressNiiFiles::usage =
"CompressNiiFiles[] promts for a folder. It then compresses all nii files to .nii.gz files in the selected folder.
CompressNiiFiles[folder] compresses all nii files to .nii.gz files in folder."


(* ::Subsection::Closed:: *)
(*Options*)


NiiMethod::usage = "NiiMethod is an option for ImportNIi. Values can be \"data\", \"dataTR\", \"header\", \"scaling\", \"headerMat\", \"rotation\", \"all\"."

NiiScaling::usage = "NiiScaling is an option for ImportNii. It scales the nii values with scale slope and offset for quantitative data."

CompressNii::usage = "CompressNii is an option for DcmToNii and ExportNii. If set True .nii.gz files will be created."

NiiDataType::usage = "NiiDataType is an option of Export Nii. The number type of Nii file can be \"Integer\", \"Real\", \"Complex\", or \"Automatic\"."

RotateGradients::usage = "RotateGradients is an option for ImportNiiDiff."

FlipBvec::usage = "FlipBvec is an option for ImportBvalvec."


(* ::Subsection::Closed:: *)
(*Error Messages*)


DcmToNii::notfount = "dcm2nii.exe not found in $UserBaseDirectory or $BaseDirectory please install QMRITools in correct directory."

DcmToNii::type = "Input should be \"file\" or \"folder\"."

ImportNii::wht = "should be \"data\", \"dataTR\", \"header\", \"scaling\", \"headerMat\", \"rotation\", \"all\"."

ImportNii::notfount = "the file `1` does not exist."

Import::niihdr = "The file `1` has an invalid header.";

Export::niitype= "The type of the specified data is incompatible with the specified `1`.";

Export::niiran="The range of the specified data `2` is incompatible with the specified `1`";

Export::niidat="The data should be an array of numbers.";

Export::niihdr="The given header is invalid.";

Export::niidim="`1` compatible with the dimensions of the input data";

ExportNii::type = "NiiDataType should be \"Integer\", \"Real\", \"Complex\", or \"Automatic\"."


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection::Closed:: *)
(*DcmToNii*)


Options[DcmToNii]={CompressNii->True, Method->Automatic}

SyntaxInformation[DcmToNii] = {"ArgumentsPattern" -> {_.,_.,OptionsPattern[]}};

DcmToNii[opt:OptionsPattern[]]:=DcmToNii[{"",""},opt];

DcmToNii[{infol_?StringQ,outfol_?StringQ},OptionsPattern[]] := Module[{filfolin,folout, log,command,compress,dcm2nii},
	
	dcm2nii=FindDcm2Nii[];
	If[dcm2nii==$Failed,Return[$Failed,Module]];
	
	Print["Using Chris Rorden's dcm2niix.exe (https://github.com/rordenlab/dcm2niix)"];

	(*generate a popup to select the file or folder*)
	filfolin=If[infol=="",FileSelect["Directory",WindowTitle->"Select direcotry containig the dcm files"],infol];
	If[filfolin==Null||folout==Null,Return[]];
	folout=If[outfol=="",FileSelect["Directory",WindowTitle->"Select directory to put nii files in"],outfol];
	If[filfolin==Null||folout==Null,Return[]];
	
	compress=If[OptionValue[CompressNii],"i","n"];
	
	(*create the cmd window command to run dcm2niix*)
	log=FileNameJoin[{folout,"output.txt"}];
	
	command = Switch[$OperatingSystem,
		"Windows",
		First@FileNameSplit[dcm2nii]<>"\ncd " <> dcm2nii <>"\ndcm2niix.exe  -f %f_%s_%t_%i_%m_%n_%p_%q -z "<>
		compress<>" -o \""<>folout<>"\" \""<> filfolin<>"\" > \""<>log<>"\nexit\n"
		,
		"Unix",
		dcm2nii<>"dcm2niix -f %f_%s_%t_%i_%m_%n_%p_%q -z "<>
		compress<>" -o '"<>folout<>"' '"<>filfolin<>"' > '"<>log<>"'\nexit\n"
		,
		"MacOSX",
		dcm2nii<>"dcm2niix -f %f_%s_%t_%i_%m_%n_%p_%q -z "<>
		compress<>" -o '"<>folout<>"' '"<>filfolin<>"' > '"<>log<>"'\nexit\n"
	];
	
	If[OptionValue[Method]=!=Automatic,Print[command]];
	
	(*perform teh conversion*)
	Monitor[
		RunProcess[$SystemShell,"StandardOutput",command],
		ProgressIndicator[Dynamic[Clock[Infinity]], Indeterminate]];
]


FindDcm2Nii[]:=Module[{fil1,fil2},
	Switch[$OperatingSystem,
		"Windows",
		fil1=$UserBaseDirectory <>"\\Applications\\QMRITools\\Applications\\windows-x86-64\\dcm2niix.exe";
		fil2=$BaseDirectory <>"\\Applications\\QMRITools\\Applications\\windows-x86-64\\dcm2niix.exe";
		,
		"Unix",
		fil1=$UserBaseDirectory <>"/Applications/QMRITools/Applications/Linux-x86-64/bin/dcm2niix";
		fil2=$BaseDirectory <>"/Applications/QMRITools/Applications/Linux-x86-64/bin/dcm2niix";
		,
		"MacOSX",
		fil1=$UserBaseDirectory <>"/Applications/QMRITools/Applications/MacOSX-x86-64/bin/dcm2niix";
		fil2=$BaseDirectory <>"/Applications/QMRITools/Applications/MacOSX-x86-64/bin/dcm2niix";
	];
	
	If[FileExistsQ[fil1],DirectoryName[fil1],If[FileExistsQ[fil2],DirectoryName[fil2], Message[DcmToNii::notfount];$Failed]]
]


(* ::Subsection:: *)
(*General Nii Functions*)


(* ::Subsubsection::Closed:: *)
(*General Nii Definitions*)


unitsNii = {0 -> "DimensionlessUnit", 1 -> "Meters", 
   2 -> "Millimeters", 3 -> "Micrometers", 8 -> "Seconds", 
   16 -> "Milliseconds", 24 -> "Microseconds", 32 -> "Hertz", 
   40 -> "PartsPerMillion", 48 -> "Radians/Seconds"};


directionsNii = {0 -> Undefined, 1 -> "x", 2 -> "y", 3 -> "z"};


coordinateNii = {0 -> "None", 1 -> "Scanner Posistion", 
   2 -> "Coregistration", 3 -> "Normalized Tal", 
   4 -> "Normalized MNI152ach" 5 -> "Normalized MNI152"};


dataTypeNii = {
   0 -> Undefined,
   1 -> "Bit",
   2 -> "Byte",
   4 -> "Integer16",
   8 -> "Integer32",
   16 -> "Real32",
   32 -> "Complex64",
   64 -> "Real64",
   128 -> {"Byte", "Byte", "Byte"}(*RGBColor*),
   256 -> "Integer8",
   512 -> "UnsignedInteger16",
   768 -> "Unsignedinteger32",
   1024 -> "Integer64",
   1280 -> "Unsignedinteger64",
   1536 -> "Real128",
   1792 -> "Complex128",
   2048 -> "Complex256"
   };


typeSizeNii = {0 -> Undefined, 2->8, 4->16, 8->32, 16->32, 32->64, 64->64, 128->24, 256->8, 512->16, 768->32, 1024->32, 1280->32, 1536->128, 1792->128, 2048->256}; 


rangeNii = {
   	Undefined -> Undefined,
   	"Byte" -> {0, 255},
   	"Integer16" -> {-32768, 32767},
   	"Integer32" -> {-2147483648, 2147483647},
   	"Real32" -> Undefined,
   	"Complex64" -> Undefined,
   	"Real64" -> Undefined, 
   	{"Byte", "Byte", "Byte"} -> {0, 255},(*RGBColor*)
   	"Integer8" -> {-128, 127},
   	"UnsignedInteger16" -> {0, 65535},
   	"UnsignedInteger32" -> {0, 4294967295},
   	"Integer64" -> {-9223372036854775808, 9223372036854775807},
   	"UnsignedInteger64" -> {0, 18446744073709551615},
   	"Real128" -> Undefined,
   	"Complex128" -> Undefined,
   	"Complex256" -> Undefined
   };


typeCheckNii = Join[
   Thread[{Undefined, "Complex64", "Complex128", "Complex256"} ->  NumberQ],
   Thread[{"Byte" , "Integer8", "Integer16", "Integer32", "Integer64", "UnsignedInteger8", "UnsignedInteger16", "UnsignedInteger32", "UnsignedInteger64"} -> IntegerQ],
   Thread[{"Real32", "Real64", "Real128"} -> Internal`RealValuedNumberQ]
   ];


(* ::Subsubsection::Closed:: *)
(*IntegerChop*)


IntegerChop = # + Chop[#2 - #] &[Round@#, #] &;


(* ::Subsubsection::Closed:: *)
(*JoinCharacters*)


JoinCharacters = StringTrim@StringJoin@StringCases[Cases[#, _?StringQ], character_ /; MemberQ[CharacterRange[" ", "~"], character]]&;


(* ::Subsubsection::Closed:: *)
(*GetNiiHeaderValues*)


GetNiiHeaderValues[verNii_]:=GetNiiHeaderValues[verNii,""]
GetNiiHeaderValues[verNii_,file_] := Switch[verNii,
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
	    _, 
	    If[file==="", Message[Export::niiver], Message[Export::niihdr, file]]; $Failed
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
(*ReplaceHeaderRule*)


ReplaceHeaderVal[hdr_, val_, rep_] := hdr /. {(val -> ___) -> (val -> rep)};
ReplaceHeaderRule[hdr_, val_, rep_] := hdr /. {(val -> ___) -> (val -> ((val /. hdr) /. rep))};


(* ::Subsubsection::Closed:: *)
(*DimInfo*)


DimInfo[inp_?IntegerQ] := ((FromDigits[Reverse[#], 2] & /@ Partition[Reverse[IntegerDigits[inp, 2, 8]], 2]) /. directionsNii)
DimInfo[inp_?ListQ] := FromDigits[Flatten@Reverse[IntegerDigits[#, 2, 2] & /@ (inp /. Reverse[directionsNii, 2])], 2]


(* ::Subsubsection::Closed:: *)
(*XyztUnits*)


XyztUnits[inp_?IntegerQ] := (FromDigits[Reverse[#], 2] & /@ Partition[Reverse[IntegerDigits[inp, 2, 8]][[;; 6]], 3] {1, 8} /. unitsNii)
XyztUnits[inp_?ListQ] := FromDigits[PadLeft[Flatten[Reverse[IntegerDigits[#, 2, 3] & /@ ((inp /. Reverse[unitsNii, 2])/{1, 8})]], 8], 2]


(* ::Subsubsection::Closed:: *)
(*ArrangeData*)


ArrangeData[data_] := ArrangeData[data, ArrayDepth[data]]

ArrangeData[data_, depth_] := 
 Flatten[If[depth == 4, 
   Transpose[Reverse[Reverse[data, depth], depth - 1]], 
   Reverse[Reverse[data, depth], depth - 1]]]


(* ::Subsubsection::Closed:: *)
(*RemoveExtention*)


RemoveExtention[fil_] := Module[{file},
  file = FileNameSplit[fil];
  file[[-1]] = FileBaseName[file[[-1]]];
  FileNameJoin[file]
  ]


(* ::Subsection:: *)
(*ImportNii*)


(* ::Subsubsection::Closed:: *)
(*Register Import*)


ImportExport`RegisterImport["Nii",
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
	"AvailableElements" -> {"Data", "Header", "VoxelSize", "TR", "RotationMatrix", "Units", "Scaling"},
	"OriginalChannel" -> True,
	"Options" -> {NiiScaling}
];


(* ::Subsubsection::Closed:: *)
(*ImportNii*)


Options[ImportNii] = {NiiMethod -> "default", NiiScaling -> False};

SyntaxInformation[ImportNii] = {"ArgumentsPattern" -> {_., OptionsPattern[]}};

ImportNii [opts : OptionsPattern[]] := ImportNii["", opts];

ImportNii[fil_String: "", OptionsPattern[]] := Module[{file,what, out,output,opt},
	
	what = OptionValue[NiiMethod];
	
	If[! MemberQ[{"default", "data", "dataTR", "header", "scaling", "headerMat", "rotation", "all"}, what], 
		Return[Message[ImportNii::wht];$Failed,Module]
		];

	(*select a file if none was given*)
	file = If[fil == "",
		FileSelect["FileOpen", {"*.nii.gz", "*.nii"}, "nifti files ", WindowTitle -> "Select the nii file to import"],
		(*chekc if given file exists if not check for .gz version*)
		If[FileExistsQ[fil], 
			fil, 
			file = fil <> ".gz";
			If[FileExistsQ[file], file, $Failed]
		]
	];
	
	(*stop if ther is no file*)
	If[file == Null || file === $Canceled || file === $Failed, Message[Import::nffil,fil];Return[$Failed,Module]];
	
	opt = NiiScaling->OptionValue[NiiScaling];
	
	output=Switch[what,
		"data", Import[file, {"nii", "Data"}, opt],
		"header", Import[file, {"nii", {"Data", "VoxelSize", "Header"}}, opt],
		"headerMat", Import[file, {"nii", {"Data", "VoxelSize", "Header", "RotationMatrix"}}, opt],
		"dataTR", Import[file, {"nii", {"Data", "VoxelSize", "TR"}}, opt], 
		"rotation", Import[file, {"nii", {"Data", "VoxelSize", "RotationMatrix"}}, opt],
		"scaling", Import[file, {"nii", {"Data", "VoxelSize", "Scaling"}}, opt],
		"all", out = Import[file, {"nii", {"Data", "VoxelSize", "Header", "TR", "Scaling", "RotationMatrix", "Units"}}, opt];
		{out[[1]], out[[2]], out[[3 ;;]]},
		_, Import[file, "nii", opt]
	];
	
	
	(*remove temp files for gz format*)
	Quiet[DeleteFile[FileNames["*" <> FileBaseName[Last[FileNameSplit[file]]], $TemporaryDirectory]]];
	
	output
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
  
  vox = Reverse[("pixDim" /. hdr)[[2 ;; 4]]];
  TR = ("pixDim" /. hdr)[[5]];
  {slope, intercept} = {"scaleSlope", "scaleInteger"} /. hdr;
  
  rotmat = ({1, 1, -1} {"sRowx", "sRowy", "sRowz"} /. hdr)[[All, 1 ;; 3]]/ConstantArray[Reverse[vox], 3];
  rotmat = DiagonalMatrix[{1, -1, 1}].ConstantArray[Diagonal[Sign[Sign[rotmat] + 0.0000001]], 3] rotmat;
  
  {voxU, TRU} = "xyztUnits" /. hdr;
  
  {{type, size, dim, offSet},{vox, voxU, TR, TRU, slope, intercept, rotmat, ddim}}
]


(* ::Subsubsection::Closed:: *)
(*ImportNiiDefault*)


Options[ImportNiiDefault] = {"Channel" -> Null, "ExtensionParsing" -> False, NiiScaling -> False, NiiInfo -> False, NiiFlip->True};

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
  If[OptionValue[NiiFlip],
  	If[Positive[("sRowx" /. hdr)[[1]]+10^-10], 
  		data = Map[Reverse[#] &, data, {info[[8]] - 1}]
  		];
  	,
  	data = Map[Reverse[#] &, data, {info[[8]] - 1}];
  ];
  
  (*scale data*)
  scaling = If[info[[1]] =!= {"Byte", "Byte", "Byte"}, info[[6]] + info[[5]] # &, Identity];
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


ImportNiiData[file_, {type_, size_, dim_, off_}, byteorder_, OptionsPattern[ImportNiiDefault]] := Module[{datafile, strm, data},
  
  datafile = ConvertNiiExtention[file, OptionValue["Channel"], "img"];
  
  strm = OpenRead[datafile, BinaryFormat -> True];
  SetStreamPosition[strm, off];
  data = BinaryReadList[strm, type, size, ByteOrdering -> byteorder];
  Close[strm];
  
  data = ArrayReshape[data, Reverse@dim];
  
  ToPackedArray[N@data]
  ]


(* ::Subsubsection::Closed:: *)
(*ImportNiiHeader*)


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
	 
	 hdrValues = GetNiiHeaderValues[verNii,file];
	 If[hdrValues===$Failed, Return[$Failed, Module]];
	 
	 (*read the header and convert all values to readable text and numbers*)
	 hdr = (
	 	val = BinaryReadList[strm, #[[2]], #[[3]], ByteOrdering -> byteorder];
	    #[[1]] -> If[#[[3]] === 1, First@val, If[#[[2]] === "Character8", JoinCharacters@val, val]]
	  ) & /@ hdrValues;
	 
	 (*check if header is valid*)
	 If[hdr[[-1, -1]] === EndOfFile, Message[Import::niihdr, file]; Return[$Failed, Module]];
	 
	 (*replace values with readable text*)
	 hdr = ReplaceHeaderVal[hdr,"dimInfo",DimInfo["dimInfo" /. hdr]];
	 hdr = ReplaceHeaderVal[hdr,"xyztUnits",XyztUnits["xyztUnits" /. hdr]];
	 hdr = ReplaceHeaderRule[hdr,"dataType",dataTypeNii];
	 hdr = ReplaceHeaderRule[hdr,"qformCode",coordinateNii];
	 hdr = ReplaceHeaderRule[hdr,"sformCode",coordinateNii];
	 
	 Close[strm];
	 
	 {"Header" -> hdr, "ByteOrder" -> byteorder}
 ]


(* ::Subsection:: *)
(*Various Special Nii import functions*)


(* ::Subsubsection::Closed:: *)
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
	{bvec, grad}=ImportBvalvec[RemoveExtention[file],FlipBvec->OptionValue[FlipBvec]];
	{data,Round[If[OptionValue[RotateGradients],grad.Inverse[mat], grad],.0001],bvec,vox}
]

ImportNiiDiff[fnii_String,fvec_String,fval_String,OptionsPattern[]]:=Module[{data,grad,bvec,vox,hdr,mat},
	{data,vox,hdr,mat}=ImportNii[fnii,NiiMethod -> "headerMat"];
	{bvec, grad} = ImportBvalvec[fval, fvec,FlipBvec->OptionValue[FlipBvec]];
	{data,Round[If[OptionValue[RotateGradients],grad.Inverse[mat], grad],.0001],bvec,vox}
]


(* ::Subsubsection::Closed:: *)
(*ImportNiiDix*)


ImportNiiDix[]:=Module[{Dix, vox, scale, B0, real, imag, mag, phase}, 
	{Dix, vox, scale} = ImportNii[NiiScaling -> False, NiiMethod -> "scaling"];
	{Dix, B0, {real, imag, mag, phase}} = CorrectDixonData[Dix, scale];
	{Dix, B0, {real, imag, mag, phase}, vox}
]

ImportNiiDix[file_String]:=ImportNiiDix[file,False]

ImportNiiDix[file_String, new_]:=Module[{Dix, vox, scale, B0, real, imag, mag, phase},  
	{Dix, vox, scale} = ImportNii[file, NiiScaling -> False, NiiMethod -> "scaling"];
	{Dix, B0, {real, imag, mag, phase}} =If[new, CorrectDixonDataNew[Dix, scale], CorrectDixonData[Dix, scale]];
	{Dix, B0, {real, imag, mag, phase}, vox}
]


CorrectDixonData[data_, scale_] := Block[{data0, B0, echos, phase, mag, real ,imag},
  (*fat,inphase,outphase,water*)
  data0 = data[[All, 1 ;; 4]];
  B0 = (scale[[1]] (data[[All, -1]] + 0.5) + scale[[2]]) /. (scale[[2]] + 0.5 scale[[1]]) -> 0.;
  (*I,M,P,R*)
  echos = data[[All, 5 ;; -2]];
  echos = Partition[Flatten[Transpose[echos], 1], 4*Length[data]];
  echos = Partition[#, 4] & /@ echos;
  (*convert I,P and R to radians*)
  phase = N[(2 Pi echos[[3]] / 4094.) - Pi] /. N[-Pi] -> 0.;
  mag = echos[[2]];
  {real, imag} = (echos[[{4, 1}]] - 2047.) / 4094.;
  (*Table[echos[[i]] = ((2 Pi echos[[i]]/4094.) - Pi) /. N[-Pi] -> 0., {i, {1, 3, 4}}];*)
  (*output*)
  {data0, B0, {real, imag, mag, phase}}
]


CorrectDixonDataNew[data_, scale_] := Block[{data0, B0, echos, phase, mag, real ,imag},
 (*fat,inphase,outphase,water*)
 data0 = data[[All, -5 ;; -2]][[All, {4, 2, 3, 1}]];
 B0 = (scale[[1]] (data[[All, -1]] + 0.5) + scale[[2]]) /. (scale[[2]] + 0.5 scale[[1]]) -> 0.;
 
 (*I,M,P,R*)
 echos = data[[All, ;; -6]];
 echos = Partition[Flatten[Transpose[echos], 1], 4*Length[data]];
 echos = Partition[#, 4] & /@ echos;
 (*convert I,P and R to radians*)
 phase = N[(2 Pi echos[[4]] / 4094.) - Pi] /. N[-Pi] -> 0.;
 mag = echos[[1]];
 {real, imag} = (echos[[{2, 3}]] - 2047.) / 4094.;
 
 (*output*)
 {data0, B0, {real, imag, mag, phase}}
 ]


(* ::Subsubsection::Closed:: *)
(*ImportNiiT2*)


ImportNiiT2[]:=Module[{T2, T2vox, T2cor, fit},
	{T2, T2vox} = ImportNii[NiiScaling -> False];
	{T2cor, fit} = CorrectMapData[T2,1];
	{T2cor, fit, T2vox}
]

ImportNiiT2[file_]:=Module[{T2, T2vox, T2cor, fit},
	{T2, T2vox} = ImportNii[file,NiiScaling -> False];
	{T2cor, fit} = CorrectMapData[T2,1];
	{T2cor, fit, T2vox}
]


(* ::Subsubsection::Closed:: *)
(*ImportNiiT1*)


ImportNiiT1[] := Module[{T1, T1vox, T1cor, fit},
  {T1, T1vox} = ImportNii[NiiScaling -> False];
  {T1cor, fit} = CorrectMapData[T1];
  {T1cor, fit, T1vox}]

ImportNiiT1[file_] := Module[{T1, T1vox, T1cor, fit},
  {T1, T1vox} = ImportNii[file,NiiScaling -> False];
  {T1cor, fit} = CorrectMapData[T1, 2];
  {T1cor, fit, T1vox}]


(* ::Subsubsection::Closed:: *)
(*CorrectMapData*)


CorrectMapData[datai_, maps_: 1] := 
 Block[{slices, echos, data, map, mask},
  {slices, echos} = Dimensions[datai][[1 ;; 2]] - {0, maps};
  (*Flatten the data and remove the T2map*)
  data = Flatten[Transpose[datai], 1];
  (**)
  map = data[[-maps slices ;;]];
  map = Partition[map, maps];
  (*partition to correct number of slices*)
  data = Partition[Drop[data, -maps slices], echos];
  data = Clip[NormalizeData[data], {0, Infinity}];
  
  {data, map}
  ]


(* ::Subsubsection::Closed:: *)
(*ImportExploreDTItens*)


ImportExploreDTItens[fil_String] := Module[{tens, vox, file},
  file = If[fil == "", FileSelect["FileOpen", {"*.nii"}, WindowTitle -> "Select the nii tensor file to import"], fil];
  If[file == Null, Return[]];
  {tens, vox} = ImportNii[file];
  tens = ((tens // Transpose)[[{4, 1, 6, 2, 5, 3}]]);
  {tens, vox}
  ]


(* ::Subsection:: *)
(*Importing val, vec and mat*)


(* ::Subsubsection::Closed:: *)
(*LineToList*)


LineToList[file_] := Module[{grads,tmp},
  grads = Import[file, "Lines"];
  grads = Transpose[(
  	tmp = StringReplace[#, {" " -> ",", "\t" -> ",", "E" -> "*10^", "e" -> "*10^"}];
	tmp = If[StringTake[tmp, -1] == ",", StringDrop[tmp, -1], tmp];
	tmp = ToExpression["{" <> tmp <> "}"]
  	) & /@ grads]
  ]


(* ::Subsubsection::Closed:: *)
(*ImportBval*)


SyntaxInformation[ImportBval] = {"ArgumentsPattern" -> {_.}};

ImportBval[]:=ImportBval[FileSelect["FileOpen", {".bval"}, WindowTitle -> "Select *.bval"]]

ImportBval[file_?StringQ]:=Flatten[LineToList[file] ]//N


(* ::Subsubsection::Closed:: *)
(*ImportBvec*)


Options[ImportBvec]={FlipBvec->False};

SyntaxInformation[ImportBvec] = {"ArgumentsPattern" -> {_.,OptionsPattern[]}};

ImportBvec[opts:OptionsPattern[]]:=ImportBvec[FileSelect["FileOpen", {".bvec"}, WindowTitle -> "Select *.bvec"], opts]

ImportBvec[file_?StringQ, OptionsPattern[]]:=Module[{grads},
	grads=Round[LineToList[file],0.0001];
	grads=If[OptionValue[FlipBvec],
		{1, -1, 1}#&/@grads,
		{1, -1, 1}RotateLeft[#]&/@grads
		];
	If[Negative[#[[3]]],-#,#]&/@grads
]


(* ::Subsubsection::Closed:: *)
(*ImportBvalvec*)


Options[ImportBvalvec]={FlipBvec->False};

SyntaxInformation[ImportBvalvec] = {"ArgumentsPattern" -> {_.,_.,OptionsPattern[]}};

ImportBvalvec[file_?StringQ,opts:OptionsPattern[]] := Module[{valf, vecf},
  valf = RemoveExtention[file]<>".bval";
  vecf = RemoveExtention[file]<>".bvec";
  ImportBvalvec[valf,vecf,opts]
]

ImportBvalvec[valf__?StringQ,vecf__?StringQ,opts:OptionsPattern[]] := {ImportBval[valf],ImportBvec[vecf,opts]}

ImportBvalvec[___,opts:OptionsPattern[]] := Module[{valf, vecf},
  valf = FileSelect["FileOpen", {"*.bval"}, WindowTitle -> "Select *.bval"];
  vecf = FileSelect["FileOpen", {"*.bvec"}, WindowTitle -> "Select *.bvec"];
  ImportBvalvec[valf,vecf,opts]
]



(* ::Subsubsection::Closed:: *)
(*ImportBmat*)


SyntaxInformation[ImportBmat] = {"ArgumentsPattern" -> {_.,_.}};

ImportBmat[] := ImportBmat[FileSelect["FileOpen", {"*.txt"}, WindowTitle -> "Select *.txt"]]

ImportBmat[fil_String] := Module[{bmati},
  If[fil == Null, Return[]];
  bmati = DeleteCases[#, ""] & /@ Import[fil, "Data"];
  Append[{-1, -1, -1, -1, -1, -1} #, 1] & /@ bmati[[All, {4, 1, 6, 2, 5, 3}]]
  ]


(* ::Subsection:: *)
(*ExportNii*)


(* ::Subsubsection::Closed:: *)
(*Register Export*)


ImportExport`RegisterExport["Nii",
	ExportNiiDefault,
	"DefaultElement" -> "Data",
	"AvailableElements" -> {
		"Data", "Header", "VoxelSize"
		},
	"OriginalChannel" -> True,
	"Options" -> {
		NiiDataType,
		NiiVersion
    	}
  ];


(* ::Subsubsection::Closed:: *)
(*ExportNii*)


SyntaxInformation[ExportNii] = {"ArgumentsPattern" -> {_,_,_., OptionsPattern[]}};

Options[ExportNii]={NiiDataType->Automatic,CompressNii->True}

ExportNii[dato_, voxi_, opts:OptionsPattern[]] := ExportNii[dato, voxi, "" ,opts]

ExportNii[dato_, voxi_, fil_, OptionsPattern[]] := Block[{fileo,data,type},
	
	fileo = If[fil == "", FileSelect["FileSave",{"*.nii"},"nifti",WindowTitle->"Select the destination file"], fil];
	If[fileo == Null || fileo === $Canceled || fileo === $Failed, Return[$Failed,Module]];
	
	(*if numbertyp is integer, Round data*)
	type=OptionValue[NiiDataType];
	data = ToPackedArray@Switch[type,"Integer",Round[dato],_,N[dato]];
	(*for lagecy reasons still allow Integer and Real*)
	type = type/.{"Integer"->"Integer16","Real"->"Real32"};
	
	(*compress the file*)
	If[OptionValue[CompressNii],
		fileo=fileo<>".gz";
		Export[fileo, {data, voxi}, {"GZIP", "Nii", {"Data", "VoxelSize"}}, NiiDataType->type],
		Export[fileo, {data, voxi}, {"Nii", {"Data", "VoxelSize"}}, NiiDataType->type]
		]; 
]


(* ::Subsubsection::Closed:: *)
(*ExportNiiDefault*)


Options[ExportNiiDefault] = {NiiDataType -> Automatic, NiiVersion -> 1, "Channel" -> Null, "ExtensionParsing" -> False}

ExportNiiDefault[file_, rule_, opts : OptionsPattern[]] := 
 Module[{ver, data, header, type, strm},
  
  ver = Clip[OptionsPattern[NiiVersion], {1, 2}];
  ver = If[NumberQ[ver], ver, 1];
  (*make the nii header*)
  header = MakeNiiHeader[rule, ver, opts];
  If[header === $Failed, Return[$Failed, Module]];
  {header, type} = header;
  (*get the data*)
  data = "Data" /. rule;
  (*write to file*)
  strm = OpenWrite[file, BinaryFormat -> True];
  BinaryWrite[strm, #[[1]], #[[2]]] & /@ header;
  BinaryWrite[strm, ArrangeData[data], type];
  Close[strm];
  ]


(* ::Subsubsection::Closed:: *)
(*MakeNiiHeader*)


MakeNiiHeader[rule_, ver_, OptionsPattern[ExportNiiDefault]] := Module[
  {vox, dim, ndim, type, range, data, header,
   headerInp, voxInp, headerDef, xoff, yoff, zoff},
  type = OptionValue[NiiDataType];
  
  (*get the data*)
  data = "Data" /. rule;
  (*check if data if number array*)
  If[! ArrayQ[data, _, NumberQ], Message[Export::niidat]; 
   Return[$Failed, Module]];
  (*get data properties*)
  ndim = ArrayDepth[data];
  dim = Dimensions[data];
   
  type = type /. (Automatic :> DetectDataType[data]);
  range = type /. rangeNii;
  
  (*Check data type and range*)
  If[! ArrayQ[data, _, type /. typeCheckNii], If[! ArrayQ[IntegerChop[data], _, type /. typeCheckNii], Message[Export::niitype, type]; Return[$Failed, Module]]];
  If[ListQ[range] && (Min[data] < range[[1]] || Max[data] > range[[2]]), Message[Export::niiran, MinMax[data]]; Return[$Failed, Module]];
  If[type == {"Byte", "Byte", "Byte"} && Last[dim] != 3, Message[Export::niidim, type]; Return[$Failed, Module]];
  
  (*is header given as input and is header a list of rules, 
  if given an valid use header*)
  header = "Header" /. rule;
  headerInp = MatchQ[header, {_Rule ..}];
  
  (*check if valid header*)
  If[headerInp,
   headerInp = 
    AllTrue[{Length[header] === 44 | Length[header] === 38}];
   If[! headerInp, Message[Export::niihdr]; Return[$Failed, Module]]
   ];
  
  (*is vox given as input and is header a list 3 number, 
  vox is given use vox (override in header) else use default 1x1x1*)
  vox = "VoxelSize" /. rule;
  voxInp = MatchQ[vox, {_?NumberQ, _?NumberQ, _?NumberQ}];
  (*if no voxel is given default*)
  vox = vox /. "VoxelSize" -> {1., 1., 1.};
  
  headerDef = {
    "size" -> Switch[ver, 1, 348, 2, 540],
    "data_Type" -> StringPadRight["", 10],
    "dbName" -> StringPadRight["", 18],
    "extents" -> 0,
    "sessionError" -> 0,
    "regular" -> "r",
    "dimInfo" -> DimInfo[{"x", "y", "z", Undefined}],(*input*)
    
    "dim" -> PadRight[Flatten[{ndim, If[ndim == 4, dim[[{4, 3, 1, 2}]], Reverse[dim]]}], 8, 1],(*input*)
    "intentP1" -> 0.,
    "intentP2" -> 0.,
    "intentP3" -> 0.,
    "intentCode" -> 0,
    "dataType" -> type /. Reverse[dataTypeNii, 2],(*input*)
    "bitPix" -> type /. Reverse[dataTypeNii, 2] /. typeSizeNii,(*input*)
    "sliceStart" -> 0,
    "pixDim" -> PadRight[Flatten[{1., Reverse[vox]}], 8, 0],(*input*)
    "voxOffset" -> Switch[ver, 1, 352, 2, 544],
    "scaleSlope" -> 1.,
    "scaleInteger" -> 0.,
    "sliceEnd" -> 0,
    "sliceCode" -> 0,
    "xyztUnits" -> XyztUnits[{"Millimeters", "Seconds"}],(*input*)
    "calMax" -> 0.`,
    "calMin" -> 0.`,
    "sliecDuration" -> 0.`,
    "tOffset" -> 0.`,
    "glMax" -> 0,
    "glMin" -> 0,
    
    "descrip" -> StringPadRight["Created with QMRITools", 80, FromCharacterCode[0]],(*input*)
    "auxFile" -> StringPadRight["None", 24, FromCharacterCode[0]],
    "qformCode" -> "Coregistration" /. Reverse[coordinateNii, 2],
    "sformCode" -> "Scanner Posistion" /. Reverse[coordinateNii, 2],
    "quaternB" -> 0,
    "quaternC" -> 0,
    "quaternD" -> 0,
    "qOffsetX" -> (xoff = -N[vox[[3]] dim[[-1]]/2]),
    "qOffsetY" -> (yoff = -N[vox[[2]] dim[[-2]]/2]),
    "qOffsetZ" -> (zoff = -N[vox[[1]] dim[[1]]/2]),
    "sRowx" -> {vox[[3]], 0., 0., xoff},
    "sRowy" -> {0., vox[[2]], 0., yoff},
    "sRowz" -> {0., 0., vox[[1]], zoff},
    "intentName" -> StringPadRight["", 16],
    "magic" -> StringPadRight[Switch[ver, 1, "n+1", 2, "n+2"], 4, FromCharacterCode[0]],
    "ECode" -> {0, 0, 0, 0}
    };
  
  header = GetNiiHeaderValues[ver] /. headerDef;
  
  {header, type}
  ]


(* ::Subsubsection::Closed:: *)
(*DetectDataType*)


DetectDataType[input_] := Module[{data, min, max},
  data = Flatten[input];
  Switch[Head[Total[data]],
   Integer,
   {min, max} = MinMax[data];
   If[min >= 0 && max <= 18446744073709551615,
   (*positive integer*)
   Switch[Total@Boole[# >= max & /@ {255, 65535, 4294967295, 18446744073709551615}], 
    	4, "Byte", 
    	3, "UnsignedInteger16", 
    	2, "UnsignedInteger32", 
    	1, "UnsignedInteger64"
    ],
    If[max <= 9223372036854775807 && min >= -9223372036854775807,
    (*signed integer*)
    Switch[Total@Boole[(max <= # && min >= -#) & /@ {127, 32767, 2147483647, 9223372036854775807}], 
     	4, "Integer8", 
     	3, "Integer16", 
     	2, "Integer32", 
     	1, "Integer64"
     ],
     (*else real*)
     "Real32"
     ]
    ],
   Real, "Real32",
   Rational, "Real32",
   Complex, "Complex64",
   _, "Real32"
   ]
  ]


(* ::Subsection:: *)
(*Exporting val, vec and mat*)


(* ::Subsubsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
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
(*Compress and Extract*)


(* ::Subsubsection::Closed:: *)
(*ExtractNiiFiles*)


SyntaxInformation[ExtractNiiFiles] = {"ArgumentsPattern" -> {_.,_.}};

ExtractNiiFiles[lim_:Infinity] := ExtractNiiFiles[FileSelect["Directory", WindowTitle -> "Select direcotry containig the nii files"],lim]

ExtractNiiFiles[folder_,lim_:Infinity] := Module[{files},
  files = FileNames["*.nii.gz", folder,lim];
  ExtractNiiFile/@ files;
]

ExtractNiiFile[file_] := Quiet[ExtractArchive[file, StringDrop[DirectoryName[file], -1]];DeleteFile[file];]


(* ::Subsubsection::Closed:: *)
(*CompressNiiFiles*)


SyntaxInformation[CompressNiiFiles] = {"ArgumentsPattern" -> {_.,_.}};

CompressNiiFiles[]:=CompressNiiFiles[Infinity]

CompressNiiFiles[folder_?StringQ]:=CompressNiiFiles[folder, Infinity]

CompressNiiFiles[lim_] := CompressNiiFiles[FileSelect["Directory", WindowTitle -> "Select direcotry containig the nii files"],lim]

CompressNiiFiles[folder_?StringQ,lim_] := Module[{files,file},
   files = FileNames["*.nii", folder,lim];
   PrintTemporary[Dynamic[file]];
   (file=#;CompressNiiFile[#])&/@files;
   ]

CompressNiiFile[file_] := Quiet[DeleteFile[file <> ".gz"];CreateArchive[file, file <> ".gz"];DeleteFile[file];]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
