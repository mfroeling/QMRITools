(* ::Package:: *)

(* ::Title:: *)
(*QMRITools GeneralTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`GeneralTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`GeneralTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


GetAssetLocation::usage = 
"GetAssetLocation[name] Gets the location of the executable assets of the package for the highest installed version.
Current assess are \"Elastix\", \"Transformix\" and \"DcmToNii\"." 

ExtractDemoData::usage = 
"ExtractDemoData[] Extracts the demo data archive."

OpenDemonstrationNotebook::usage =
"OpenDemonstrationNotebook[] Opens the demonstration notebook."

OpenQMRIToolsLocation::usage =
"OpenQMRIToolsLocation[] Opens the QMRITools location in the file explorer."

SetDemoDirectory::usage =
"SetDemoDirectory[] Sets the directory to the demo data directory."


StringPadInteger::usage = 
"StringPadInteger[num] converts the integer num to a string and pads it with zeros to length 3.
StringPadInteger[{num, len}] converts the integer num to a string and pads it with zeros to length len.
StringPadInteger[pre, num] the same but with prefix pre.
StringPadInteger[pre, {num, len}] the same but with prefix pre.
StringPadInteger[num, post] the same but with postfix post.
StringPadInteger[{num, len}, post] the same but with postfix post.
StringPadInteger[pre, num, post] the same but with pre and postfix pre and post.
StringPadInteger[post, {num, len}, post] the same but with pre and postfix pre and post."

DateName::usage =
"DateName[] gives the current date and time in the format \"YYYYMMDD-HH\"."

FileSelect::usage = 
"FileSelect[action] creates a system dialog which returns file/foldername action can be \"FileOpen\", \"FileSave\" or \"Directory\".
FileSelect[action, {type}] same but allows the definition of filetypes for \"FileOpen\" and \"FileSave\" e.g. \"jpg\" or \"pdf\"."

ConvertExtension::usage = 
"ConvertExtension[filename, extension] checks if file has correct extension. Removes .gz or changes the extension or adds extension if not present."

EmptyDirectoryQ::usage = 
"EmptyDirectoryQ[dir] checks if directory dir is empty."

NiiFileExistQ::usage = 
"NiiFileExistQ[file] checks if the *.nii or *.nii.gz file exists."


SaveImage::usage = 
"SaveImage[image] exports graph to image, ImageSize, FileType and ImageResolution can be given as options.
SaveImage[image, \"filename\"] exports graph to image with \"filename\", ImageSize, FileType and ImageResolution can be given as options."


GridData::usage = 
"GridData[{data1,data2,...}, part] makes a grid of multiple datasets with part sets on each row."

GridData3D::usage = 
"GridData3D[{data1,data2,...}, part] same as grid data, but only works on 4D data where the data is gridded in axial, coronal and sagittal."


RotateDimensionsLeft::usage = 
"RotateDimensionsLeft[data] rotates the dimensions of the data one to the left.
RotateDimensionsLeft[data, i] rotates the dimensions of the data i to the left."

RotateDimensionsRight::usage = 
"RotateDimensionsRight[data] rotates the dimensions of the data one to the right.
RotateDimensionsRight[data, i] rotates the dimensions of the data i to the right."

ReverseDimensions::usage = 
"ReverseDimensions[data] reverses the dimensions of the data."


FindMaxDimensions::usage = 
"FindMaxDimensions[{data1, data2, ..}] finds the maximal dimensions of all datasets. Each dataset is 3D."

PadToDimensions::usage = 
"PadToDimensions[data] pads the data to the max dimensions of data, using FindMaxDimensions.
PadToDimensions[data, dim] pads the data to dimensions dim." 

RescaleData::usage = 
"RescaleData[data,dim] rescales image/data to given dimensions.
RescaleData[data,{vox1, vox2}] rescales image/data from size vox1 to size vox2."


DataToVector::usage = 
"DataToVector[data] converts the non zero data to vector.
DataToVector[data, mask] converts the data within the mask to vector.

the data can be reconstructed using VectorToData.

output is the vectorized data and a list containing the original data dimensions and a list with the data coordinates. {vec, {dim,pos}}."

VectorToData::usage = 
"VectorToData[vec, {dim,pos}] converts the vectorized data from DataToVector back to its original Dimensions."

MakeCoordinates::usage = 
"MakeCoordinates[data, vox] gives the coordinates of every voxel.
MakeCoordinates[dim, vox] gives the coordinates of every voxel for a dataset with dimensions dim."

Squeeze::usage =
"Squeeze[data] Removes the singleton dimensions from data."


TensMat::usage=
"TensMat[tensor] transforms tensor form vector format {xx,yy,zz,xy,xz,yz} to matrix format {{xx,xy,xz},{xy,yy,yz},{xz,yz,zz}}."

TensVec::usage=
"TensVec[tensor] transforms tensor form matrix format {{xx,xy,xz},{xy,yy,yz},{xz,yz,zz}} to vector format {xx,yy,zz,xy,xz,yz}."


CropData::usage =
"CropData[data] creates a dialog window to crop the data (assumes voxel size (1,1,1)).
CropData[data,vox] creates a dialog window to crop the data."

FindCrop::usage = 
"FindCrop[data] finds the crop values of the data by removing all zeros surrounding the data."

AutoCropData::usage = 
"AutoCropData[data] crops the data by removing all background zeros."

ReverseCrop::usage = 
"ReverseCrop[data, dim, crop] reverses the crop on the cropped data with crop values crop to the original size dim.
ReverseCrop[data, dim, crop, {voxorig, voxnew}] reverses the crop on the cropped data with crop values crop to the original size dim."

ApplyCrop::usage =
"ApplyCrop[data,crop] applies the cropped region obtained form CropData to the data.
ApplyCrop[data,crop,{voxorig,voxnew}] applies the cropped region obtained form CropData to the data." 


CutData::usage = 
"CutData[data] splits the data in two equal sets left and right."

StichData::usage =
"StichData[datal,datar] joins left and right part of the data generated by CutData."


DivideNoZero::usage = 
"DivideNoZero[a, b] divides a/b but when b=0 the result is 0. a can be a number or vector."

MeanNoZero::usage = 
"MeanNoZero[data] calculates the mean of the data ignoring the zeros."

StandardDeviationNoZero::usage = 
"StandardDeviationNoZero[data] calculates the mean of the data ignoring the zeros."

MedianNoZero::usage = 
"MedianNoZero[data] calculates the Median of the data ignoring the zeros."

LogNoZero::usage = 
"LogNoZero[val] return the log of the val which can be any dimension array. if val=0 the output is 0."

ExpNoZero::usage = 
"ExpNoZero[val] return the Exp of the val which can be any dimension array. if val=0 the output is 0."

RMSNoZero::usage = 
"RMSNoZero[vec] return the RMS error of the vec which can be any dimension array. if vec={0...} the output is 0. Zeros are ignored."

MADNoZero::usage = 
"MADNoZero[vec] return the MAD error of the vec which can be any dimension array. if vec={0...} the output is 0. Zeros are ignored."

SignNoZero::usage = 
"SignNoZero[val] gives the sign of the val, where the sign of val > 0 is 1 and val < 0 is -1." 

SumOfSquares::usage = 
"SumOfSquares[{data1, data2, .... datan}] calculates the sum of squares of the datasets.
Output is the SoS and the weights, or just the SoS."

LLeastSquares::usage = 
"LLeastSquares[A, y] = performs a Linear Linear Least Squares fit.
It uses a compiled version of the Pseudo inverse of A."

NNLeastSquares::usage = 
"NNLeastSquares[A, y] performs a Non Negative Linear Least Squares fit.
finds an x that solves the linear least-squares problem for the matrix equation A.x==y.

output is the solution x."

BSplineCurveFit::usage = 
"BSplineCurveFit[points] fits a b-spline to the points. Output is a list of same size as points."


LapFilter::usage = 
"LapFilter[data] Laplacian filter of data with kernel size 0.8.
LapFilter[data, ker] Laplacian filter of data with kernel ker."

MedFilter::usage = 
"MedFilter[data] Median filter of data with kernel size 1.
MedFilter[data, ker] Median filter of data with kernel ker."

StdFilter::usage =
"StdFilter[data] StandardDeviation filter of data using gaussian kernel 2. 
StdFilter[data, ker] StandardDeviation filter of data using kernel with size ker."


GyromagneticRatio::usage =
"GyromagneticRatio[] gives the gyromagnetic ratio for \"1H\" in MHz/T.
GyromagneticRatio[nuclei] gives the gyromagnetic ratio for the nuclei, e.g. \"31P\" of \"1H\"."


DynamicPartition::usage = 
"DynamicPartition[data, {part}] partitions the data into parts which is a list of integers. The remainders is los. 
DynamicPartition[data,part,last] partitions the data into parts which is a list of integers. The remainders is partitioned into equal parts defined by last.
If last is All, the remainders is just one partition."


MakeIntFunction::usage = 
"MakeIntFunction[data] makes an interpolation function of the data using voxel size {1, 1, 1} and interpolation order 1
MakeIntFunction[data, int] makes an interpolation function of the data using voxel size {1, 1, 1} and interpolation order int.
MakeIntFunction[data, vox ,int] makes an interpolation function of the data using voxel size vox and interpolation order int."


DecomposeScaleMatrix::usage = 
"DecomposeScaleMatrix[mat] decomposes the affine matrix in T, R, S and Q."

DecomposeAffineMatrix::usage = 
"DecomposeAffineMatrix[S] decomposes the scale matrix in S1, S2 and S3."

QuaternionToRotationMatrix::usage =
"QuaternionToRotationMatrix[{a, b,c,d}] converts quaternion to rotation matrix R."

QuaternionVectorToRotationMatrix::usage =
"QuaternionVectorToRotationMatrix[{b,c,d}] converts quaternion to rotation matrix R."

RotationMatrixToQuaternion::usage = 
"RotationMatrixToQuaternion[R] converts rotation matrix to quaternions {a, b,c,d}."

RotationMatrixToQuaternionVector::usage =
"RotationMatrixToQuaternionVector[R] converts rotation matrix to quaternions {b,c,d}."


QMRIToolsPackages::usage =
"QMRIToolsPackages[] give list of all the QMRITools packages."

QMRIToolsFunctions::usage = 
"QMRIToolsFunctions[] give list of all the QMRITools packages, functions and options.
QMRIToolsFunctions[p] print a table with length p of all the QMRITools functions and options.
QMRIToolsFunctions[\"toolbox\"] gives a list of all the functions and options in toolbox.
QMRIToolsFunctions[\"toolbox\", p] gives a table op length p of all the functions and options in toolbox. If toolbox is \"All\" it will list all toolboxes."

QMRIToolsFuncPrint::usage = 
"QMRIToolsFuncPrint[] gives a list of all the QMRITools functions with their usage information."

CompilableFunctions::usage = 
"CompilableFunctions[] generates a formatted table of all compilable functions generated by Compile`CompilerFunctions."


MakeFunctionGraph::usage = 
"MakeFunctionGraph[function] makes a function dependency graph of the function."


MemoryUsage::usage = 
"MemoryUsage[] gives a table of which definitions use up memory.
MemoryUsage[n] gives a table of which definitions use up memory, where n is the amount of definitions to show."

MBCount::usage =
"MBCount[expr] gives the memory usage of the expression in MB."

ClearTemporaryVariables::usage = 
"ClearTemporaryVariables[] Clear temporary variables."


(* ::Subsection::Closed:: *)
(*General Options*)


PadValue::usage = 
"PadValue is an option for PadToDimensions. It specifies the value of the padding."

PadDirection::usage = 
"PadDirection is an option for PadToDimensions. It specifies the direction of padding, \"Center\", \"Left\" or \"Right\"."


CropOutput::usage = 
"CropOutput is an option for CropData, can be \"All\",\"Data\" or \"Crop\"."

CropInit::usage = 
"CropInit is an option for CropData. By default the crop is not initialized bu can be with {{xmin,xmax},{ymin,ymax},{zmin,zmax}}."

CropPadding::usage = 
"CropPadding is an option for AutoCropData or FindCrop. It specifies how much padding to use around the data."

CropAlways::usage = 
"CropAlways is an optin for ApplyCrop. If set True is will always crop even if outside the data."

OutputWeights::usage = 
"OutputWeights is an option for SumOfSquares. If True it also output the SoS weights."


SplineKnotsNumber::usage =
"SplineKnotsNumber is an option for BSplineCurveFit and defines how many knots the b-spline has."

SplineRegularization::usage =
"SplineRegularization is an option for BSplineCurveFit and defines the amount of regularization for the linear fit."


CenterVoxel::usage = 
"CenterVoxel is an option for MakeIntFunction. If set True the centers of the voxels are interpolated else its the corners."

CenterRange::usage = 
"CenterRange is an option for MakeIntFunction. If set True the centers of the dataset is the origin else its the corner."


MonitorCalc::usage = 
"MonitorCalc is an option for many processing functions. When true the process of the calculation is shown."


LabelPlacement::usage =
"LabelPlacement is an option for MakeFunctionGraph. Defines where to place the label of the function graph. Accepts values that can be used in Placed."

AllowSelfDependencies::usage = 
"AllowSelfDependencies is and option for MakeFunctionGraph. Can be True or False. If True a function that calls itself is also shown."


(* ::Subsection::Closed:: *)
(*Error Messages*)


RescaleData::dim = "Given dimensions `1` not the same depth as that of the given data `2`."

RescaleData::data = "Error: Input must be 2D with {xdim,ydim} input or 3D dataset with {xdim,ydim} or {zdim, xdim, ydim} input."


DataToVector::dim = "`1` should be 2D, 3D or 4D, data is `2`D.";

DataToVector::mask = "Data and mask should have the same dimensions: data `1` and mask `2`";


ApplyCrop::dim = "Crop region lies outside data range."


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*Asset Functions*)


GetAssetLocation[name_] := Block[{file},
	file = Last[SortBy[PacletFind["QMRITools"], #["Version"]&]]["AssetLocation", name];
	If[FileExistsQ[file],file,$Failed]
]


ExtractDemoData[] := Block[{file},
	file = GetAssetLocation["DemoData"];
	If[! DirectoryQ[FileNameJoin[{DirectoryName[file], "DemoData"}]],
		If[FileExistsQ[file], Quiet@ExtractArchive[file, DirectoryName[file]],
		Print["DemoData archive does not exist"]]
	];
]


OpenDemonstrationNotebook[] := NotebookOpen[GetAssetLocation["Demo"]];


OpenQMRIToolsLocation[] := SystemOpen[First[PacletFind["QMRITools"]]["Location"]];


SetDemoDirectory[]:=SetDirectory[FileNameJoin[{DirectoryName[GetAssetLocation["DemoData"]], "DemoData"}]];


(* ::Subsection:: *)
(*General Functions*)


(* ::Subsubsection::Closed:: *)
(*StringPad*)


StringPadInteger[x_?IntegerQ]:=StringPadInteger["", {x, 3}, ""]

StringPadInteger[{x_?IntegerQ, n_?IntegerQ}]:=StringPadInteger["", {x, n}, ""]

StringPadInteger[pre_?StringQ, x_?IntegerQ]:=StringPadInteger[pre, {x, 3}, ""]

StringPadInteger[pre_?StringQ, {x_?IntegerQ, n_?IntegerQ}]:=StringPadInteger[pre, {x, n}, ""]

StringPadInteger[x_?IntegerQ, post_?StringQ]:=StringPadInteger["", {x, 3}, post]

StringPadInteger[{x_?IntegerQ, n_?IntegerQ}, post_?StringQ]:=StringPadInteger["", {x, n}, post]

StringPadInteger[pre_?StringQ, {x_?IntegerQ, n_?IntegerQ}, post_?StringQ]:=pre<>StringPadLeft[ToString[x], n, "0"]<>post


(* ::Subsubsection::Closed:: *)
(*DateName*)


DateName[]:=StringReplace[DateString[{"YearShort", "Month", "Day", "-", "Hour"}], ":" -> ""]


(* ::Subsubsection::Closed:: *)
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
		Print["Canceled!"];$Canceled
		,
		If[action == "Directory",
			StringDrop[input,-1],
			input
		]
	]
]


(* ::Subsubsection::Closed:: *)
(*CheckExtension*)


SyntaxInformation[ConvertExtension] = {"ArgumentsPattern" -> {_, _}};


ConvertExtension[fileIn : {_?StringQ ..}, ext_?StringQ] := ConvertExtension[#, ext] & /@ fileIn

ConvertExtension[fileIn_?StringQ, ext_?StringQ] := Block[{extOld, file, extNew},
	(*get filename and extension*)
	file = StringReplace[fileIn,".gz"->""];
	extOld = "." <> FileExtension[file];

	(*make new extension*)
	extNew = If[StringTake[ext, 1] === ".", ext, "." <> ext];

	(*add or replace *)
	If[extOld === ".",
		file <> extNew,
		StringReplace[file, extOld -> extNew]
	]
]


(* ::Subsubsection::Closed:: *)
(*EmptyDirectoryQ*)


SyntaxInformation[EmptyDirectoryQ] = {"ArgumentsPattern" -> {_}};

EmptyDirectoryQ[dir_] := FileNames[All, dir] === {}


(* ::Subsubsection::Closed:: *)
(*NiiFileExistQ*)


NiiFileExistQ[file_] := FileExistsQ[ConvertExtension[file, ".nii"]] || FileExistsQ[ConvertExtension[file, ".nii.gz"]]


(* ::Subsubsection::Closed:: *)
(*RotateDimensionsLeft*)


SyntaxInformation[RotateDimensionsLeft] = {"ArgumentsPattern" -> {_, _.}};

RotateDimensionsLeft[dat_, i_ : 1] := RotateDimensions[dat, i, RotateRight]


(* ::Subsubsection::Closed:: *)
(*RotateDimensionsLeft*)


SyntaxInformation[RotateDimensionsRight] = {"ArgumentsPattern" -> {_, _.}};

RotateDimensionsRight[dat_, i_ : 1] := RotateDimensions[dat, i, RotateLeft]


(* ::Subsubsection::Closed:: *)
(*RotateDimensionsLeft*)


SyntaxInformation[ReverseDimensions] = {"ArgumentsPattern" -> {_, _.}};

ReverseDimensions[dat_] := RotateDimensions[dat, 1, Reverse]


(* ::Subsubsection::Closed:: *)
(*RotateDimensions*)


RotateDimensions[dat_, 0, dir_] := dat

RotateDimensions[dat_, i_, dir_] := Transpose[ToPackedArray[dat], dir[Range[ArrayDepth[dat]], i]]


(* ::Subsubsection::Closed:: *)
(*SaveImage*)


Options[SaveImage] = {ImageSize -> 6000, FileType -> ".jpg", ImageResolution -> 300};

SyntaxInformation[SaveImage] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

SaveImage[exp_, opts : OptionsPattern[]] := Module[{input, type},
	type = OptionValue[FileType];
	type = If[StringTake[type, 2] === "*.", StringDrop[type,1], If[StringTake[type, 1] === ".",type, "."<>type]];

	input = FileSelect["FileSave", {"*"<>type}];

	If[input != "Canceled!",
		input=If[StringTake[input,-4]===type,input,input<>type];
		SaveImage[exp, input, opts]
	];
]

SaveImage[exp_, filei_String, OptionsPattern[]] := Module[{file,imsize,res,type},
	type = OptionValue[FileType];
	type = If[StringTake[type, 2] === "*.", StringDrop[type,1], If[StringTake[type, 1] === ".",type, "."<>type]];

	file=If[StringTake[filei,-4]===type||StringTake[filei,-5]===type,filei,filei<>type];

	imsize=OptionValue[ImageSize];
	res=OptionValue[ImageResolution];

	If[OptionValue[FileType]===".tiff"||OptionValue[FileType]===".tif",
		Export[file, exp , ImageSize->imsize, ImageResolution -> res,"ImageEncoding"->"LZW"],
		Export[file, exp , ImageSize->imsize, ImageResolution -> res]
	];

	Print["File was saved to: " <> file];
]


(* ::Subsection:: *)
(*Reshaping and Resizing*)


(* ::Subsubsection::Closed:: *)
(*PadToDimensions*)


SyntaxInformation[FindMaxDimensions] = {"ArgumentsPattern" -> {_}};

FindMaxDimensions[data_] := Max /@ Transpose[(Dimensions /@ data)]


(* ::Subsubsection::Closed:: *)
(*PadToDimensions*)


Options[PadToDimensions]={PadValue->0., PadDirection -> "Center"}

SyntaxInformation[PadToDimensions] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

PadToDimensions[data_, opts:OptionsPattern[]]:=PadToDimensions[data, FindMaxDimensions[data], opts]

PadToDimensions[data_, dim_, opts:OptionsPattern[]] := Block[{diffDim, padval, pad,dir,zer},
	If[Length[Dimensions[First[data]]]===Length[dim],
		(*pad a series of datasets*)		
		PadToDimensions[#, dim, opts]&/@data
		,

		(*Pad an individual dataset*)
		padval = OptionValue[PadValue];
		diffDim = dim - Dimensions[data];

		dir = OptionValue[PadDirection];
		dir = If[StringQ[dir], ConstantArray[dir, Length[dim]], dir];
		zer = ConstantArray[0, Length[dim]];

		pad = MapThread[
			Switch[#1, "Left", #2, "Right", #3, _, #4] &, {
				dir, 
				Transpose@{zer, diffDim}, Transpose@{diffDim, zer}, 
				Transpose@{Floor[diffDim/2], Ceiling[diffDim/2]}
			}
		];

		ToPackedArray[N@ArrayPad[data,pad,padval]]
	]
]


(* ::Subsubsection::Closed:: *)
(*RescaleData*)


Options[RescaleData] = {InterpolationOrder -> 3};

SyntaxInformation[RescaleData] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

RescaleData[data_?ArrayQ, {v1_?VectorQ, v2_?VectorQ}, opts : OptionsPattern[]] := RescaleDatai[data, v1/v2, "v", opts]

RescaleData[data_?ArrayQ, dim_?VectorQ, opts : OptionsPattern[]] := RescaleDatai[data, dim, "d", opts]

Options[RescaleDatai] = {InterpolationOrder -> 3};

RescaleDatai[data_?ArrayQ, sc_?VectorQ, met_, opts : OptionsPattern[]] := Block[{type, dim, int, dataOut},
	dim = Dimensions[data];
	int = OptionValue[InterpolationOrder];

	dataOut = Switch[ArrayDepth[data],
		(*rescale an image*)
		2,
		If[Length[sc] != 2,
			Return[Message[RescaleData::dim, sc, Dimensions[data]]];,
			RescaleImgi[data, {sc, met}, int]
		],
		3 (*rescale a 3D dataset*),
		Switch[Length[sc],
			2 (*rescale a stac of 2D images*),
			RescaleImgi[#, {sc, met}, int] & /@ data,
			3 (*rescale 3D data*),
			RescaleImgi[data, {sc, met}, int],
			_,
			Return[Message[RescaleData::dim, sc, Dimensions[data]]];
		],
		4 (*rescale a 4D dataset, treat data as multiple 3D sets*),
		Transpose[RescaleDatai[#, sc, met, opts] & /@ Transpose[data]],
		_,
		Return[Message[RescaleData::data]];
	];

	ToPackedArray[N@Chop[Clip[dataOut,MinMax[data]]]]
]


RescaleImgi[dat_, {sc_, met_}, n_] := Block[{type, im, dim},
	(*data type*)
	type = If[ArrayQ[dat, _, IntegerQ], "Bit16", "Real32"];
	dim = If[met == "v", Round[sc Dimensions[dat]], sc];
	(*convert to 2D or 3D image*)
	im = Switch[ArrayDepth[dat], 2, Image[dat, type], 3, Image3D[dat, type]];
	ImageData[ImageResize[im, Reverse[dim], Resampling ->{"Spline", n}], type]
]


(* ::Subsubsection::Closed:: *)
(*GridData*)


Options[GridData] = {Padding-> None}

SyntaxInformation[GridData] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

GridData[dati_, opts:OptionsPattern[]]:=GridData[dati, Ceiling[Sqrt[Length@dati]], opts]

GridData[dati_, part_, opts:OptionsPattern[]] := Block[{dim, data, adepth, pad, val},
	adepth = ArrayDepth[dati[[1]]];

	(*pad the images with zeros*)
	pad = OptionValue[Padding];
	data = If[pad =!= None,
		{pad, val} = If[IntegerQ[pad], {pad, 0.}, pad];
		pad = PadLeft[{pad, pad}, adepth];
		ArrayPad[#, pad, val] & /@ dati
		,
		ToPackedArray@N@dati
	];

	(*make the first dimension such that it is dividable by part*)	
	dim = Dimensions[data];
	dim[[1]] = dim[[1]] + (part - (Mod[dim[[1]], part] /. 0 -> part));
	data = Transpose[Partition[PadRight[data, dim], part]];

	(*make the grid*)
	data = MapThread[Join, #, adepth - 2] & /@ data;
	data = MapThread[Join, data, adepth - 1];
	ToPackedArray@N@data
]


(* ::Subsubsection::Closed:: *)
(*GridData3D*)


SyntaxInformation[GridData3D] = {"ArgumentsPattern" -> {_, _}};

GridData3D[data_, part_] := Block[{AX, COR, SAG},
	AX = GridData[data, part];
	COR = GridData[Transpose[Reverse@#, {2, 1, 3}] & /@ data, part];
	SAG = GridData[Transpose[Reverse@#, {2, 3, 1}] & /@ data, part];
	If[Dimensions[AX]===Dimensions[COR]===Dimensions[SAG],
		Transpose[{AX, COR, SAG}],
		{AX, COR, SAG}
	]
]


(* ::Subsubsection::Closed:: *)
(*DataToVector*)


SyntaxInformation[DataToVector] = {"ArgumentsPattern" -> {_, _.}};

DataToVector[datai_] := DataToVector[datai, 1]

DataToVector[datai_, maski_] := Module[{data, sp, mask, depthd, depthm, depth, dimm, dimd},
	data = ToPackedArray[N[datai]];
	dimd = Dimensions[data];
	depthd = ArrayDepth[data];

	(*data dimensions are not correct, mask must be 2D or 3D*)
	If[! (depthd == 2 || depthd == 3 || depthd == 4), 
	Return@Message[DataToVector::dim, "Data", depthd]];
	If[maski =!= 1,
		(*masks is given*)
		mask = maski;
		dimm = Dimensions[mask];
		depthm = ArrayDepth[mask];
		
		(*check mask*)
		depth = depthd - depthm;
		If[! (depthm == 2 || depthm == 3), Message[DataToVector::dim, "Mask", depthm]];
		dimd = If[depth == 0,
			(*mask and data are same dimensions*)
			dimd,
			(*data is one dimensions larger than mask either 2D and 3D or 3D and 4D*)
			If[depth == 1 && depthd == 3, Drop[dimd, 1], Drop[dimd, {2}]]
		];
		(*Dimensions must be equal*)
		If[dimd =!= dimm, Return@Message[DataToVector::mask, dimd, dimm]];
		,
		(*mask is not given make mask*)
		mask = If[depthd == 4, Unitize[Mean@Transpose@data], 1]
	];

	If[mask === 1,
		sp = SparseArray[data];
		{sp["ExplicitValues"], {dimd, sp["ExplicitPositions"]}}
		,
		(*Flatten the data*)
		data = If[depthd == 4, 
		Flatten[RotateDimensionsLeft[Transpose[data]], 2], 
		If[depthd == 3 && depth == 1, Flatten[RotateDimensionsLeft[data], 1], Flatten[data]]];
		(*get the data and positions there mask is 1*)
		{Pick[data, Round[Flatten[mask]], 1], {dimd, Position[mask, 1]}}
	]
]


(* ::Subsubsection::Closed:: *)
(*VectorToData*)


SyntaxInformation[VectorToData] = {"ArgumentsPattern" -> {_, {_, _}}};

VectorToData[vec_, {dim_, pos_}] := ToPackedArray@N@If[VectorQ[vec],
	Normal[SparseArray[pos -> vec, dim]],
	If[Length[dim] == 2,
		Normal[SparseArray[pos -> #, dim]] & /@ Transpose[vec],
		Transpose[Normal[SparseArray[pos -> #, dim]] & /@ Transpose[vec]]
	]
]


(* ::Subsubsection::Closed:: *)
(*FitGradientMap*)


SyntaxInformation[FitGradientMap] = {"ArgumentsPattern" -> {_, _., _.}};

FitGradientMap[data_]:=FitGradientMap[data, 2, 1]

FitGradientMap[data_, ord_]:=FitGradientMap[data, ord, 1]

FitGradientMap[data_, ord_, smp_]:=FitGradientMap[{data, 1}, ord, smp]

FitGradientMap[{data_, msk_}, ord_, smp_] := Block[{val, dim, coor, fit, x, y, z},
	Clear[x, y, z];
	{val, {dim, coor}} = DataToVector[msk data];
	fit = Fit[
		Thread[{coor[[;;;;smp, 1]], coor[[;;;;smp, 2]], coor[[;;;;smp, 3]], val[[;;;;smp]]}], 
		DeleteDuplicates[Times @@ # & /@ Tuples[{1, x, y, z}, ord]], 
		{x, y, z}
	];
	If[msk=!=1, coor = DataToVector[data][[2,2]]];
	{x, y, z} = Transpose[coor];
	VectorToData[fit, {dim, coor}]
]


(* ::Subsubsection::Closed:: *)
(*MakeCoordinates*)


SyntaxInformation[MakeCoordinates] = {"ArgumentsPattern" -> {_, _}};

MakeCoordinates[dim_?VectorQ,vox_]:=vox RotateDimensionsRight@Array[{##}&, dim]

MakeCoordinates[dat_?ArrayQ,vox_]:=MakeCoordinates[Dimensions@dat,vox]


(* ::Subsubsection::Closed:: *)
(*TensMat*)


SyntaxInformation[TensMat] = {"ArgumentsPattern" -> {_}};

TensMat[tens : {_?ArrayQ ..}] := RotateDimensionsLeft[TensMati[tens], 2];

TensMat[tens_?VectorQ] := TensMati[tens]


TensMati[{xx_, yy_, zz_, xy_, xz_, yz_}] := {{xx, xy, xz}, {xy, yy, yz}, {xz, yz, zz}};


(* ::Subsubsection::Closed:: *)
(*TensVec*)


SyntaxInformation[TensVec] = {"ArgumentsPattern" -> {_}};

TensVec[{{xx_, xy_, xz_}, {_, yy_, yz_}, {_, _, zz_}}] := {xx, yy, zz, xy, xz, yz};

TensVec[tens : {_?ArrayQ ..}] := TensVec[RotateDimensionsRight[tens, 2]]


(* ::Subsection:: *)
(*Cropping*)


(* ::Subsubsection::Closed:: *)
(*CropData*)


Options[CropData] = {CropOutput -> "All", CropInit -> Automatic};

SyntaxInformation[CropData] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

CropData[data_, opts:OptionsPattern[]] := CropData[data, {1,1,1}, opts]

CropData[data_, vox:{_?NumberQ, _?NumberQ, _?NumberQ}, OptionsPattern[]] := Block[
	{a, b, c, d, e, f, clipall, dd, dataout, output, init},

	NotebookClose[cropwindow];
	dd = ArrayDepth[data];

	DynamicModule[{dat, zd, xd, yd, outp, size,  r1, r2, r3},

		dat = Switch[dd, 4,  Mean@Transpose@data, 3, dat = data, _, Return[]];
		{zd, xd, yd} = Dimensions[dat];

		clipall = Ceiling[{0.5, zd - 0.5, 0.5, xd - .5, 0.5, yd - .5}];

		r1 = (vox[[2]]*xd)/(vox[[3]]*yd);
		r2 = (vox[[1]]*zd)/(vox[[3]]*yd);
		r3 = (vox[[1]]*zd)/(vox[[2]]*xd);

		size = Min[{r1, r2, r3}] 400;

		init = OptionValue[CropInit];
		init = If[ListQ[init] && Length[init]==6,
			{0, 0, xd+1, xd+1, 0, 0} + {1, 1, -1, -1, 1, 1} init[[{1,2,4,3,5,6}]],
			{1,zd,1,xd,1,yd}
		];

		cropwindow = DialogInput[{
			DefaultButton[],

			Manipulate[
				outp = Ceiling[{zmin, zmax, xd - xmax, xd - xmin, ymin, ymax}];

				Grid[{
					{Dynamic[Row[{"size: ", Ceiling[{zmax-zmin,xmax-xmin,ymax-ymin}]},"   "]]},
					{
						LocatorPane[Dynamic[{{ymin, xmax}, {ymax, xmin}}],
							Show[ArrayPlot[dat[[z]], ColorFunction -> "GrayTones", Frame -> False, AspectRatio -> r1, ImageSize -> size/r2],
								Graphics[{
									Red, Thick, Dynamic[Line[{{ymin, xmin}, {ymin, xmax}, {ymax, xmax}, {ymax, xmin}, {ymin, xmin}}]], 
									Green, Line[{{y - 0.5, -10}, {y - 0.5, xd + 10}}], 
									Blue, Line[{{-10, xd - x + 0.5}, {yd + 10, xd - x + 0.5}}], 
									Red, Dynamic[Circle[Mean[{{ymin, xmin}, {ymax, xmax}}], 2]]
								}], 
								PlotRange -> {{0, yd}, {0, xd}}
							], {{0.5, 0.5}, {yd - 0.5, xd - 0.5}}, Appearance -> Graphics[{Red, Disk[]}, ImageSize -> 10]
						]
					}, {
						LocatorPane[Dynamic[{{ymin, zmax}, {ymax, zmin}}],
							Show[ArrayPlot[Reverse[dat[[All, x]]], ColorFunction -> "GrayTones", Frame -> False, AspectRatio -> r2, ImageSize -> size/r2], 
								Graphics[{
									Blue, Thick, Dynamic[Line[{{ymin, zmin}, {ymin, zmax}, {ymax, zmax}, {ymax, zmin}, {ymin, zmin}}]], 
									Green, Line[{{y - 0.5, -10}, {y - 0.5, zd + 10}}], 
									Red, Line[{{-10, z - 0.5}, {yd + 10, z - 0.5}}], 
									Blue,  Dynamic[Circle[Mean[{{ymin, zmin}, {ymax, zmax}}], 2]]
								}], 
								PlotRange -> {{0, yd}, {0, zd}}
							], {{0.5, 0.5}, {yd - 0.5, zd - 0.5}}, Appearance -> Graphics[{Blue, Disk[]}, ImageSize -> 10]
						]
						,
						LocatorPane[Dynamic[{{xmin, zmax}, {xmax, zmin}}],
							Show[ArrayPlot[Reverse /@ Reverse[dat[[All, All, y]]], ColorFunction -> "GrayTones", Frame -> False, AspectRatio -> r3, ImageSize -> size/r3], 
								Graphics[{
									Green, Thick, Dynamic[ Line[{{xmin, zmin}, {xmin, zmax}, {xmax, zmax}, {xmax, zmin}, {xmin, zmin}}]],
									Blue, Line[{{x - 0.5, -10}, {x - 0.5, zd + 10}}], 
									Red, Line[{{-10, z - 0.5}, {xd + 10, z - 0.5}}], 
									Green, Dynamic[Circle[Mean[{{xmin, zmin}, {xmax, zmax}}], 2]]
								}], 
								PlotRange -> {{0, xd}, {0, zd}}
						], {{0.5, 0.5}, {xd - 0.5, zd - 0.5}}, Appearance -> Graphics[{Green, Disk[]}, ImageSize -> 10]]
					}
				}, Spacings -> 0]
					
					,
				{{z, Round[zd/2], "slice"}, 1, zd, 1},
				{{x, Round[xd/2], "row"}, 1, xd, 1},
				{{y, Round[yd/2], "column"}, 1, yd, 1},
				{{xmin, init[[3]] - 0.5}, 1, xmax - 1, ControlType -> None},
				{{xmax, init[[4]] - 0.5}, xmin + 1, xd, ControlType -> None},
				{{ymin, init[[5]] - 0.5}, 1, ymax - 1, ControlType -> None},
				{{ymax, init[[6]] - 0.5}, ymin + 1, yd, ControlType -> None},
				{{zmin, init[[1]] - 0.5}, 1, zmax - 1, ControlType -> None},
				{{zmax, init[[2]] - 0.5}, zmin + 1, zd, ControlType -> None},
				SynchronousUpdating->True

			]
		}, WindowTitle -> "Crop the data and press done", WindowFloating -> True, Modal -> True
	];

	dataout =If[!(OptionValue[CropOutput] === "Clip"),
	{a, b, c, d, e, f} = outp;
	ToPackedArray@N@If[dd == 3, 
		data[[a ;; b, c ;; d, e ;; f]], 
		data[[a ;; b, All, c ;; d, e ;; f]]]
	];

	output=Switch[OptionValue[CropOutput],
		"All",{dataout, outp},
		"Data",dataout,
		"Crop",outp];
	];

	Return[output]
]


(* ::Subsubsection::Closed:: *)
(*FindCrop*)


Options[FindCrop] = {CropPadding->5}

SyntaxInformation[FindCrop] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

FindCrop[dat_, OptionsPattern[]] := Block[{add, data, dim, d1, d2, unit, crp},
	add = {-1, 1} OptionValue[CropPadding];

	data = Unitize@Switch[ArrayDepth[dat],
		3, dat,
		4, dat[[All, 1]],
		_, Return[$Failed]
	];
	dim = Dimensions[data];

	d1 = Total[data];
	d2 = Total[data, {2}];
	unit = Unitize[{Total[d2, {2}], Total[d1, {2}], Total[d1]}];

	crp = MinMax[DeleteCases[# Range[Length[#]], 0]] + add & /@ unit;
	Flatten[MapThread[Clip[#1, {1, #2}] &, {crp, dim}]]
]


(* ::Subsubsection::Closed:: *)
(*AutoCropData*)


Options[AutoCropData] = {CropPadding->5}

SyntaxInformation[AutoCropData] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

AutoCropData[data_, opts:OptionsPattern[]] := Module[{crp},
	crp = FindCrop[data, opts];
	{ApplyCrop[data,crp], crp}
]


(* ::Subsubsection::Closed:: *)
(*ReverseCrop*)


SyntaxInformation[ReverseCrop] = {"ArgumentsPattern" -> {_, _, _, _.}};

ReverseCrop[data_, dim_, crop_] := ReverseCrop[data, dim, crop, {0, 0}]

ReverseCrop[data_, dimI_, crop_, {v1_, v2_}] := Module[{datac, pad, dim, p},

	dim = If[Length@Dimensions@data===4, dimI[[{1,3,4}]], dimI];

	pad = If[v1 === 0 && v2 === 0,
		(*use original crop*)
		Partition[Abs[{1, dim[[1]], 1, dim[[2]], 1, dim[[3]]} - crop], 2],
		(*use other voxel size*)
		Floor[(v1/v2) Partition[Abs[{1, dim[[1]], 1, dim[[2]], 1, dim[[3]]} - crop], 2]]
	];

	p = If[IntegerQ[Min@data], 0, 0.];
	datac = Switch[ArrayDepth[data],
		3, ArrayPad[data, pad, p],
		4, ArrayPad[data, Insert[pad, {0, 0}, 2], p],
		_, Return[$Failed, Module]
	];

	ToPackedArray@datac
]


(* ::Subsubsection::Closed:: *)
(*ApplyCrop*)


Options[ApplyCrop]={CropAlways->False}

SyntaxInformation[ApplyCrop] = {"ArgumentsPattern" -> {_, _, _.,OptionsPattern[]}};

ApplyCrop[data_, crop_?VectorQ, opts:OptionsPattern[]] := ApplyCrop[data, crop, {{1,1,1}, {1,1,1}}, opts]

ApplyCrop[data_, crop_?VectorQ , {v1_,v2_}, opts:OptionsPattern[]] := Module[{z1, z2, x1, x2, y1, y2, dim, out},

	out = ToPackedArray@data;
	dim = Dimensions[out];
	dim = If[Length[dim]==4, dim[[{1, 3, 4}]], dim];

	(*get crops coors*)
	{z1, z2, x1, x2, y1, y2} = Round[crop Flatten[{#, #} & /@ (v1/v2)]];

	If[OptionValue[CropAlways],
		{z1,z2} = Clip[{z1, z2}, {1, dim[[1]]}];
		{x1,x2} = Clip[{x1, x2}, {1, dim[[2]]}];
		{y1,y2} = Clip[{y1, y2}, {1, dim[[3]]}];
		,
		If[z1<1||z2>dim[[1]]||x1<1||x2>dim[[2]]||y1<1||y2>dim[[3]], Return[Message[ApplyCrop::dim]]]
	];

	ToPackedArray@Switch[ArrayDepth[out],
		4, out[[z1 ;; z2, All, x1 ;; x2, y1 ;; y2]],
		3, out[[z1 ;; z2, x1 ;; x2, y1 ;; y2]],
		2, out[[x1 ;; x2, y1 ;; y2]]
	]
]


(* ::Subsection:: *)
(*Split and merge*)


(* ::Subsubsection::Closed:: *)
(*CutData*)


SyntaxInformation[CutData] = {"ArgumentsPattern" -> {_,_.}}

CutData[data_] := CutData[data, FindMiddle[data, False]]

CutData[data_, print_?BooleanQ] := CutData[data, FindMiddle[data, print]]

CutData[data_, cut_?IntegerQ] := Switch[ArrayDepth[data],
	4,{data[[All, All, All, ;; cut]],data[[All, All, All, (cut + 1) ;;]],cut},
	3,{data[[All, All, ;; cut]], data[[All, All, (cut + 1) ;;]],cut}
]


FindMiddle[dati_, print_] := Module[{dat, fdat, len, datf,peaks,mid,peak,center,mask,ran,blur,i, max},
	(*flatten mean and normalize data*)
	dat=dati;
	fdat = Flatten[dat];
	dat = Clip[dat, {0, Quantile[Pick[fdat, Unitize[fdat], 1], .75]}];
	dat = N@Nest[Mean, dat, ArrayDepth[dat] - 1];
	max = 0.85 len;
	len = Length[dat];
	dat = max dat/Max[dat];
	mask = UnitStep[dat - .1 len];
	ran = Flatten[Position[mask, 1][[{1, -1}]]];

	peaks = {};
	blur = 20;
	i = 0;
	While[peaks === {} && i < 5,
		(*smooth the data a bit*)
		datf = max - GaussianFilter[mask dat, len/blur];
		(*find the peaks*)
		peaks = FindPeaks[datf];
		peaks = If[Length[peaks] >= 3, peaks[[2 ;; -2]], peaks];
		peaks = Select[peaks, (ran[[1]] < #[[1]] < ran[[2]]) &];
		blur += 10;
		i++;
	];

	If[peaks==={},
		Print["could not find the center."];
		$Failed,
		(*find the most middle peak*)
		mid = Round[Length[dat]/2];
		center = {mid, .75 max};
		peak = Nearest[peaks, center];

		If[print,Print[Show[
			ListLinePlot[{max-dat,datf}, PlotStyle->{Black,Orange}],
			ListPlot[{peaks,{center},peak},PlotStyle->(Directive[{PointSize[Large],#}]&/@{Blue,Gray,Green})]
			,ImageSize->75, Ticks -> None]
		]];
		(*output*)
		Round[First@First@peak]
	]
]


(* ::Subsubsection::Closed:: *)
(*StichData*)


SyntaxInformation[StichData] = {"ArgumentsPattern" -> {_,_}}

StichData[datal_, datar_] := RotateDimensionsLeft[Join[RotateDimensionsRight[datal], RotateDimensionsRight[datar]]];


(* ::Subsection::Closed:: *)
(*MakeIntFunction*)


Options[MakeIntFunction]={CenterVoxel->True, CenterRange->False}

SyntaxInformation[MakeIntFunction] = {"ArgumentsPattern" -> {_,_.,_.,OptionsPattern[]}};

MakeIntFunction[dat_, opts:OptionsPattern[]] := MakeIntFunction[dat, {1,1,1}, 1, opts]

MakeIntFunction[dat_, vox_?VectorQ, opts:OptionsPattern[]] := MakeIntFunction[dat, vox, 1, opts]

MakeIntFunction[dat_, int_?IntegerQ, opts:OptionsPattern[]] := MakeIntFunction[dat, {1,1,1}, int, opts]

MakeIntFunction[dat_, vox_?VectorQ, int_?IntegerQ, opts:OptionsPattern[]] := Block[{dim, def, range},
	dim = Dimensions[dat][[;;3]];
	range = Thread[{vox, vox dim}] - If[OptionValue[CenterVoxel]&&int>0, 0.5, 0] vox - If[OptionValue[CenterRange], 0.5, 0] dim;
	def = 0. dat[[1,1,1]];
	def =If[ListQ[def], Flatten@def, def];

	With[{
			ex = def,
			fdat = Flatten[dat, 3]
		},
		InterpolatingFunction[
			range,
			{5, If[ArrayDepth[dat]===3, 6, 2], 0, dim,
			{int, int, int} + 1, 0, 0, 0, 0, ex&, {}, {}, False},
			Range[range[[#, 1]], range[[#, 2]], vox[[#]]]& /@ {1, 2, 3},
			If[ArrayDepth[dat] === 3 && $VersionNumber >= 13.3, 
				{PackedArrayForm, Range[0,Length[fdat]], fdat}, 
				ToPackedArray@N@dat],
			{Automatic, Automatic, Automatic}
		]
	]
]


(* ::Subsection:: *)
(*Package Functions*)


(* ::Subsubsection::Closed:: *)
(*QMRIToolsPackages*)


SyntaxInformation[QMRIToolsPackages] = {"ArgumentsPattern" -> {}};

QMRIToolsPackages[] := DeleteDuplicates[(StringSplit[#, "`"] & /@ Contexts["QMRITools`*`"])[[All, 2]]]


(* ::Subsubsection::Closed:: *)
(*QMRIToolsFunctions*)


SyntaxInformation[QMRIToolsFunctions] = {"ArgumentsPattern" -> {_.,_.}};

QMRIToolsFunctions[] := QMRIToolsFunctions[QMRIToolsPackages[], 0];

QMRIToolsFunctions[toolb_?StringQ] := QMRIToolsFunctions[{toolb}, 0]

QMRIToolsFunctions[toolb : {_?StringQ ..}] := QMRIToolsFunctions[toolb, 0]

QMRIToolsFunctions[p_?IntegerQ] := QMRIToolsFunctions[QMRIToolsPackages[], p]

QMRIToolsFunctions[toolb_?StringQ, p_?IntegerQ] := QMRIToolsFunctions[{toolb}, p]

QMRIToolsFunctions[toolb : {_?StringQ ..} | {"All"}, p_?IntegerQ] := Block[{
		pack, ind, contexts, allNames, func, opts, names, out
	},

	{pack, ind} = If[toolb === {"All"}, {QMRIToolsPackages[], False}, {toolb, True}];
	contexts = Select["QMRITools`" <> # <> "`" & /@ pack, MemberQ[QMRITools`$Contexts, #] &];
	allNames = Names[# <> "*"] & /@ contexts;

	{func, opts} = Transpose[(
		names = #;
		opts = ToString /@ Sort[DeleteDuplicates[Flatten[Options[ToExpression[#]][[All, 1]] & /@ names]]];
		func = Complement[names, opts];
		{func, opts}
	) & /@ allNames];

	out = If[ind,
		Transpose[{pack, func, opts}],
		{{"QMRITools", Sort@DeleteDuplicates@Flatten@func, Sort@DeleteDuplicates@Flatten@opts}}
	];

	If[p > 0,
		PrintFuncList[#, p] & /@ out;,
		If[ind, out, First@out]
	]
]


PrintFuncList[{toolbox_, functions_, options_}, p_] := Block[{i, func, opt},
	func = If[Length[functions] > 0,
		i = Ceiling[Length[functions]/p];
		Transpose@Partition[functions, i, i, 1, ""],
		functions
	];

	opt = If[Length[options] > 0,
		i = Ceiling[Length[options]/p];
		Transpose@Partition[options, i, i, 1, ""],
		options
	];

	Print@Column[{
		"",
		Style[toolbox, Bold, 24],
		"",
		Style["Functions", Bold, 16],
		Grid[func, Alignment -> Left, ItemSize -> 18],
		"",
		Style["Options", Bold, 16],
		Grid[opt, Alignment -> Left, ItemSize -> 18],
		""
	}];
]


(* ::Subsubsection::Closed:: *)
(*QMRIToolsFuncPrint*)


SyntaxInformation[QMRIToolsFuncPrint] = {"ArgumentsPattern" -> {_.}};

QMRIToolsFuncPrint[]:=QMRIToolsFuncPrint[""]

QMRIToolsFuncPrint[toolb_String]:=If[toolb=="",PrintAll/@QMRIToolsFunctions[];,PrintAll[QMRIToolsFunctions[toolb]];]

PrintAll[{name_, functions_, options_}]:=(
	Print[Style[name, Bold, 24]];
	Print[Style["Functions", {Bold, 16}]];
	Print[Information[#]]& /@ functions;
	Print[Style["Options", {Bold, 16}]];
	Print[Information[#]]& /@ options;
);


(* ::Subsubsection::Closed:: *)
(*Compilable functions*)


SyntaxInformation[CompilableFunctions] = {"ArgumentsPattern" -> {}};

CompilableFunctions[]:=Block[{list1, list2, grids},
	(*based on https://mathematica.stackexchange.com/questions/1096/list-of-compilable-functions*)
	Internal`CompileValues[]; (*to trigger auto-load*)
	ClearAttributes[Internal`CompileValues, ReadProtected];

	list1 = Select[Compile`CompilerFunctions[], ! StringContainsQ[ToString[#], "`"] &];

	list2 = DownValues[Internal`CompileValues] /. HoldPattern[Verbatim[HoldPattern][Internal`CompileValues[sym_]] :> _] :> sym;
	list2 = Select[Complement[list2, list1], ! StringContainsQ[ToString[#], "`"] &];

	grids=Grid[Transpose@Partition[# // Sort, i=Ceiling[Length[#]/4], i, 1, ""], Alignment -> Left, ItemSize -> 15]&/@{list1, list2};

	Column[{
		Style["CompilerFunctions",Bold,16],
		"",
		grids[[1]],
		"",
		Style["CompileValues",Bold,16],
		"",
		grids[[2]]
	}]
]


(* ::Subsection:: *)
(*NoZeroFunctions*)


(* ::Subsubsection::Closed:: *)
(*DivideNoZero*)


SyntaxInformation[DivideNoZero] = {"ArgumentsPattern" -> {_,_,_.}};

DivideNoZero[numi_, deni_, "Comp"] := N@DivideNoZeroi[numi, deni]

DivideNoZero[numi_, deni_] := Re@N@DivideNoZeroi[numi, deni]

DivideNoZeroi = Compile[{{num, _Complex, 0}, {den, _Complex, 0}}, If[Abs[den] == 0., 0., num/den], 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> True];


(* ::Subsubsection::Closed:: *)
(*LogNoZero*)


SyntaxInformation[LogNoZero] = {"ArgumentsPattern" -> {_}};

LogNoZero[val_] := N[LogNoZeroi[val]]

LogNoZeroi = Compile[{{val, _Real, 0}},If[val == 0., 0., Log[val]], 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> False]


(* ::Subsubsection::Closed:: *)
(*ExpNoZero*)


SyntaxInformation[ExpNoZero] = {"ArgumentsPattern" -> {_}};

ExpNoZero[val_] := N[ExpNoZeroi[val]]

ExpNoZeroi = Compile[{{val, _Real, 0}},If[val == 0., 0., Exp[val]],
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> False]


(* ::Subsubsection::Closed:: *)
(*SignNoZero*)


SyntaxInformation[SignNoZero] = {"ArgumentsPattern" -> {_}};

SignNoZero[val_] := N[SignNoZeroi[val]]

SignNoZeroi = Compile[{{val, _Real, 0}}, Sign[Sign[val] + 0.0001],
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> False]


(* ::Subsubsection::Closed:: *)
(*MeanNoZero*)


SyntaxInformation[MeanNoZero] = {"ArgumentsPattern" -> {_}};

MeanNoZero[vec_] := MeanNoZeroi[If[ArrayDepth[vec] > 1, RotateDimensionsLeft[vec], vec]]

MeanNoZeroi = Compile[{{vec, _Real, 1}}, If[AllTrue[vec, # === 0. &], 0., Mean[Pick[vec, Unitize[vec], 1]]], 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*StandardDeviationNoZero*)


SyntaxInformation[StandardDeviationNoZero] = {"ArgumentsPattern" -> {_}};

StandardDeviationNoZero[vec_] := StandardDeviationNoZeroi[If[ArrayDepth[vec] > 1, RotateDimensionsLeft[vec], vec]]

StandardDeviationNoZeroi = Compile[{{vec, _Real, 1}}, If[AllTrue[vec, # === 0. &], 0., StandardDeviation[Pick[vec, Unitize[vec], 1]]], 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*MedianNoZero*)


SyntaxInformation[MedianNoZero] = {"ArgumentsPattern" -> {_}};

MedianNoZero[vec_] := MedianNoZeroi[If[ArrayDepth[vec] > 1, RotateDimensionsLeft[vec], vec]]

MedianNoZeroi = Compile[{{vec, _Real, 1}}, If[AllTrue[vec, # === 0. &], 0., Median[Pick[vec, Unitize[vec], 1]]], 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*RMSNoZero*)


SyntaxInformation[RMSNoZero] = {"ArgumentsPattern" -> {_}};

RMSNoZero[vec_] := RMSNoZeroi[If[ArrayDepth[vec] > 1, RotateDimensionsLeft[vec], vec]]

RMSNoZeroi = Compile[{{vec, _Real, 1}}, If[Total[vec] === 0., 0.,
	Sqrt[Mean[Pick[vec, Unitize[vec], 1]^2]]
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*MADNoZero*)


SyntaxInformation[MADNoZero] = {"ArgumentsPattern" -> {_}};

MADNoZero[vec_] := MADNoZeroi[If[ArrayDepth[vec] > 1, RotateDimensionsLeft[vec], vec]]

MADNoZeroi = Compile[{{vec, _Real, 1}}, Block[{vec2}, If[Total[vec] === 0.,	0.,
	vec2 = Pick[vec, Unitize[vec], 1];
	Median[Abs[vec2 - Median[vec2]]]
]],	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsection:: *)
(*Memory functions*)


(* ::Subsubsection::Closed:: *)
(*MBCount*)


MBCount[exp_] := Round[UnitConvert[N@Quantity[ByteCount[exp], "Byte"], "Megabytes"], .1];


(* ::Subsubsection::Closed:: *)
(*MemoryUsage*)


SyntaxInformation[MemoryUsage] = {"ArgumentsPattern" -> {_.}}

MemoryUsage[size_: 1] := (
	NotebookClose[memwindow];
	memwindow = CreateWindow[DialogNotebook[{
		CancelButton["Close", DialogReturn[]],
		Manipulate[
			If[listing === {},
				"Noting found",
				TableForm[Join[	#, 
						{
							Row[Dimensions[ToExpression["Global`" <> #[[1]]]], "x"], Head[ToExpression["Global`" <> #[[1]]]], 
							ClearButton["Global`" <> #[[1]]]
						}
					] & /@ listing,
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
)


SetAttributes[myByteCount, Listable ]

myByteCount[symbolName_String] := Replace[
	ToExpression[symbolName, InputForm, Hold], 
	Hold[x__] :> If[MemberQ[Attributes[x], Protected | ReadProtected],
		Sequence @@ {},
		(*output size in MB and name*)
		{
			StringDelete[symbolName, "Global`"],
			Round[ByteCount[Through[{OwnValues, DownValues, UpValues, SubValues, DefaultValues, FormatValues, NValues}[Unevaluated@x, Sort -> False]]]/1000000., .01]
		}
	]
];


(* ::Subsubsection::Closed:: *)
(*ClearTemporaryVariables*)


SyntaxInformation[ClearTemporaryVariables] = {"ArgumentsPattern" -> {_.}}

ClearTemporaryVariables[] := Block[{names, attr},
	names = Names["QMRITools`*`Private`*"];
	attr = Attributes /@ names;
	MapThread[If[#1 === {Temporary}, ClearAll[#2]] &, {attr, names}];
]


(* ::Subsection::Closed:: *)
(*SumOfSquares*)


Options[SumOfSquares] = {OutputWeights -> True}

SyntaxInformation[SumOfSquares] = {"ArgumentsPattern" -> {_,OptionsPattern[]}};

SumOfSquares[data_, OptionsPattern[]] := Block[{sos, weights, dataf},
	dataf = RotateDimensionsLeft[data];
	sos = SumOfSquaresi[dataf];

	If[OptionValue[OutputWeights],
		weights = DivideNoZero[dataf, sos];
		{sos, RotateDimensionsRight[weights]},
		sos
	]
]


SumOfSquaresi = Compile[{{sig, _Real, 1}}, Sqrt[Total[sig^2]], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> True];


(* ::Subsection::Closed:: *)
(*LLeastSquares*)


SyntaxInformation[LLeastSquares] = {"ArgumentsPattern" -> {_, _}};

LLeastSquares[ai_, y_]:=Block[{a},
	a = If[Length[y] == Length[ai], ai, Transpose[ai]];
	If[RealValuedNumberQ[Total[Flatten[y]]], 
		LLeastSquaresC[a, y], 
		If[RealValuedNumberQ[Total[Flatten[a]]],
			LLeastSquaresCC[a, y],
			LLeastSquaresCCC[a, y]
		]
	]
]	


LLeastSquaresC = Compile[{{A, _Real, 2}, {y, _Real, 1}}, 
	Inverse[Transpose[A] . A] . Transpose[A] . y,
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];

LLeastSquaresCC = Compile[{{A, _Real, 2}, {y, _Complex, 1}}, 
	Inverse[Transpose[A] . A] . Transpose[A] . y,
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];

LLeastSquaresCCC = Compile[{{A, _Complex, 2}, {y, _Complex, 1}}, 
	Inverse[ConjugateTranspose[A] . A] . ConjugateTranspose[A] . y,
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsection::Closed:: *)
(*NNLeastSquares*)


SyntaxInformation[NNLeastSquares] = {"ArgumentsPattern" -> {_, _}};

(*Main function*)
NNLeastSquares[a_, y_] := Block[{at, x, zeroed, w, zerow, pos, sp, xp, neg, xi, si, alpha, i, j ,l},
	(*Initialze values*)
	at = Transpose[a];(*already define Transpose A for speed up CalcW*)

	(*initialize values*)
	x = 0. a[[1]];
	l = Length[x];
	zeroed = x + 1.;

	(*initial vector w*)
	w = CalcW[a, at, y, x];

	(*first while loop: select highest positive solution in the zero set as long as the zero set is not empty*)
	j = 1; 
	While[j < l && Total[zeroed] > 0. && Max[zerow = zeroed w] > 0,
		j++;
		(*add the index of max(zerow) to the positive set (1-zeroed is positive zet)*)
		zeroed[[Last[Ordering[zerow]]]] = 0.;
		(*Calculate the LLS solution of the positive set*)
		pos = PosInd[zeroed];
		sp = LLSC[at[[pos]], y];

		(*recalculate the solutions of sp untill all values of sp are positive*)
		i = 1;
		While[i < l && Min[sp] < 0.,
			i++;
			(*calculated alpha, which is the minimal values of the rations xi/(xi-si) for all 
			negative values of s*)
			xp = x[[pos]];
			neg = UnitStep[-sp];
			xi = Pick[xp, neg, 1];
			si = Pick[sp, neg, 1];
			alpha = Min[xi/(xi - si)];

			(*removed the lowest value of sp from the posetive set*)
			zeroed[[Pick[pos, Unitize[Chop[(xp + alpha (sp - xp))]], 0]]] = 1.;
			(*recalculate the solution sp*)
			pos = PosInd[zeroed];
			sp = LLSC[at[[pos]], y];
		];

		(*set xp to sp and recalculate w*)
		x = 0. x; x[[pos]] = sp;
		w = CalcW[a, at, y, x];
	];
	x
]

PosInd = Compile[{{v, _Real, 1}}, Block[{z = Round@Total[1 - v]}, Ordering[v][[;; z]]], RuntimeOptions -> "Speed"];

CalcW = Compile[{{A, _Real, 2}, {At, _Real, 2}, {y, _Real, 1}, {x, _Real, 1}}, Chop[At . (y - A . x)], RuntimeOptions -> "Speed"];

LLSC = Compile[{{A, _Real, 2}, {y, _Real, 1}}, Inverse[A . Transpose[A]] . A . y, RuntimeOptions -> "Speed"];


(* ::Subsection::Closed:: *)
(*Filters*)


(* ::Subsubsection::Closed:: *)
(*LapFilter*)


LapFilter[data_, fil_:0.5] := Clip[Chop[ImageData[TotalVariationFilter[
	If[ArrayDepth[data]===3,Image3D[N@data, "Real"],Image[N@data, "Real"]],
	fil, Method -> "Laplacian", MaxIterations -> 30]]], MinMax[data]
]


(* ::Subsubsection::Closed:: *)
(*MedFilter*)


MedFilter[data_, fil_:1] := Clip[Chop[ImageData[MedianFilter[
	If[ArrayDepth[data]===3, Image3D[N@data, "Real"], Image[N@data, "Real"]],
	Round[fil]]]], MinMax[data]
]


(* ::Subsubsection::Closed:: *)
(*StdFilter*)


StdFilter[data_, ker_:2] := Abs[Sqrt[GaussianFilter[data^2, ker] - GaussianFilter[data, ker]^2]]


(* ::Subsection::Closed:: *)
(*GyromagneticRatio*)


GyromagneticRatio[nuc_]:=(nuc/.{"1H"->42.57747892,"2H"-> 6.536,"3He"-> -32.434,"7Li"->16.546,"13C"->10.7084,"14N"->3.077,"15N"-> -4.316,"17O"-> -5.772,
"19F"->40.052,"23Na"->11.262,"27Al"->11.103,"29Si"-> -8.465,"31P"->17.235,"57Fe"->1.382,"63Cu"->11.319,"67Zn"->2.669,"129Xe"-> 11.777})


(* ::Subsection::Closed:: *)
(*Squeeze*)


SyntaxInformation[Squeeze] = {"ArgumentsPattern" -> {_}}

Squeeze[data_] := ToPackedArray@First@Flatten[data, Flatten[Position[Dimensions[data], 1]]]


(* ::Subsection::Closed:: *)
(*DynamicPartition*)


SyntaxInformation[DynamicPartition] = {"ArgumentsPattern" -> {_,_,_.}}

(*partition data in lists of arbitrary length*)
DynamicPartition[l_, p : {__Integer}, x___] := dPcore[l, Accumulate@p, x] /; ! Negative@Min@p && Length@l >= Tr@p

(*Partition function*)
dPcore[l_, p : {q___, _}] := Inner[l[[# ;; #2]] &, {0, q} + 1, p, Head@l]
dPcore[l_, p_, All] := Append[dPcore[l, p], Drop[l, Last@p]]
dPcore[l_, p_, n__] := Join[dPcore[l, p], Partition[Drop[l, Last@p], n]]


(* ::Subsection:: *)
(*BSplineCurveFit*)


(* ::Subsubsection::Closed:: *)
(*BSplineCurveFit*)


Options[BSplineCurveFit] = {SplineDegree -> 2, SplineKnotsNumber -> 50, SplineRegularization -> 0};

SyntaxInformation[BSplineCurveFit] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

BSplineCurveFit[pts_, opts : OptionsPattern[]] := Block[{paras, knots, coeffMat, ctrlpts, cpn, sd, reg, len, Amat, ptsP},
	len = Length[pts];
	cpn = Min[{len - 2, OptionValue[SplineKnotsNumber]}];

	{coeffMat, Amat} = BSplineBasisFunctions[len, opts];
	ptsP = PadRight[pts, Length[Amat]];
	ctrlpts = LLeastSquares[Amat, ptsP];
	(Amat . ctrlpts)[[;; (-cpn - 1)]]
	]


(* ::Subsubsection::Closed:: *)
(*BSplineBasisFunctions*)


Options[BSplineBasisFunctions] = {SplineDegree -> 2, SplineKnotsNumber -> 50, SplineRegularization -> 0};

BSplineBasisFunctions[nPts_, opts : OptionsPattern[]] := BSplineBasisFunctions[nPts, opts] = Block[{
	cpn, sd, reg, len, paras, knots, coeffMat, coeffMatR, coeffMatDD, smooth
	},
	(*Get the options*)
	cpn = OptionValue[SplineKnotsNumber];
	sd = OptionValue[SplineDegree];
	reg = OptionValue[SplineRegularization];

	(*define the bpline points x = [0,1]*)
	paras = Range[0, 1, 1/(nPts - 1 + 2)] // N;
	(*define the knots for order sd and cpn degrees of freedome*)
	knots = Join[ConstantArray[0., sd], N@Range[0, 1, 1/(cpn - sd)], ConstantArray[1., sd]];

	(*generate the coefficient matrix*)
	coeffMat = Basis[sd, knots, paras];

	(*maker reg coefficient matirx*)
	coeffMatDD = ListConvolve[{1, -2, 1}, #] & /@ coeffMat;
	coeffMat = Transpose[coeffMat[[All, 2 ;; -2]]];
	smooth = reg (coeffMatDD . Transpose[coeffMatDD]);
	coeffMatR = Join[coeffMat, smooth];
	(*output*)
	{coeffMat, coeffMatR}
]


(* ::Subsubsection::Closed:: *)
(*Basis*)


(*generate b-spline basis functions with order p and knots, for x points*)
Basis[p_, knots_, x_] := Basis[p, knots, x] = Block[{kn, ui, ui1, uip, uip1, bi, bi1, d1, d2},
	(*function cashes the basis function already calculated*)
	If[p == 0,
		(*first order splines*)
		kn = knots;
		kn[[-1]] = kn[[-1]] + 1;
		(*get the 0th order basis function*)
		UnitComp[#[[1]], #[[2]], x] & /@ Partition[kn, 2, 1]
		,
		(*higher order splines, first partition the knots*)
		{ui, ui1, uip, uip1} = Transpose[Partition[knots, 2 + p, 1][[All, {1, 2, -2, -1}]]];
		(*get the basis functions of order p-1 and partition*)
		{bi, bi1} = Transpose[Partition[Basis[p - 1, knots, x], 2, 1]];
		(*get the basis functions of order p*)
		DivComp1[uip - ui, ui, x] bi + DivComp2[uip1 - ui1, uip1, x] bi1
	]
]


(*0th order b-spline basis function*)
UnitComp = Compile[{{min, _Real, 0}, {max, _Real, 0}, {x, _Real, 0}}, If[min <= x < max, 1, 0], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];

(*b-spline division i*)
DivComp1 = Compile[{{d, _Real, 0}, {u, _Real, 0}, {x, _Real, 1}}, If[d == 0., x, (x - u)/d], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];

(*b-spline division i+1*)
DivComp2 = Compile[{{d, _Real, 0}, {u, _Real, 0}, {x, _Real, 1}}, If[d == 0., x, (u - x)/d], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsection:: *)
(*RotationMatrix*)


(* ::Subsubsection::Closed:: *)
(*DecomposeAffineMatrix*)


SyntaxInformation[DecomposeAffineMatrix] = {"ArgumentsPattern" -> {_}}

DecomposeAffineMatrix[mat_] := Block[{T, L, R, S, Q},
	{T, L} = GetTranslation[mat];
	{R, S} = PolarDecomposition[L];
	{S, Q} = GetScaleSkew[S];
	N@{T, R, S, Q}
]


GetTranslation[mat_] := Block[{out = IdentityMatrix[4]},
	out[[All, 4]] = mat[[All, 4]];
	N@{out, Inverse[out] . mat}
]


PolarDecomposition[mat_] := Block[{R, S},
	{R, S} = {# . ConjugateTranspose[#3], #3 . #2 . ConjugateTranspose[#3]} & @@ SingularValueDecomposition[mat];
	If[Det[R] < 0, 
		R[[;; 3, ;; 3]] = -R[[;; 3, ;; 3]];
		S[[;; 3, ;; 3]] = -S[[;; 3, ;; 3]]];
	N@{R, S}
]


GetScaleSkew[s_] := Block[{sc},
	sc = (Norm /@ Transpose[s]) /. 0. -> 1.;
	N@{DiagonalMatrix[sc], Transpose[Transpose[s]/sc]}
]


(* ::Subsubsection::Closed:: *)
(*DecomposeScaleMatrix*)


SyntaxInformation[DecomposeScaleMatrix] = {"ArgumentsPattern" -> {_}}

DecomposeScaleMatrix[s_] := DeleteCases[N@MapThread[IdentityMatrix[4] + Transpose[{#2}] . {#2} (#1 - 1) &, Eigensystem[s]], N@IdentityMatrix[4]]


(* ::Subsubsection::Closed:: *)
(*QuaternionToRotationMatrix*)


SyntaxInformation[QuaternionToRotationMatrix] = {"ArgumentsPattern" -> {_}}

QuaternionToRotationMatrix[{a_?NumericQ, b_?NumericQ, c_?NumericQ, d_?NumericQ}] := N@{
	{a^2 + b^2 - c^2 - d^2, 2 b c - 2 a d, 2 b d + 2 a c},
	{2 b c + 2 a d, a^2 + c^2 - b^2 - d^2, 2 c d - 2 a b},
	{2 b d - 2 a c, 2 c d + 2 a b, a^2 + d^2 - c^2 - b^2}
};

SyntaxInformation[QuaternionToRotationMatrix] = {"ArgumentsPattern" -> {_}}

QuaternionVectorToRotationMatrix[{b_?NumericQ, c_?NumericQ, d_?NumericQ}] := QuaternionToRotationMatrix[{Sqrt[1 - b^2 - c^2 - d^2], b, c, d}];


(* ::Subsubsection::Closed:: *)
(*RotationMatrixToQuaternion*)


SyntaxInformation[RotationMatrixToQuaternion] = {"ArgumentsPattern" -> {_}}

RotationMatrixToQuaternion[{
	{r11_?NumericQ, r12_?NumericQ, r13_?NumericQ}, 
	{r21_?NumericQ, r22_?NumericQ, r23_?NumericQ}, 
	{r31_?NumericQ, r32_?NumericQ, r33_?NumericQ}}] := Block[{trace, a, b, c, d, x, y, z},

	trace = 1. + r11 + r22 + r33;
	If[trace > .5,
		a = .5 Sqrt[trace];
		b = .25 (r32 - r23)/a;
		c = .25 (r13 - r31)/a;
		d = .25 (r21 - r12)/a
		,
		x = 1. + r11 - (r22 + r33);
		y = 1. + r22 - (r11 + r33);
		z = 1. + r33 - (r11 + r22);

		Which[
			x > 1.,
			b = .5 Sqrt[x];
			a = .25 (r32 - r23)/b;
			c = .25 (r12 + r21)/b;
			d = .25 (r13 + r31)/b,

			y > 1.,
			c = .5 Sqrt[y];
			a = .25 (r13 - r31)/c;
			b = .25 (r12 + r21)/c;
			d = .25 (r23 + r32)/c,

			True,
			d = .5 Sqrt[z];
			a = .25 (r21 - r12)/d;
			b = .25 (r13 + r31)/d;
			c = .25 (r23 + r32)/d;
		];

		If[a < 0., {a, b, c, d} *= -1.]
	];
	N@{a, b, c, d}
];

SyntaxInformation[RotationMatrixToQuaternionVector] = {"ArgumentsPattern" -> {_}}

RotationMatrixToQuaternionVector[r : {{_?NumericQ, _?NumericQ, _?NumericQ}, {_?NumericQ, _?NumericQ, _?NumericQ}, {_?NumericQ, _?NumericQ, _?NumericQ}}] := Rest[RotationMatrixToQuaternion[r]];


(* ::Subsection:: *)
(*MakeFunctionGraph*)


(* ::Subsubsection::Closed:: *)
(*MakeFunctionGraph*)


Options[MakeFunctionGraph] = {
	LabelPlacement -> Tooltip,
	AllowSelfDependencies -> False
}

SyntaxInformation[MakeFunctionGraph] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

MakeFunctionGraph[func_, opts:OptionsPattern[]] := Block[{
		fName, lab, self, in, cont, flist, names, types, contexts, edges, 
		vertCol, vertFunc, vertLab
	},

	fName = StringSplit[ToString[#], "`"][[-1]] &;

	(*get options*)
	{lab, self} = OptionValue[{LabelPlacement, AllowSelfDependencies}];

	(*get list of all dependant functions*)
	in = {ToString@func};
	cont = True;
	While[cont, flist = DeleteDuplicates[Flatten[{in, GetDependencies /@ in}]];
		cont = If[in === flist, False, in = flist; True]];

	(*get properties of function list*)
	names = fName /@flist;

	types = Which[
		Head[ToExpression@#] === CompiledFunction, "Compiled",
		Head[ToExpression@#] === Function, "Function",
		True, "SetDelayed"
	]& /@flist;

	contexts = 	(cont = ToString[Context[#]];
		If[StringSplit[cont, "`"][[2]] === StringSplit[Context[func], "`"][[2]], "Internal", "External"] <> "_" <> If[StringContainsQ[cont, "Private"], "Private", "Global"]
	)& /@ flist;

	edges = DeleteDuplicates[Flatten[(f = #; DirectedEdge[fName[f], fName[#]] & /@ GetDependencies[f]) & /@ flist]];

	If[!self, edges = Select[edges, #[[1]] =!= #[[2]] &]];

	(*make graph properties*)
	vertCol = Thread[names -> (Directive[#, EdgeForm[None]] & /@ (contexts /. Thread[
		{"Internal_Global", "Internal_Private", "External_Global", "External_Private"} -> 
		{RGBColor[{52, 168, 83}/256], RGBColor[{66, 133, 244}/256],RGBColor[{251, 188, 5}/256], RGBColor[{235, 67, 53}/256]}]))
	];
	vertFunc = Thread[names -> types /. {"SetDelayed" -> "Circle", "Function" -> "Triangle", "Compiled" -> "Star"}];
	vertLab = Thread[names -> (Placed[#, lab] & /@ names)];

	Graph[edges, 
		VertexLabels -> vertLab, VertexShapeFunction -> vertFunc, VertexStyle -> vertCol,
		VertexLabelStyle -> Directive[Black, Bold, Automatic], EdgeStyle -> Directive[Black, Thick], 
		VertexSize -> Automatic, ImageSize -> {Automatic, 600}]
]


(* ::Subsubsection::Closed:: *)
(*GetDependencies*)


GetDependencies[sym_String] := If[
	Head[ToExpression[sym]] === Symbol,
	(*SymbolQ[ToExpression[sym]],*) 
	SelectFunctions[(Union@Level[(Hold @@ DownValues[sym])[[All, 2]], {-1}, Hold, Heads -> True])], Nothing]


SelectFunctions[func_] := Select[StringSplit[ToString /@ Pick[List @@ Defer /@ func, List@ReleaseHold[FunctionQC /@ func], 1],"[" | "]"][[All, 2]], Context[#] =!= "System`" &]

FunctionQC[f_Symbol] := If[(DownValues[f] =!= {}) && (OwnValues[f] === {}), 1, 0]
FunctionQC[f_] := If[Head[f] === Function || Head[f] === CompiledFunction, 1, 0]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
