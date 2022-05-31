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
Current assests are \"Elastix\", \"Transformix\" and \"DcmToNii\"." 

FileSelect::usage = 
"FileSelect[action] creates a systemdialog wicht returs file/foldername action can be \"FileOpen\", \"FileSave\" or \"Directory\".
FileSelect[action, {type}] same but allows the definition of filetypes for \"FileOpen\" and \"FileSave\" e.g. \"jpg\" or \"pdf\"."

TransData::usage = 
"TransData[data,dir] Rotates the dimesions of the data to left or rigthg. For example {z,x,y} to {x,y,z} dir is \"l\" or \"r\"."

RotateDimensionsLeft::usage = 
"RotateDimensionsLeft[data] rotates the dimensions of the data one to the left.
RotateDimensionsLeft[data, i] rotates the dimensions of the data i to the left."

RotateDimensionsRight::usage = 
"RotateDimensionsRight[data] rotates the dimensions of the data one to the right.
RotateDimensionsRight[data, i] rotates the dimensions of the data i to the right."

SaveImage::usage = 
"SaveImage[image] exports graph to image, ImageSize, FileType and ImageResolution can be given as options.
SaveImage[image, \"filename\"] exports graph to image with \"filname\", ImageSize, FileType and ImageResolution can be given as options."

FindMaxDimensions::usage = 
"FindMaxDimensions[{data1, data2, ..}] finds the maximal dimensions of all datasets. Each dataset is 3D."

PadToDimensions::usage = 
"PadToDimensions[data, dim] pads the data to dimensions dim." 

RescaleData::usage = 
"RescaleData[data,dim] rescales image/data to given dimensions.
RescaleData[data,{vox1, vox2}] rescales image/data from size vox1 to size vox2."

GridData::usage = 
"GridData[{data1,data2,...}, part] makes a grid of multiple datasets with part sets on each row."

GridData3D::usage = 
"GridData3D[{data1,data2,...}, part] same as grid data, but only works on 4D data where the data is gridded in axial, coronal and sagital."

Data2DToVector::usage = 
"Data2DToVector[data] converst the data to vector.
Data2DToVector[data,mask] converst the data within the mask to vector.

the data can be reconstructed using VectorToData.

output is the vecotrized data and a list contining the original data dimensions and a list with the data coordinates. {vec, {dim,pos}}."

Data3DToVector::usage = 
"Data3DToVector[data] converst the data to vector.
Data3DToVector[data,mask] converst the data within the mask to vector.

the data can be reconstructed using VectorToData.

output is the vecotrized data and a list contining the original data dimensions and a list with the data coordinates. {vec, {dim,pos}}."

DataToVector::usage = "
DataToVector[data] converst the non zero data to vector.
DataToVector[data,mask] converst the data within the mask to vector.

the data can be reconstructed using VectorToData.

output is the vecotrized data and a list contining the original data dimensions and a list with the data coordinates. {vec, {dim,pos}}."

VectorToData::usage = 
"VectorToData[vec, {dim,pos}] converts the vectroized data, using Data2DToVector or Data3DToVector, back to its original Dimensoins."

TensMat::usage=
"TensMat[tensor] transforms tensor form vector format {xx,yy,zz,xy,xz,yz} to matrix format {{xx,xy,xz},{xy,yy,yz},{xz,yz,zz}}."

TensVec::usage=
"TensVec[tensor] transforms tensor form matrix format {{xx,xy,xz},{xy,yy,yz},{xz,yz,zz}} to vector format {xx,yy,zz,xy,xz,yz}."


CropData::usage =
"CropData[data] creates a dialog window to crop the data (assumes voxsize (1,1,1)).
CropData[data,vox] creates a dialog window to crop the data."

FindCrop::usage = 
"FindCrop[data] finds the crop values of the data by removing all zeros surrounding the data."

AutoCropData::usage = 
"AutoCropData[data] crops the data by removing all background zeros.
AutoCropData[data,pad] crops the data by removing all background zeros with padding of pad."

ReverseCrop::usage = 
"ReverseCrop[data,dim,crop] reverses the crop on the cropped data with crop values crop to the original size dim.
ReverseCrop[data,dim,crop,{voxorig,voxnew}] reverses the crop on the cropped data with crop values crop to the original size dim."

ApplyCrop::usage =
"ApplyCrop[data,crop] aplies the corpped region obtained form CropData to the data.
ApplyCrop[data,crop,{voxorig,voxnew}] aplies the corpped region obtained form CropData to the data." 


CutData::usage = 
"CutData[data] splits the data in two equal sets left and right."

StichData::usage =
"StichData[datal,datar] joins left and right part of the data generated by CutData."


QMRIToolsPackages::usage =
"QMRIToolsPackages[] give list of all the QMRITools pacakges."

QMRIToolsFunctions::usage = 
"QMRIToolsFunctions[] give list of all the QMRITools packages, functions and options.
QMRIToolsFunctions[p] print a table with length p of all the QMRITools functions and options.
QMRIToolsFunctions[\"toobox\"] gives a list of all the functions and options in toolbox.
QMRIToolsFunctions[\"toobox\", p] gives a table op length p of all the functions and options in toolbox. If toolbox is \"All\" it will list all toolboxes."

QMRIToolsFuncPrint::usage = 
"QMRIToolsFuncPrint[] gives a list of all the QMRITools functions with their usage infomation."

CompilebleFunctions::usage = 
"CompilebleFunctions[] generates a formatted table of all compilable functions generated by Compile`CompilerFunctions."


DevideNoZero::usage = 
"DevideNoZero[a, b] devides a/b but when b=0 the result is 0. a can be a number or vector."

MeanNoZero::usage = 
"MeanNoZero[data] calculates the mean of the data ignoring the zeros."

MedianNoZero::usage = 
"MedianNoZero[data] calculates the Median of the data ignoring the zeros."

LogNoZero::usage = 
"LogNoZero[val] return the log of the val which can be anny dimonsion array. if val=0 the output is 0."

ExpNoZero::usage = 
"ExpNoZero[val] return the Exp of the val which can be anny dimonsion array. if val=0 the output is 0."

RMSNoZero::usage = 
"RMSNoZero[vec] return the RMS error of the vec which can be anny dimonsion array. if vec={0...} the output is 0. Zeros are ignored."

MADNoZero::usage = 
"MADNoZero[vec] return the MAD error of the vec which can be anny dimonsion array. if vec={0...} the output is 0. Zeros are ignored."


MemoryUsage::usage = 
"MemoryUsage[] gives a table of which definitions use up memory.
MemoryUsage[n] gives a table of which definitions use up memory, where n is the amout of definitions to show."

ClearTemporaryVariables::usage = 
"ClearTemporaryVariables[] Clear temporary variables."


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

LapFilter::usage = 
"LapFilter[data] Laplacian filter of data with kernel size 0.8.
LapFilter[data, ker] Laplacian filter of data with kernel ker."

StdFilter::usage =
"StdFilter[data] StandardDeviation filter of data using gaussian kernel 2. 
StdFilter[data, ker] StandardDeviation filter of data using kernel with size ker."


GyromagneticRatio::usage =
"GyromagneticRatio[] gives the gyromagnetic ratio for \"1H\" in MHz/T.
GyromagneticRatio[nucle] gives the gyromagnetir ratio for the nuclei, e.g. \"31P\" of \"1H\"."


Squeeze::usage =
"Squeeze[data] Revomes the singelton dimensions from data."

DynamicPartition::usage = 
"DynamicPartition[data, {part}] patitions the data into parts which is a list of integers. The remainders is los. 
DynamicPartition[data,part,last] patitions the data into parts which is a list of integers. The remainders is partitioned into equal parts defined by last.
If last is All, the remainders is just one partition."


BSplineCurveFit::usage = 
"BSplineCurveFit[points] fits a bspline to the points. Output is a list of same size as points."


DecomposeScaleMatrix::usage = 
"DecomposeScaleMatrix[mat] decomposes the affine matirx in T, R, S and Q."

DecomposeAffineMatrix::usage = 
"DecomposeAffineMatrix[S] decomposes the scale matrix in S1, S2 and S3."

QuaternionToRotationMatrix::usage =
"QuaternionToRotationMatrix[{a, b,c,d}] converts quarternion to rotation matrix R."

QuaternionVectorToRotationMatrix::usage =
"QuaternionVectorToRotationMatrix[{b,c,d}] converts quarternion to rotation matrix R."

RotationMatrixToQuaternion::usage = 
"RotationMatrixToQuaternion[R] converts rotation matrix to quarternions {a, b,c,d}."

RotationMatrixToQuaternionVector::usage =
"RotationMatrixToQuaternionVector[R] converts rotation matrix to quarternions {b,c,d}."


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


OutputWeights::usage = 
"OutputWeights is an option for SumOfSqares. If True it also output the SoS weights."


SplineKnotsNumber::usage =
"SplineKnotsNumber is an option for BSplineCurveFit and defines how many knots the bspline has."

SplineRegularization::usage =
"SplineRegularization is an option for BSplineCurveFit and defines the amount of regularization for the linear fit."


(* ::Subsection::Closed:: *)
(*Error Messages*)


RescaleData::dim = "Given dimensions `1` not the same depth as that of the given data `2`."

RescaleData::data = "Error: Inpunt must be 2D with {xdim,ydim} input or 3D dataset with {xdim,ydim} or {zdim, xdim, ydim} input."


Data2DToVector::dim = "Data should be 2D or 3D, data is `1`D."

Data2DToVector::mask = "Data and mask should have the same dimensions: data `1` and mask `2`."

Data3DToVector::dim = "Data should be 3D or 4D, data is `1`D."

Data3DToVector::mask = "Data and mask should have the same dimensions: data `1` and mask `2`."


ApplyCrop::dim = "Crop region lies outside data range."


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection:: *)
(*General Functions*)


GetAssetLocation[name_] := Block[{file},
	file=Last[SortBy[PacletFind["QMRITools"], #["Version"]]]["AssetLocation", name];
	If[FileExistsQ[file],file,$Failed]
]


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
  	Print["Canceled!"], 
  	If[action == "Directory",
  		StringDrop[input,-1],
  		input
  		]
  	]
]


(* ::Subsubsection::Closed:: *)
(*TransData*)


SyntaxInformation[TransData] = {"ArgumentsPattern" -> {_, _}};

TransData[dat_,dir_]:=Block[{data,ran,dep,fun},
	data = ToPackedArray[dat];
	ran = Range[dep=ArrayDepth[data]];
	fun = Switch[dir,"r",RotateLeft[ran],"l",RotateRight[ran]];
	Transpose[data,fun]
]



(* ::Subsubsection::Closed:: *)
(*RotateDimensionsLeft*)


SyntaxInformation[RotateDimensionsLeft] = {"ArgumentsPattern" -> {_, _.}};

RotateDimensionsLeft[dat_, i_ : 1] := RotateDimensions[dat, i, RotateRight]


(* ::Subsubsection::Closed:: *)
(*RotateDimensionsLeft*)


SyntaxInformation[RotateDimensionsRight] = {"ArgumentsPattern" -> {_, _.}};

RotateDimensionsRight[dat_, i_ : 1] := RotateDimensions[dat, i, RotateLeft]


(* ::Subsubsection::Closed:: *)
(*RotateDimensions*)


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
  	  Export[file, exp(*Rasterize[exp, ImageResolution -> 2*res],RasterSize -> imsize,*) , ImageSize->imsize, ImageResolution -> res,"ImageEncoding"->"LZW"],
  	  Export[file, exp(*Rasterize[exp, ImageResolution -> 2*res],RasterSize -> imsize,*) , ImageSize->imsize, ImageResolution -> res]
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

SyntaxInformation[PadToDimensions] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

PadToDimensions[data_, dim_, OptionsPattern[]] := Block[{diffDim, padval, pad,dir,zer},
  padval = OptionValue[PadValue];
  diffDim = dim - Dimensions[data];
  
  dir = OptionValue[PadDirection];
  dir = If[StringQ[dir], ConstantArray[dir, Length[dim]], dir];
  zer = ConstantArray[0, Length[dim]];
  
  pad = MapThread[
  Switch[#1, "Left", #2, "Right", #3, _, #4] &, {dir, 
   Transpose@{zer, diffDim}, Transpose@{diffDim, zer}, 
   Transpose@{Floor[diffDim/2], Ceiling[diffDim/2]}}];
  
  ToPackedArray[N@ArrayPad[data,pad,padval]]
	 
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
   3(*rescale a 3D dataset*),
   Switch[Length[sc],
    2(*rescale a stac of 2D images*),
    RescaleImgi[#, {sc, met}, int] & /@ data
    ,
    3(*rescale 3D data*),
    RescaleImgi[data, {sc, met}, int],
    _,
    Return[Message[RescaleData::dim, sc, Dimensions[data]]];
    ],
   4(*rescale a 4D dataset, treat data as multiple 3D sets*),
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

SyntaxInformation[GridData] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

GridData[dati_, part_, OptionsPattern[]] := Block[{dim, data, adepth, pad, val},
	adepth = ArrayDepth[dati[[1]]];
	
	(*pad the images with zeros*)
	pad = OptionValue[Padding];
	data = If[pad =!= None,
		{pad, val} = If[IntegerQ[pad], {pad, 0.}, pad];
		pad = PadLeft[{pad, pad}, adepth];
		ArrayPad[#, pad, val] & /@ dati
		,
		dati
	];
	
	(*make the first dimention such that it is devidable by part*)	
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
		{AX, COR, SAG}]
]


(* ::Subsubsection::Closed:: *)
(*Data2DToVector*)


SyntaxInformation[Data2DToVector] = {"ArgumentsPattern" -> {_, _.}};

Data2DToVector[datai_]:=Data2DToVector[datai,1]
Data2DToVector[datai_, maski_]:=Block[{data,depth,mask,pos,vecdata,dimd,dimm},

depth = ArrayDepth[datai];
If[!(depth==2||depth==3), Return@Message[Data2DToVector::dim,depth]];

dimd = If[depth==2,Dimensions[datai],Drop[Dimensions[datai],1]];
dimm = Dimensions[maski];
If[!(mask===1) && dimm!=dimd, Return@Message[Data2DToVector::mask,dimd,dimm]];

data = N@If[depth==3, maski #&/@datai, maski datai];

mask = Unitize@If[depth==3, Total[data], data];
pos = Position[mask,1];

vecdata = If[depth==3,
	DeleteCases[Flatten[TransData[data,"l"], 1], ConstantArray[0., Length[datai]]],
	DeleteCases[Flatten[data], 0.]
];

{vecdata, {dimd,pos}}
];


(* ::Subsubsection::Closed:: *)
(*Data3DToVector*)


SyntaxInformation[Data3DToVector] = {"ArgumentsPattern" -> {_, _.}};

Data3DToVector[datai_]:=Data3DToVector[datai,1]
Data3DToVector[datai_, maski_]:=Module[{data,depth,mask,pos,vecdata,dimd,dimm},

depth=ArrayDepth[datai];
If[!(depth==3||depth==4), Message[Data3DToVector::dim,depth]];

dimd = If[depth==3, Dimensions[datai], Drop[Dimensions[datai],1]];
dimm = Dimensions[maski];
If[!(mask===1) && dimm!=dimd, Return@Message[Data3DToVector::mask,dimd,dimm]];

data = N@If[depth==4,maski #&/@datai,maski datai];

mask = Unitize@If[depth==4,Total[data],data];
pos = Position[mask,1];

vecdata = If[depth==4,
	DeleteCases[Flatten[TransData[data,"l"], 2], ConstantArray[0., Length[datai]]],
	DeleteCases[Flatten[data], 0.]
];

{vecdata, {dimd,pos}}
];


(* ::Subsubsection::Closed:: *)
(*DataToVector*)


Clear[DataToVector]
DataToVector::dim = "`1` should be 2D, 3D or 4D, data is `2`D.";
DataToVector::mask = "Data and mask should have the same dimensions: data `1` and mask `2`";

SyntaxInformation[DataToVector] = {"ArgumentsPattern" -> {_, _.}};

DataToVector[datai_] := DataToVector[datai, 1]

DataToVector[datai_, maski_] := Module[{data, mask, depthd, depthm, depth, dimm, dimd},
	depthd = ArrayDepth[datai];
	If[! (depthd == 2 || depthd == 3 || depthd == 4), Return@Message[DataToVector::dim, "Data", depthd]];
	
	data = N[datai];	
	mask = If[maski === 1, Unitize[data], maski];
	
	depthm = ArrayDepth[mask];
	depth = depthd - depthm;
	
	(*data dimensions are not correct, mask must be 2D or 3D*)
	If[! (depthm == 2 || depthm == 3), Message[DataToVector::dim, "Mask", depthm]];
	
	dimm = Dimensions[mask];
	dimd = Dimensions[data];
	
	dimd = If[depth == 0, 
		(*mask and data are same dimensions*)
		dimd,
		(*data is one dimensions larger than mask either 2D and 3D or 3D and 4D*)
		If[depth == 1 && depthd == 3, Drop[dimd, 1], Drop[dimd, {2}]]
	];
	
	(*Dimensions must be equal*)
	If[ dimd =!= dimm, Return@Message[DataToVector::mask, dimd, dimm]];
	
	(*Flatten the data*)
	data = If[depthd == 4,
		Flatten[RotateDimensionsLeft[Transpose[data]], 2],
		If[depthd == 3 && depth == 1,
			Flatten[RotateDimensionsLeft[data], 1],
			Flatten[data]
		]
	];

	(*get the data and positions there mask is 1*)
	{Pick[data, Round[Flatten[mask]], 1] , {dimd, Position[mask, 1]}}
]


(* ::Subsubsection::Closed:: *)
(*VectorToData*)


SyntaxInformation[VectorToData] = {"ArgumentsPattern" -> {_, {_, _}}};

VectorToData[vec_, {dim_, pos_}] := If[VectorQ[vec],
	Normal[SparseArray[pos -> vec, dim]],
	If[Length[dim] == 2,
		Normal[SparseArray[pos -> #, dim]] & /@ Transpose[vec],
		Transpose[Normal[SparseArray[pos -> #, dim]] & /@ Transpose[vec]]
	]
]


(* ::Subsubsection::Closed:: *)
(*TensMat*)


SyntaxInformation[TensMat] = {"ArgumentsPattern" -> {_}};

TensMat[tens : {_?ArrayQ ..}] := RotateDimensionsLeft[TensMati[tens], 2];

TensMat[tens_?ListQ] := TensMati[tens]


TensMati[{xx_, yy_, zz_, xy_, xz_, yz_}] := {{xx, xy, xz}, {xy, yy, yz}, {xz, yz, zz}};


(* ::Subsubsection::Closed:: *)
(*TensVec*)


SyntaxInformation[TensVec] = {"ArgumentsPattern" -> {_}};

TensVec[tens : {_?ArrayQ ..}] := TensVeci[RotateDimensionsRight[tens, 2]]

TensVec[tens_?MatrixQ] := TensVeci[tens]


TensVeci[{{xx_, xy_, xz_}, {_, yy_, yz_}, {_, _, zz_}}] := {xx, yy, zz, xy, xz, yz};


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
     
    cropwindow = DialogInput[
      {
       DefaultButton[],

       Manipulate[
        outp = Ceiling[{zmin, zmax, xd - xmax, xd - xmin, ymin, ymax}];
        
        Grid[
         {
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
				], {{0.5, 0.5}, {yd - 0.5, zd - 0.5}}, Appearance -> Graphics[{Blue, Disk[]}, ImageSize -> 10]]
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
			}}, Spacings -> 0]
			
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
       
       }, WindowTitle -> "Crop the data and press done", 
      WindowFloating -> True, Modal -> True
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

FindCrop[data_, OptionsPattern[]] := Block[{unit, crp, p, dim, add},
	add= OptionValue[CropPadding];
	unit = Unitize[Total[Total[data, {#[[1]]}], {#[[2]]}] & /@ {{2, 2}, {1, 2}, {1, 1}}];
	dim = Dimensions[data];
	crp = (p = Position[#, 1]; Flatten@{First[p] - add, Last[p] + add}) & /@ unit;
	Flatten[MapThread[Clip[#1, {1, #2}] &, {crp, dim}]]
]


(* ::Subsubsection::Closed:: *)
(*AutoCropData*)


Options[AutoCropData] = {CropPadding->5}

SyntaxInformation[AutoCropData] = {"ArgumentsPattern" -> {_, OptionsPattern[]}}

AutoCropData[data_, opts:OptionsPattern[]] := Module[{datac,crp, add},
	add= OptionValue[CropPadding];
	datac = Switch[ArrayDepth[data],
		3, data,
		4, data[[All, 1]],
		_, Return[$Failed]
    ];
  
    crp = FindCrop[datac, opts];
    
    {ApplyCrop[data,crp],crp}
  ]


(* ::Subsubsection::Closed:: *)
(*ReverseCrop*)


SyntaxInformation[ReverseCrop] = {"ArgumentsPattern" -> {_, _, _, _.}};

ReverseCrop[data_, dim_, crop_] := ReverseCrop[data, dim, crop, {0, 0}]

ReverseCrop[data_, dim_, crop_, {v1_, v2_}] := Module[{datac, pad},
  
  pad = If[v1 === 0 && v2 === 0,
    (*use original crop*)
    Partition[Abs[{1, dim[[1]], 1, dim[[2]], 1, dim[[3]]} - crop], 2]
    ,
    (*use other voxel size*)
    Floor[(v1/v2) Partition[
       Abs[{1, dim[[1]], 1, dim[[2]], 1, dim[[3]]} - crop], 2]]
    ];

  datac = Switch[ArrayDepth[data],
    3, ArrayPad[data, pad],
    4, Transpose[ArrayPad[#, pad] & /@ Transpose[data]],
    _, Return[$Failed, Module]
    ];
    
    ToPackedArray@N@datac
]


(* ::Subsubsection::Closed:: *)
(*ApplyCrop*)


SyntaxInformation[ApplyCrop] = {"ArgumentsPattern" -> {_, _, _.}};

ApplyCrop[data_, crop_] := ApplyCrop[data, crop , {0,0}]

ApplyCrop[data_, crop_ , {v1_,v2_}] := Module[{z1, z2, x1, x2, y1, y2,dim},
	
	dim=Dimensions[data];
	dim=If[Length[dim]==4,dim[[{1,3,4}]],dim];
	
	{z1, z2, x1, x2, y1, y2} = If[v1===0&&v2===0,
		crop,
		Round[(crop - 1) Flatten[Transpose[ConstantArray[v1/v2, 2]]] + 1]
	];
	
	If[z1<1||z2>dim[[1]]||x1<1||x2>dim[[2]]||y1<1||y2>dim[[3]],Return[Message[ApplyCrop::dim]]];
		
  out = If[ArrayDepth[data] === 4,
   data[[z1 ;; z2, All, x1 ;; x2, y1 ;; y2]],
   If[ArrayDepth[data] === 3,
    data[[z1 ;; z2, x1 ;; x2, y1 ;; y2]],
    If[ArrayDepth[data] === 2,
     data[[x1 ;; x2, y1 ;; y2]]
     ]]];
     
     
     ToPackedArray@N@out
     
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
		3,{data[[All, All, ;; cut]], data[[All, All, (cut + 1) ;;]],cut}]

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

StichData[datal_, datar_] := TransData[Join[TransData[datal, "r"], TransData[datar, "r"]], "l"];


(* ::Subsection:: *)
(*Package Functions*)


(* ::Subsubsection::Closed:: *)
(*QMRIToolsPackages*)


SyntaxInformation[QMRIToolsPackages] = {"ArgumentsPattern" -> {}};

QMRIToolsPackages[] := DeleteDuplicates[(StringSplit[#, "`"] & /@ Contexts["QMRITools`*`"])[[All, 2]]]


(* ::Subsubsection::Closed:: *)
(*QMRIToolsFunctions*)


SyntaxInformation[QMRIToolsFunctions] = {"ArgumentsPattern" -> {_.,_.}};

QMRIToolsFunctions[]:=QMRIToolsFunctions[""];

QMRIToolsFunctions[toolb_String]:= Block[{packages,names,functions,options,allNames,output},
	packages = If[toolb === "", QMRIToolsPackages[], {toolb}];
	allNames = Sort[Flatten[Names["QMRITools`" <> # <> "`*"]]] & /@ packages;
		
	{functions, options} = Transpose[(names = #;
		options = ToString /@ Sort[DeleteDuplicates[Flatten[Options[ToExpression[#]][[All, 1]] & /@ names]]];
		functions = Complement[names, options];
		{functions, options}) & /@ allNames];
	
	output = Transpose[{packages, functions, options}];
	
	If[toolb === "", output, First@output]
]

QMRIToolsFunctions[p_Integer]:=Block[{toolbox,functions,options},
	{toolbox,functions,options}=Transpose[QMRIToolsFunctions[]];
	
	functions = Sort@DeleteDuplicates@Flatten[functions];
	options = Sort@DeleteDuplicates@Flatten[options];
	
	functions = If[Length[functions]<=p,{functions},Partition[functions, p, p, 1, ""]];
	options = If[Length[options]<=p,{options},Partition[options, p, p, 1, ""]];
	
	Print[Column[{"",Style["Functions", Bold, 16], "",functions // Transpose // TableForm,""}]];
	Print[Column[{"",Style["Options", Bold, 16], "",options // Transpose // TableForm,""}]];
]

QMRIToolsFunctions[toolb_String,p_Integer]:=Block[{toolbox,functions,options, output},
	If[toolb == "All",
		QMRIToolsFunctions[#, p] & /@ QMRIToolsPackages[];
		,
		{toolbox,functions,options}=QMRIToolsFunctions[toolb];
		
		functions = If[Length[functions]<=p,{functions},Partition[functions, p, p, 1, ""]];
		options = If[Length[options]<=p,{options},Partition[options, p, p, 1, ""]];
	
		output = Column[{"",
			Style[toolbox, Bold, 24], "",
			Style["Functions", Bold, 16], functions // Transpose// TableForm,"",
			Style["Options", Bold, 16], options// Transpose // TableForm,""}];
		
		Print[output];
	]
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


SyntaxInformation[CompilebleFunctions] = {"ArgumentsPattern" -> {}};

CompilebleFunctions[]:=(Partition[Compile`CompilerFunctions[] // Sort, 50, 50, 1, 1] // Transpose) /. 1 -> {} // TableForm


(* ::Subsection:: *)
(*NoZeroFunctions*)


(* ::Subsubsection::Closed:: *)
(*DevideNoZero*)


SyntaxInformation[DevideNoZero] = {"ArgumentsPattern" -> {_,_}};

DevideNoZero[numi_, deni_] := N@Chop@DevideNoZeroi[numi, deni]

DevideNoZeroi = Compile[{{num, _Complex, 0}, {den, _Complex, 0}}, If[Abs[den] == 0., 0., num/den], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> True];


(* ::Subsubsection::Closed:: *)
(*LogNoZero*)


SyntaxInformation[LogNoZero] = {"ArgumentsPattern" -> {_}};

LogNoZero[val_] := N[LogNoZeroi[val]]

LogNoZeroi = Compile[{{val, _Real, 0}},If[val == 0., 0., Log[val]], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> False]


(* ::Subsubsection::Closed:: *)
(*ExpNoZero*)


SyntaxInformation[ExpNoZero] = {"ArgumentsPattern" -> {_}};

ExpNoZero[val_] := N[ExpNoZeroi[val]]

ExpNoZeroi = Compile[{{val, _Real, 0}},If[val == 0., 0., Exp[val]],RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> False]


(* ::Subsubsection::Closed:: *)
(*MeanNoZero*)


SyntaxInformation[MeanNoZero] = {"ArgumentsPattern" -> {_}};

MeanNoZero[vec_] := MeanNoZeroi[If[ArrayDepth[vec] > 1, RotateDimensionsLeft[vec], vec]]

MeanNoZeroi = Compile[{{vec, _Real, 1}}, If[AllTrue[vec, # === 0. &], 0., Mean[Pick[vec, Unitize[vec], 1]]], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*MedianNoZero*)


SyntaxInformation[MedianNoZero] = {"ArgumentsPattern" -> {_}};

MedianNoZero[vec_] := MedianNoZeroi[If[ArrayDepth[vec] > 1, RotateDimensionsLeft[vec], vec]]

MedianNoZeroi = Compile[{{vec, _Real, 1}}, If[AllTrue[vec, # === 0. &], 0., Median[Pick[vec, Unitize[vec], 1]]], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];



(* ::Subsubsection::Closed:: *)
(*RMSNoZero*)


SyntaxInformation[RMSNoZero] = {"ArgumentsPattern" -> {_}};

RMSNoZero[vec_] := RMSNoZeroi[If[ArrayDepth[vec] > 1, RotateDimensionsLeft[vec], vec]]

RMSNoZeroi = Compile[{{vec, _Real, 1}}, If[AllTrue[vec, # === 0. &], 0., RootMeanSquare[Pick[vec, Unitize[vec], 1]]], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*MADNoZero*)


SyntaxInformation[MADNoZero] = {"ArgumentsPattern" -> {_}};

MADNoZero[vec_] := MADNoZeroi[If[ArrayDepth[vec] > 1, RotateDimensionsLeft[vec], vec]]

MADNoZeroi = Compile[{{vec, _Real, 1}}, If[AllTrue[vec, # === 0. &], 0., MedianDeviation[Pick[vec, Unitize[vec], 1]]], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


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
  dataf = TransData[data, "l"];
  sos = SumOfSquaresi[dataf];
  
  If[OptionValue[OutputWeights],
   weights = DevideNoZero[dataf, sos];
   {sos, TransData[weights,"r"]},
   sos
   ]
  ]

SumOfSquaresi = Compile[{{sig, _Real, 1}}, Sqrt[Total[sig^2]], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> True];


(* ::Subsection::Closed:: *)
(*LLeastSquares*)


SyntaxInformation[LLeastSquares] = {"ArgumentsPattern" -> {_, _}};

LLeastSquares[Ai_,y_]:=Block[{A},
	A = If[Length[y] == Length[Ai], Ai, Transpose[Ai]];
	If[RealQ[Total[Flatten[y]]], 
		LLeastSquaresC[A, y], 
		If[RealQ[Total[Flatten[A]]],
			LLeastSquaresCC[A, y],
			LLeastSquaresCCC[A, y]
		]
	]
]	


LLeastSquaresC = Compile[{{A, _Real, 2}, {y, _Real, 1}}, 
	Inverse[Transpose[A].A].Transpose[A].y,
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];

LLeastSquaresCC = Compile[{{A, _Real, 2}, {y, _Complex, 1}}, 
	Inverse[Transpose[A].A].Transpose[A].y,
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];

LLeastSquaresCCC = Compile[{{A, _Complex, 2}, {y, _Complex, 1}}, 
	Inverse[ConjugateTranspose[A].A].ConjugateTranspose[A].y,
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];

(* ::Subsection::Closed:: *)
(*NNLeastSquares*)


SyntaxInformation[NNLeastSquares] = {"ArgumentsPattern" -> {_, _}};

(*Main function*)
NNLeastSquares[A_, y_] := Block[{At, x, zeroed, w, zerow, pos, sp, xp, neg, xi, si, alpha, i, j ,l},
	(*Initialze values*)
	At = Transpose[A];(*already define Transpose A for speed up CalcW*)
	
	(*initialize values*)
	x = 0. A[[1]];
	l = Length[x];
	zeroed = x + 1.;
	
	(*initial vector w*)
	w = CalcW[A, At, y, x];
	
	(*first while loop: select highest positive solution in the zero set as long as the zero set is not empty*)
	j = 1; 
	While[j < l && Total[zeroed] > 0. && Max[zerow = zeroed w] > 0,
		j++;
		(*add the index of max(zerow) to the positive set (1-zeroed is positive zet)*)
		zeroed[[Last[Ordering[zerow]]]] = 0.;
		(*Calculate the LLS solution of the positive set*)
		pos = PosInd[zeroed];
		sp = LLSC[At[[pos]], y];
		
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
			sp = LLSC[At[[pos]], y];
		];
		
		(*set xp to sp and recalculate w*)
		x = 0. x; x[[pos]] = sp;
		w = CalcW[A, At, y, x];
	];
	x
]

PosInd = Compile[{{v, _Real, 1}}, Block[{z = Round@Total[1 - v]}, Ordering[v][[;; z]]], RuntimeOptions -> "Speed"];
CalcW = Compile[{{A, _Real, 2}, {At, _Real, 2}, {y, _Real, 1}, {x, _Real, 1}}, Chop[At.(y - A.x)], RuntimeOptions -> "Speed"];
LLSC = Compile[{{A, _Real, 2}, {y, _Real, 1}},Inverse[A.Transpose[A]].A.y, RuntimeOptions -> "Speed"];


(* ::Subsection::Closed:: *)
(*Filters*)


(* ::Subsubsection::Closed:: *)
(*LapFilter*)


LapFilter[data_, fil_:0.5] := Clip[Chop[ImageData[TotalVariationFilter[
	If[ArrayDepth[data]===3,Image3D[N@data, "Real"],Image[N@data, "Real"]],
	 fil, Method -> "Laplacian", MaxIterations -> 30]]], MinMax[data]]


(* ::Subsubsection::Closed:: *)
(*StdFilter*)


StdFilter[data_, ker_:2] := Abs[Sqrt[GaussianFilter[data^2, ker] - GaussianFilter[data, ker]^2]]


(* ::Subsection::Closed:: *)
(*GyromagneticRatio*)


GyromagneticRatio[nuc_]:=(nuc/.{"1H"->42.57747892,"2H"-> 6.536,"3He"-> -32.434,"7Li"->16.546,"13C"->10.7084,"14N"->3.077,"15N"-> -4.316,"17O"-> -5.772,
"19F"->40.052,"23Na"->11.262,"27Al"->11.103,"29Si"-> -8.465,"31P"->17.235,"57Fe"->1.382,"63Cu"->11.319,"67Zn"->2.669,"129Xe"-> 11.777})
gyro


(* ::Subsection::Closed:: *)
(*Squeeze*)


SyntaxInformation[Squeeze] = {"ArgumentsPattern" -> {_}}

Squeeze[data_] := Block[{single},
  single = ((1 - Unitize[Dimensions[data] - 1]) /. 0 -> All);
  While[single[[-1]] === All && Length[single] > 1, 
   single = Drop[single, -1]];
  ToPackedArray[data[[##]] & @@ single]
  ]


(* ::Subsection::Closed:: *)
(*DynamicPartition*)


SyntaxInformation[DynamicPartition] = {"ArgumentsPattern" -> {_,_,_.}}

(*partition data in lists of arbitrary length*)
DynamicPartition[L_, p : {__Integer}, x___] := dPcore[L, Accumulate@p, x] /; ! Negative@Min@p && Length@L >= Tr@p

(*Partition function*)
dPcore[L_, p : {q___, _}] := Inner[L[[# ;; #2]] &, {0, q} + 1, p, Head@L]
dPcore[L_, p_, All] := Append[dPcore[L, p], Drop[L, Last@p]]
dPcore[L_, p_, n__] := Join[dPcore[L, p], Partition[Drop[L, Last@p], n]]


(* ::Subsection::Closed:: *)
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
	(Amat.ctrlpts)[[;; (-cpn - 1)]]
	]


(* ::Subsubsection::Closed:: *)
(*BSplineBasisFunctions*)


Options[BSplineBasisFunctions] = {SplineDegree -> 2, SplineKnotsNumber -> 50, SplineRegularization -> 0};

BSplineBasisFunctions[Npts_, opts : OptionsPattern[]] := BSplineBasisFunctions[Npts, opts] = Block[{
	cpn, sd, reg, len, paras, knots, coeffMat, coeffMatR, coeffMatDD, smooth
	},
	(*Get the options*)
	cpn = OptionValue[SplineKnotsNumber];
	sd = OptionValue[SplineDegree];
	reg = OptionValue[SplineRegularization];
	
	(*define the bpline points x = [0,1]*)
	paras = Range[0, 1, 1/(Npts - 1 + 2)] // N;
	(*define the knots for order sd and cpn degrees of freedome*)
	knots = Join[ConstantArray[0., sd], N@Range[0, 1, 1/(cpn - sd)], ConstantArray[1., sd]];
	
	(*generate the coefficient matrix*)
	coeffMat = Basis[sd, knots, paras];
	
	(*maker reg coefficient matirx*)
	coeffMatDD = ListConvolve[{1, -2, 1}, #] & /@ coeffMat;
	coeffMat = Transpose[coeffMat[[All, 2 ;; -2]]];
	smooth = reg (coeffMatDD.Transpose[coeffMatDD]);
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


GetScaleSkew[S_] := Block[{sc},
	sc = (Norm /@ Transpose[S]) /. 0. -> 1.;
	N@{DiagonalMatrix[sc], Transpose[Transpose[S]/sc]}
]


(* ::Subsubsection::Closed:: *)
(*DecomposeScaleMatrix*)


SyntaxInformation[DecomposeScaleMatrix] = {"ArgumentsPattern" -> {_}}

DecomposeScaleMatrix[S_] := DeleteCases[N@MapThread[IdentityMatrix[4] + Transpose[{#2}] . {#2} (#1 - 1) &, Eigensystem[S]], N@IdentityMatrix[4]]


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

RotationMatrixToQuaternion[{{r11_?NumericQ, r12_?NumericQ, 
     r13_?NumericQ}, {r21_?NumericQ, r22_?NumericQ, 
     r23_?NumericQ}, {r31_?NumericQ, r32_?NumericQ, r33_?NumericQ}}] := 
  Block[
   {trace, a, b, c, d, x, y, z},
   
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


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
