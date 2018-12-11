(* ::Package:: *)

(* ::Title:: *)
(*DTITools ManipulationTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`ManipulationTools`", {"Developer`"}];

$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


RescaleData::usage = 
"RescaleData[data,dim] rescales image/data to given dimensions.
RescaleData[data,{vox1, vox2}] rescales image/data from size vox1 to size vox2."

Unwrap::usage = 
"Unwrap[data] unwraps the given dataset."

UnwrapSplit::usage = 
"UnwrapSplit[phase, data] unwarps the give phase dataset but splits the data into left and right using SplitData based in the data and performs the unwrapping seperately."

Correct::usage =
"Correct[data, phase, shiftpar] corrects the dataset data using the phasemap and the shiftpar and interpolation order 1.
Correct[data, phase, shiftpar, int] corrects the dataset data using the phasemap and the shiftpar and interpolation order int."

TensorCorrect::usage=
"TensorCorrect[tensor, phase, shift, vox] corrects the tensor based on B0 field map. Can perform both translation and rotation of tensor."

Deriv::usage = 
"Deriv[disp, vox] calculates the derivative of the displacement along the three main axes. disp is the displacement field, vox is the voxel size.
Deriv[disp, vox, mask] calculates the derivative of the displacement along the three main axes. Sharp edges between the background en disp are solved by the mask. mask is a mask delining the edge of the displacement field.";

JoinSets::usage =
"JoinSets[{dat1,dat2,...}, over] joins dat1, dat2, ... with over slices overlap.
JoinSets[{dat1,dat2,dat3...},{over1,over2,...}] joins dat1 and dat2 with over1 slices overlap, Joins dat2 and dat3 with over2 slices overlap and so on.
JoinSets[{dat1,dat2,...},{{over,drop1,drop2},...}] joins dat1, dat2 with over slices overlap and drops drop1 slices for dat1 and drop2 from drop 2."

SplitSets::usage = 
"SplitSets[data, Nsets, Nover] splits the data in Nsets with Nover slices overlap."

CorrectJoinSetMotion::usage =
"CorrectJoinSetMotion[[{dat1,dat2,...}, vox, over] motion correts multiple sets with overlap. Over is the number of slices overlap between stes. A Translation registration is performed."

TensMat::usage=
"TensMat[tensor] transforms tensor form vector format {xx,yy,zz,xy,xz,yz} to matrix format {{xx,xy,xz},{xy,yy,yz},{xz,yz,zz}}."

TensVec::usage=
"TensVec[tensor] transforms tensor form matrix format {{xx,xy,xz},{xy,yy,yz},{xz,yz,zz}} to vector format {xx,yy,zz,xy,xz,yz}."

GridDataPlot::usage = 
"GridDataPlot[{data1,data2,...}, part] makes a grid of multiple datasets with part sets on each row"

FindCrop::usage = 
"FindCrop[data] finds the crop values of the data by removing all zeros surrounding the data."

CropData::usage =
"CropData[data] creates a dialog window to crop the data (assumes voxsize (1,1,1)).
CropData[data,vox] creates a dialog window to crop the data."

ApplyCrop::usage =
"ApplyCrop[data,crop] aplies the corpped region obtained form CropData to the data.
ApplyCrop[data,crop,{voxorig,voxnew}] aplies the corpped region obtained form CropData to the data." 

ReverseCrop::usage = 
"ReverseCrop[data,dim,crop] reverses the crop on the cropped data with crop values crop to the original size dim.
ReverseCrop[data,dim,crop,{voxorig,voxnew}] reverses the crop on the cropped data with crop values crop to the original size dim."

AutoCropData::usage = 
"AutoCropData[data] crops the data by removing all background zeros.
AutoCropData[data,pad] crops the data by removing all background zeros with padding of pad."

TriggerGrid::usage =
"TriggerGrid[data, dyns, {{xmin, xmax}, {ymin, ymax}}]."

Data2DToVector::usage = 
"Data2DToVector[data] converst the data to vector.
Data2DToVector[data,mask] converst the data within the mask to vector.

the data can be reconstructed using VectorToData.

output is the vecotrized data and a list contining the original data dimensions and a list with the data coordinates. {vec, {dim,pos}}."

Data3DToVector::usage = 
"Data3DToVector[data] converst the data to vector..
Data3DToVector[data,mask] converst the data within the mask to vector.

the data can be reconstructed using VectorToData.

output is the vecotrized data and a list contining the original data dimensions and a list with the data coordinates. {vec, {dim,pos}}."

TransData::usage = 
"TransData[data,dir] Rotates the dimesions of the data to left or rigthg. For example {z,x,y} to {x,y,z} dir is \"l\" or \"r\"."

VectorToData::usage = 
"VectorToData[vec, {dim,pos}] converts the vectroized data, using Data2DToVector or Data3DToVector, back to its original Dimensoins"

DriftCorrect::usage = 
"DriftCorrect[data, bval] dirft corrects the data using the signals of the lowest bvalue that has 6 or more unique volumes.
For the function to work optimal it is best to have these volumes evenly spread througout thet data \
and for the first and last volume to have this low bvalue." 

CutData::usage = 
"CutData[data] splits the data in two equal sets left and right."

StichData::usage =
"StichData[datal,datar] joins left and right part of the data generated by CutData."

DataTranformation::usage = 
"DataTranformation[data,vox,w] transforms a 3D dataset accordint to the affine transformation vector w"

InvertDataset::usage = 
"InvertDataset[data] inverts the data along the x y and z axes. In other words it is rotated aroud the origin such that (x,y,z)=(-x,-y,-z) and (0,0,0)=(0,0,0)"

SortDiffusionData::usage = 
"SortDiffusionData[data, grad, bval] sorts the diffusion datasets grad and bval for magnitude of bvalue."

RemoveIsoImages::usage = 
"RemoveIsoImages[data, grad, bval] Romoves the ISO images from the philips scanner from the data. ISO images have g={0,0,0} and b>0."

ConcatenateDiffusionData::usage=
"ConcatenateDiffusionData[{{data1, .., dataN}, {grad1, .., gradN}, {bval, .., bvalN}, {vox, .., voxN}}] concatenates the diffusion data sets.
ConcatenateDiffusionData[{data1, .., dataN}, {grad1, .., gradN}, {bval, .., bvalN}, {vox, .., voxN}] concatenates the diffusion data sets."


(* ::Subsection::Closed:: *)
(*Options*)


Kernel::usage = 
"Kernel is an option for DeNoise. Kernel can be \"Gaussian\", \"Disk\" or \"Box\"."

ReverseData::usage =
"ReverseData is an option for JoinSets. Reverses each individual datset given as input for the JoinSets function. True by default."

ReverseSets::usage =
"ReverseSets is an option for JoinSets. Reverses the order of the datsets, False by default."

NormalizeSets::usage = 
"NormalizeSets is an option for JoinSets. True normalizes the individual stacs before joining."

MotionCorrectSets::usage = 
"MotionCorrectSets is an option for JoinSets. True motion corrects the individual stacs before joining using CorrectJoinSetMotion."

JoinSetSplit::usage = 
"JoinSetSplit is an option ofr CorrectJoinSetMotion. If True RegisterDataTransformSplit is used else RegisterDataTransform is used."

PaddOverlap::usage = 
"PaddOverlap is an option of CorrectJoinSetMotion and JoinSets. it allows for extra motion in the z direction."

MonitorUnwrap::usage = 
"MonitorUnwrap is an option for Unwrap and PhaseCalc. Monitor the unwrapping progress."

UnwrapDimension::usage = 
"UnwrapDimension is an option for Unwrap and PhaseCalc. Can be \"2D\" or \"3D\". 2D is for unwarpping 2D images or unwrapping the individual images from a 3D dataset \
(does not unwrap in the slice direction). 3D unwraps a 3D dataset in all dimensions."

RotationCorrect::usage =
"RotationCorrect is an option for TensorCorrect. Default is False. Is a tensor is deformed setting to True also the shear is accounted for by local rotation of the tensor"

DeNoiseIterations::usage = 
"DeNoiseIterations are the number of the denoising iterations."

MonitorDeNoise::usage = 
"MonitorDeNoise monitor the denoising progres."

NormalizeSignal::usage = 
"NormalizeSignal is an option for DriftCorrect."

CropOutput::usage = 
"CropOutput is an option for CropData, can be \"All\",\"Data\" or \"Crop\"."

CropInit::usage = 
"CropInit is an option for CropData. By default the crop is not initialized bu can be with {{xmin,xmax},{ymin,ymax},{zmin,zmax}}. "


(* ::Subsection::Closed:: *)
(*Error Messages*)


RescaleData::dim = "Given dimensions `1` not the same depth as that of the given data `2`."

RescaleData::data = "Error: Inpunt must be 2D with {xdim,ydim} input or 3D dataset with {xdim,ydim} or {zdim, xdim, ydim} input."

ApplyCrop::dim = "Crop region lies outside data range."

JoinSets::over = "Error: The overlap must be a number or a list which gives the overlap and how many slice must be droped. Not: `1`."

Unwrap::data2D = "Unwrapping Dimensions is 2D and this can only be preformed on 2D or 3D data. This data is `1`D."

Unwrap::data3D = "Unwrapping Dimensions is 3D and this can only be preformed on 3D data. This data is `1`D."

Unwrap::dim = "Unwrapping Dimensions can be \"2D\" or \"3D\", current value is `1`."

Data2DToVector::dim = "Data should be 2D or 3D, data is `1`D."

Data2DToVector::mask = "Data and mask should have the same dimensions: data `1` and mask `2`"

Data3DToVector::dim = "Data should be 3D or 4D, data is `1`D."

Data3DToVector::mask = "Data and mask should have the same dimensions: data `1` and mask `2`"

ConcatenateDiffusionData::dim= "data, grad and bval should be the same length:  data `1` / grad `2` / bval `2`"


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection::Closed:: *)
(*ConcatenateDiffusionData*)


SyntaxInformation[ConcatenateDiffusionData] = {"ArgumentsPattern" -> {_, _., _., _.}};

ConcatenateDiffusionData[data_?ListQ] :=If[Length[data] == 4,ConcatenateDiffusionData[data[[1]], data[[2]], data[[3]], data[[4]]]]

ConcatenateDiffusionData[data_, grad_, val_, vox_] := 
  Module[{dataout, gradout, valout, voxout},
   If[Length[data] == Length[grad] == Length[val],
    dataout = Transpose@Flatten[Transpose[NormalizeDiffData[#]] & /@ data, 1];
    gradout = Flatten[grad, 1];
    valout = Flatten[val];
    
    {dataout, gradout, valout} = RemoveIsoImages[dataout, gradout, valout];
    {dataout, gradout, valout} = SortDiffusionData[dataout, gradout, valout];
    ,
    Return[Message[ConcatenateDiffusionData::dim, Length[data],Length[grad], Length[val]]];
    ];
   
   voxout = If[ListQ[vox] && ! ListQ[vox[[1]]], vox, vox[[1]]];
   
   {dataout, gradout, valout, voxout}
   ];


(* ::Subsection::Closed:: *)
(*SortDiffusionData*)


SyntaxInformation[SortDiffusionData] = {"ArgumentsPattern" -> {_, _, _}};

SortDiffusionData[data_, grad_, val_] := Module[{pos, valu, sel},
  {valu, pos} = UniqueBvalPosition[val];
  sel = Flatten[pos];
  {data[[All, sel]], grad[[sel]], val[[sel]]}
  ]


(* ::Subsection::Closed:: *)
(*RemoveIsoImages*)


SyntaxInformation[RemoveIsoImages] = {"ArgumentsPattern" -> {_, _, _}};

RemoveIsoImages[data_, grad_, val_] := Module[{sel},
  sel = Complement[Range[Length[val]], 
    Complement[Flatten[Position[grad, {0., 0., 0.}]], 
     Flatten[Position[val, 0.]]]];
  {data[[All, sel]], grad[[sel]], val[[sel]]}
  ]


(* ::Subsection::Closed:: *)
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
    ]
]


(* ::Subsection::Closed:: *)
(*CropData*)


SyntaxInformation[AutoCropData] = {"ArgumentsPattern" -> {_,  _.}}

AutoCropData[data_, add_: 2] := Module[{datac,crp},
  datac = Switch[ArrayDepth[data],
    3, data,
    4, data[[All, 1]],
    _, Return[$Failed, Module]
    ];
  
  crp=Flatten@{
    FindCropVals[datac, add],
    FindCropVals[TransData[datac, "l"], add],
    FindCropVals[TransData[datac, "r"], add]
    };
    
    {ApplyCrop[data,crp],crp}
  ]

FindCropVals[data_, add_] := Module[{pos(*, partpos, diff, postr*)},
  (*pos = Flatten@Position[Total[Flatten[N@#]] & /@ data, 0.];*)
  pos = Unitize[Total[Flatten[N@#]] & /@ data];
  If[pos === {},
   {1, Length[data]},
   	Clip[Flatten[{FirstPosition[pos, 1], Abs[FirstPosition[Reverse@pos, 1] - Length[pos] - 1]}] + {-add, add}, {1, Length[pos]}]
	(*
	partpos = Partition[pos, 2, 1];
	diff = Subtract @@@ partpos;
	postr = DeleteCases[diff + 1, 0] - 1;
	postr = Flatten@Position[diff, #] & /@ postr;
	postr = Flatten[partpos[[#]] & /@ postr];
	postr[[{1, -1}]] + {-add, add}
	*)
   ]
  ]


(* ::Subsection::Closed:: *)
(*FindCrop*)

SyntaxInformation[FindCrop] = {"ArgumentsPattern" -> {_}};

FindCrop[data_] := Block[{unit, crp, p, dim},
  unit = Unitize[
    Total[Total[data, {#[[1]]}], {#[[2]]}] & /@ {{2, 2}, {1, 2}, {1, 
       1}}];
  dim = Dimensions[data];
  crp = (
      p = Position[#, 1]; Flatten@{First[p] - 2, Last[p] + 2}
      ) & /@ unit;
  Flatten[MapThread[Clip[#1, {1, #2}] &, {crp, dim}]]
]


(* ::Subsection::Closed:: *)
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
 If[dd == 3, 
data[[a ;; b, c ;; d, e ;; f]], 
data[[a ;; b, All, c ;; d, e ;; f]]
]
];
output=Switch[OptionValue[CropOutput],
"All",{dataout, outp},
"Data",dataout,
"Clip",outp];
];
Return[output]
]


(* ::Subsection::Closed:: *)
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
		
  If[ArrayDepth[data] === 4,
   data[[z1 ;; z2, All, x1 ;; x2, y1 ;; y2]],
   If[ArrayDepth[data] === 3,
    data[[z1 ;; z2, x1 ;; x2, y1 ;; y2]],
    If[ArrayDepth[data] === 2,
     data[[x1 ;; x2, y1 ;; y2]]
     ]]]]



(* ::Subsection::Closed:: *)
(*Rescale Data*)


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
   
   Chop[Clip[dataOut,MinMax[data]]]
  ]

RescaleImgi[dat_, {sc_, met_}, n_] := Block[{type, im, dim},
  (*data type*)
  type = If[ArrayQ[dat, _, IntegerQ], "Bit16", "Real32"];
  dim = If[met == "v", Round[sc Dimensions[dat]], sc];
  (*convert to 2D or 3D image*)
  im = Switch[ArrayDepth[dat], 2, Image[dat, type], 3, Image3D[dat, type]];
  ImageData[ImageResize[im, Reverse[dim], Resampling ->{"Spline", n},Padding->0], type]
  ]


(* ::Subsubsection::Closed:: *)
(*GridData*)


SyntaxInformation[GridDataPlot] = {"ArgumentsPattern" -> {_, _}};

GridDataPlot[data_, part_] := Block[{dim, temp, adepth},
	adepth = ArrayDepth[data[[1]]];
	dim = Dimensions[data];
	dim[[1]] = dim[[1]] + (part - (Mod[Length[data], part] /. 0 -> part));
	temp = Transpose[Partition[PadRight[data, dim], part]];
	temp = MapThread[Join, #, adepth - 2] & /@ temp;
	temp = MapThread[Join, temp, adepth - 1]
  ]


(* ::Subsubsection::Closed:: *)
(*TriggerGrid*)


SyntaxInformation[TriggerGrid] = {"ArgumentsPattern" -> {_, _, _}};

TriggerGrid[data_, dyns_, {{r11_, r12_}, {r21_, r22_}}] := 
 Module[{tmp1, tmp1b, tmp2},
  tmp1 = Partition[Transpose[data[[All, All, r11 ;; r12, r21 ;; r22]]], 4];
  tmp1b = Map[N[#/Mean[Flatten[#]]] &, tmp1, {3}];
  tmp2 = Flatten[Flatten[Transpose[tmp1b, {1, 3, 2, 4, 5}], {2, 3}],1];
  GridDataPlot[tmp2, dyns]
  ]


(* ::Subsection:: *)
(*Phase unwrap*)


(* ::Subsubsection:: *)
(*Unwrap*)


Options[Unwrap]={MonitorUnwrap->True,UnwrapDimension->"2D"};

SyntaxInformation[Unwrap] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

(* Phase unwrapping algorithem based on M.A. Herraez et al 2002 and Abdul-Rahman 2007.
unwraping algorithem. Unwraps one image, needs functions: Diff, SecDiff, EdgeReliability and PairCoor*)

Unwrap[dat_,OptionsPattern[]]:=
Module[{data,step,undim,out,mon},
	
	undim=OptionValue[UnwrapDimension];
	mon=OptionValue[MonitorUnwrap];
	
	Switch[undim,
		
		"2D",
		data = N[dat];
		Which[
			MatrixQ[data],
			If[mon,PrintTemporary["Unwrapping one image using 2D algorithm."]];
			out = Unwrapi[data];
			,
			ArrayQ[data,3],
			If[mon,PrintTemporary["Unwrapping ",Length[data]," images using 2D algorithm"]];
			Monitor[
				out = MapIndexed[( step = First[#2]; Unwrapi[#1] )&, data, {ArrayDepth[data]-2} ];
				out = UnwrapZi[out];
			 ,If[mon,ProgressIndicator[step, {0, Length[data]}],""]
			]
		];,
		
		"3D",
		data = N[ArrayPad[dat,1]];
		If[ArrayQ[data,3],
			If[mon,PrintTemporary["Unwrapping 3D data using 3D algorithm"]];
			out = ArrayPad[Unwrapi[data, mon],-1];
			,
			Message[Unwrap::data3D,ArrayDepth[data]]
			];,
		_,
		Message[Unwrap::dim,undim]
		];
		
	(*center around 0*)
	out
	]


(* ::Subsubsection:: *)
(*UnwrapSplit*)


Options[UnwrapSplit] = Options[UnwrapSplit]

SyntaxInformation[UnwrapSplit] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

UnwrapSplit[phase_, mag_,opts:OptionsPattern[]] := Module[{cutVal, phaseSplit, B0split},
  cutVal = CutData[mag][[3]];
  phaseSplit = CutData[phase, cutVal][[1 ;; 2]];
  B0split = Unwrap[#, opts] & /@ phaseSplit;
  StichData @@ B0split
  ]


(* ::Subsubsection:: *)
(*UnwrapZi*)


UnwrapZi[data_]:=
Module[{
	Roundi=Module[{num2,tresh=.75},
	If[Negative[#],
		num2=#-Ceiling[#];
		If[-1<num2<-tresh,Floor[#],Ceiling[#]],
		num2=#-Floor[#];
		If[1>num2>tresh,Ceiling[#],Floor[#]]
		]
		]&,
		mask,slice,diff,meandiff,steps,off,unwrap,dat},
	
	mask=Unitize[data];
	slice=Round[0.5Length[data]];
	diff=(#[[1]]-#[[2]]&/@ Partition[data/(2Pi),2,1]);
	
	meandiff=Median[#]&/@ Map[DeleteCases[Flatten[N[#]],0.]&,diff];
	steps=FoldList[Plus,0,Map[Roundi[#]&,meandiff]];
	off=Round[Median[DeleteCases[Flatten[N[data[[slice]]/(2Pi)]],0.]]];
	
	unwrap=steps-(steps[[slice]]+off);
	dat=(2Pi unwrap+data)mask//N;
	(dat - mask Round[MeanNoZero[Flatten[dat]],2Pi])
	]


(* ::Subsubsection:: *)
(*Unwrapi*)


Unwrapi[dat_] := Unwrapi[dat, False]

Unwrapi[dat_, mon_] := Block[{data,datai,mask, crp, dimi, sorted,groups,groupsize,groupnr,task},
	(*monitor*)
	task="Preclustering data.";
	If[mon,PrintTemporary[Dynamic[task]]];
	
	(*rescale the dat to integers allows for faster matrix replacements*)
	datai = Round[10000. * dat / (2. Pi)];
	dimi = Dimensions[datai];
	
	(*remove zeros*)
	data = If[ArrayDepth[datai] == 3, crp = FindCrop[datai]; ApplyCrop[datai, crp], datai];
	
	(*make mask to pervent unwrapping in background*)
	mask = ArrayPad[Closing[ArrayPad[Mask[Ceiling[Abs@data], 1], 5], 1], -5];
	
	(*Get the edges sotrted for reliability and precluster groups*)
	sorted = GetEdgeList[data, mask];
	{groups, groupsize, groupnr} = MakeGroups[data, mask];

	(*make 2D data 3D and define shifts in add*)
	If[ArrayDepth[data] == 2,
		groups = {groups}; data = {data};
		sorted = {#[[1]] + 1, 1, #[[2]], #[[3]]} & /@ sorted;
		];
	
	(*Unwrap the data*)
	task="Unwrapping edges.";
	data = UnWrapC[sorted, data, groups, groupsize, groupnr];
	
	(*make output in rad*)	
	If[ArrayDepth[dat] == 2, 
		(*output the 2D in rad*)
		2 Pi data[[1]]/10000.,
		(*align to zero and ouput 3D in rad*)
		data = 2 Pi (data - mask Round[MeanNoZero[Flatten[data]],10000])/10000.;
		ReverseCrop[data, dimi, crp]
		]
]


UnWrapC = Compile[{{sorted, _Integer, 2}, {datai, _Integer, 3}, {groupsi, _Integer, 3}, {groupsizei, _Integer, 1}, {groupnri, _Integer, 0}},
	Block[{data, const, dir, dim, groups, group1, group2, groupsize, groupnr, z1, z2, x1, x2, y1, y2, wrap, wrapT, pos, g1, g2, out, adds,add},
   	
    groups = groupsi;
    groupsize = groupsizei;
    data = datai;
    groupnr = groupnri;

    dim = Dimensions[data];
    adds = {{1, 0, 0}, {0, 0, 1}, {0, 1, 0}};
    z1=0;x1=0;y1=0;z2=0;x2=0;y2=0;
    group1=0;group2=0;g1=0;g2=0;

    out = Map[(

        (*get the voxel corrdinates and the neighbour, contrain to dimensions*)
        add = adds[[#[[1]]]];
        z1=#[[2]]; z2=z1+add[[1]]; If[z2>dim[[1]],z2=dim[[1]]];
        x1=#[[3]]; x2=x1+add[[2]]; If[x2>dim[[2]],x2=dim[[2]]];
        y1=#[[4]]; y2=y1+add[[3]]; If[y2>dim[[3]],y2=dim[[3]]];
                
        (*Get the group numbers*)
        group1 = groups[[z1, x1, y1]];
        group2 = groups[[z2, x2, y2]];
        
        (*unwrapping logic*)
        (*0. check if both are not in background if one is background skip*)
        If[! (group1 == 0 || group2 == 0),
         (*1. both already in same group but not zero (most cases 60%) do nothing*)
         If[(! (group1 > 1 && group1 == group2)),
          (*If not in the same group determine the wrap of the edge and determine how to unwrap*)
          wrap = Round[Round[data[[z1, x1, y1]] - data[[z2, x2, y2]],5000], 10000];
          wrapT = (wrap != 0);
          
          (*get group sizes*)
          g1 = groupsize[[group1]];
          g2 = groupsize[[group2]];
          
          Which[
           (*2. one of two pixels already in group othter not in group (second most cases 30%)*)
           (*2A. group 1 existst add group 2 to group 1*)
           group2 == 1 && group1 > 1, 
           (*add the non group voxel to the existing group*)
           groups[[z2, x2, y2]] = group1;
           groupsize[[group1]] += 1;
           If[wrapT, data[[z2, x2, y2]] += wrap];
           ,
           (*2B. group 2 existst add group 1 to group 2*)
           group1 == 1 && group2 > 1, 
           (*add the non group voxel to the existing group*)
           groups[[z1, x1, y1]] = group2;
           groupsize[[group2]] += 1;
           If[wrapT, data[[z1, x1, y1]] -= wrap];
           ,
           (*3. both belong to no group (third most cases 6%)*)
           group1 == group2 == 1, 
           (*unwrap right or botom pixel and assign both to group*)
           groupnr++;
           groups[[z1, x1, y1]] = groups[[z2, x2, y2]] = groupnr;
           groupsize[[groupnr]] = 2;
           If[wrapT, data[[z2, x2, y2]] += wrap];
           ,
           (*4. both already in a group (least cases 4%)*)
           (*4A. group 1 is greather than or eaqual to group 2,
           add group 2 to group 1*)
           g1 >= g2, 
           (*unwrap group 2 with respect to group 1,
           only do if wrap\[NotEqual]0*)
           pos = 1 - Unitize[groups - group2];
           If[wrapT, data += wrap pos];
           groups += ((group1 - group2) pos);
           groupsize[[group1]] += g2;
           groupsize[[group2]] = 0;
           ,
           (*4B. group 2 is greather than or eaqual to group 1,
           add group 1 to group 2*)
           g2 > g1, 
           (*unwrap group 1 with respect to group 2,
           only do if wrap\[NotEqual]0*)
           pos = 1 - Unitize[groups - group1];
           If[wrapT, data -= wrap pos];
           groups += ((group2 - group1) pos);
           groupsize[[group2]] += g1;
           groupsize[[group1]] = 0;
           ]]];
        1
        ) &, sorted];
    data],
    {{dim, _Integer, 1},{add, _Integer, 1}, {group1, _Integer, 0}, {group2, _Integer, 0}, {g1, _Integer, 0}, {g2, _Integer, 0}, 
    	{z1, _Integer, 0}, {x1, _Integer, 0}, {y1, _Integer, 0}, {z2, _Integer, 0}, {x2, _Integer, 0}, {y2, _Integer, 0}},
    RuntimeOptions -> "Speed", Parallelization -> True];

(* ::Subsubsection::Closed:: *)
(*GetKernels*)


GetKernels[dep_] := Block[{ker, i, j, k, keri, kers},
   Switch[dep,
    2,
    ker = ConstantArray[0, {3, 3}];
    ker[[2, 2]] = -1;
    kers = ({i, j} = #; keri = ker; keri[[i, j]] = 1; keri) & /@ {
       {2, 1}, {2, 3}, {1, 2}, {3, 2}, {1, 1}, {3, 3}, {1, 3}, {3, 1}
       }(*H,V,D1,D2*),
    3,
    ker = ConstantArray[0, {3, 3, 3}];
    ker[[2, 2, 2]] = -1;
    kers = ({i, j, k} = #; keri = ker; keri[[i, j, k]] = 1; 
        keri) & /@ {
       {2, 1, 2}, {2, 3, 2}, {1, 2, 2}, {3, 2, 2}, {2, 2, 1}, {2, 2, 3},(*H,V,N*)
       {1, 1, 2}, {3, 3, 2}, {1, 3, 2}, {3, 1, 2},(*plane 1*)
       {1, 2, 1}, {3, 2, 3}, {1, 2, 3}, {3, 2, 1},(*plane 2*)
       {2, 1, 1}, {2, 3, 3}, {2, 1, 3}, {2, 3, 1},(*plane 3*)
       {1, 1, 1}, {3, 3, 3}, {1, 1, 3}, {3, 3, 1}, {1, 3, 1}, {3, 1, 3}, {3, 1, 1}, {1, 3, 3}(*diagonals*)
       }];
   Total /@ Partition[kers, 2]
   ];


(* ::Subsubsection::Closed:: *)
(*GetEdgeList*)


GetEdgeList[data_] := GetEdgeList[data, 1]

GetEdgeList[data_, maski_] := Block[{dep, diff, mask, edge, coor, fedge, ord, pos},
	dep = ArrayDepth[data];
	(*maske a mask if needed*)
	mask = If[maski === 1, Closing[Mask[Ceiling[Abs@data], 1], 1], maski];
	(*calculate the second order diff*)
	diff = ListConvolve[#, data, ConstantArray[2, dep], 0] & /@ GetKernels[dep];
	diff = If[dep == 2, DiffC2[diff], DiffC3[diff]];
	
	(*get the edge reliability*)
	edge = Switch[dep,
		2(*2D data*), 
		N@{(RotateLeft[#] & /@ diff) + diff, RotateLeft[diff] + diff},
		3(*3D data*), 
		N@{RotateLeft[diff] + diff, ((RotateLeft[#] & /@ #) & /@ diff) + diff, (RotateLeft[#] & /@ diff) + diff}
	];
	edge = mask # & /@ edge;
	 
	(*sort the edges for reliability*)
	coor = MapIndexed[#2 &, edge, {dep + 1}];
	fedge = Flatten[edge, dep];
	ord = Ordering[fedge];
	pos = Position[Unitize[fedge[[ord]]], 1, 1, 1][[1, 1]];
	Flatten[coor, dep][[ord]][[pos ;;]]
  ]

DiffC2 = Compile[{{diff, _Real, 3}}, Total[(diff - Round[diff, 10000])^2]];
DiffC3 = Compile[{{diff, _Real, 4}}, Total[(diff - Round[diff, 10000])^2]];


(* ::Subsubsection:: *)
(*MakeGroups*)


MakeGroups[data_] := MakeGroups[data, 1]
MakeGroups[data_, maski_]:=Block[{dep,dim,fun,min,max,part,dat,masks,mclus,clus,groupsize,groups,groupnr,mask},
	(*get data properties*)
	dep=ArrayDepth[data];
	dim=Dimensions[data];
	fun=If[dep==2,Image,Image3D];
	
	(*maske a mask if needed*)
	mask=If[maski===1,Closing[Mask[Ceiling[Abs@data],1],1],maski];
	
	(*find mask ranges*)
	{min,max}=MinMax[data];
	part=Partition[Range[min,max,(max-min)/6]//N,2,1];
	
	(*remove background form masks, and create masks*)
	dat=data/. 0->(-2min);
	masks=Mask[dat,#]&/@part;
	
	(*make groups from masks*)
	mclus=0;
	groups=Total[(
		(*make groups from masks and only keeps groups*)
		clus=MorphologicalComponents[DeleteSmallComponents[fun[#]]];
		clus=Clip[clus,{0,1}](mclus+clus);
		(*find max group number and export*)
		mclus=Max[clus];clus
	)&/@masks];
	groups=groups+mask;
	
	(*create outputs, the size vector and group nrs*)
	groupsize=ConstantArray[0,Count[Flatten[groups],1]];
	(groupsize[[#]]=Count[groups,#,2])&/@Range[1,Max[groups]];
	groupnr=Max[groups];
	{groups,groupsize,groupnr}
]


(* ::Subsection::Closed:: *)
(*Correct Phase*)


CorrectPhase[phase_,left_,right_]:=
Module[{dim,ones,add,out},
	dim=Dimensions[phase];
	ones=ConstantArray[1,{dim[[1]],dim[[2]],dim[[3]]/2}];
	add=Join[left*ones,right*ones,3];
	out=N[phase+add];
	out/.{N[left]->0.,N[right]->0.}
	]


(* ::Subsection:: *)
(*Correct*)


(* ::Subsubsection::Closed:: *)
(*Correct*)


SyntaxInformation[Correct] = {"ArgumentsPattern" -> {_, _, _, _.}};

Correct[data_?MatrixQ,phase_?MatrixQ,shift_]:=
Correcti[data,phase,shift,1]

Correct[data_?MatrixQ,phase_?MatrixQ,shift_,int_]:=
Correcti[data,phase,shift,int]

Correct[data:{_?MatrixQ..},phase:{_?MatrixQ..},shift_]:=
MapThread[Correcti[#1,#2,shift,1]&,{data,phase}]

Correct[data:{_?MatrixQ..},phase:{_?MatrixQ..},shift_,int_]:=
MapThread[Correcti[#1,#2,shift,int]&,{data,phase}]

Correct[data:{{_?MatrixQ..}..},phase:{_?MatrixQ..},shift_]:=
Transpose[Map[MapThread[Correcti[#1,#2,shift,1]&,{#,phase}]&,Transpose[data,{2,1}]],{2,1}]

Correct[data:{{_?MatrixQ..}..},phase:{_?MatrixQ..},shift_,int_]:=
Transpose[Map[MapThread[Correcti[#1,#2,shift,int]&,{#,phase}]&,Transpose[data,{2,1}]],{2,1}]


(* ::Subsubsection::Closed:: *)
(*Correcti*)


Correcti[dat_,ph_,shift_,int_]:=
Module[{pos,acpos,shiftpx,data,phase,output},
	If[shift[[2]]=="COL",
		data=Transpose[dat];phase=Transpose[ph];,
		data=dat;phase=ph;
		];
	shiftpx=phase*shift[[1]];
	output=Round[MapThread[(
		pos=Range[Length[#1]];
		acpos=pos-#1;
		ListInterpolation[#2,InterpolationOrder->int][acpos]
		)&,{shiftpx,data}]];
	If[shift[[2]]=="COL",Return[Transpose[output]],Return[output]]
	]


(* ::Subsection:: *)
(*TensorCorrect*)


(* ::Subsubsection::Closed:: *)
(*TensorCorrect*)


Options[TensorCorrect]={RotationCorrect->False};

SyntaxInformation[TensorCorrect] = {"ArgumentsPattern" -> {_, _, _, _, _., OptionsPattern[]}};

(* zonder masker, dus met sprongen in de afgeleide by grens tussen deformatie veld en achtergrond *)
TensorCorrect[tens_,phase_,shift_,vox_,OptionsPattern[]]:=
TensorCorrect[tens,phase,0,shift,vox];

(* met masker, dus zonder sprongen in de afgeleide by grens tussen deformatie veld en achtergrond *)
TensorCorrect[tens_,phase_,mask_,shift_,vox_,OptionsPattern[]]:=
	Module[{dim,pxshift,der,F,tensM,tensC,tensCV,tensT},
	
	dim=Dimensions[phase];
	(*deformation expessed in pixels*)
	pxshift=phase*shift[[1]];
	
	If[OptionValue[RotationCorrect]==True,
		PrintTemporary["Cacluating Derivative"];
		(*local derivative of the displacement in the slice direction*)
		der=If[!ArrayQ[mask],
		Deriv[pxshift,vox],
		Deriv[pxshift,vox,mask]
		];
		F=Fmat[der,shift[[2]]];
		
		PrintTemporary["Rotation Correction"];
		(*rotation correction of matrix*)
		(*tensor to matrixform*)
		tensM=TensMat[tens];
		(*rotation correct tensor matrix*)
		tensC=MapThread[DRot[#1,#2]&,{tensM,F},3];
		(*corrected tensor back to vector form*)
		tensCV=TensVec[tensC];
		,
		tensCV=tens;
		];
	
	PrintTemporary["Translation Correction"];
	(*Translation correction of the rotation corrected Tensor*)
	tensT=Map[(
		MapThread[
			TransCorrect[#1,#2,shift[[2]],1]
			&,{#,pxshift}
			]
		)&,tensCV]
	];


(* ::Subsubsection::Closed:: *)
(*TransCorrect*)


(* Translation correct one slice*)
TransCorrect[dat_,sh_,dir_,int_]:=
Module[{data,shift,pos,acpos,out},
	(*Transpose the data zo the deformation is always in the "ROW" direction*)
	If[dir=="COL",
		data=Transpose[dat];shift=Transpose[sh];,
		data=dat;shift=sh;
		];
	(*{dims,dimx,dimy}=Dimensions[data];*)
	(*deformation Correction*)
	out=MapThread[(
		pos=Range[Length[#1]];
		acpos=pos-#1;
		ListInterpolation[#2,InterpolationOrder->int][Clip[acpos,{1,Length[acpos]}]]
		(*[{Clip[acpos[[1]],{1,dimx}],Clip[acpos[[2]],{1,dimy}]}]*)
		)&,{shift,data}];
	
	(*If deformation was in the "COL" direction rotate back*)
	If[dir=="COL",Return[Transpose[Chop[out]]],Return[Chop[out]]]
	];


(* ::Subsubsection::Closed:: *)
(*FMat*)


Fmat[der_,shift_]:=
Module[{Dx,Dy,Dz,dim,zero,ones,F},
	{Dx,Dy,Dz}=der;
	dim=Dimensions[Dx];
	zero=ConstantArray[0,dim];
	ones=ConstantArray[1,dim];
	If[shift=="COL",
		F=Transpose[{{ones,zero,zero},{Dx,Dy+1,Dz},{zero,zero,ones}},{4,5,1,2,3}],
		If[shift=="ROW",
			F=Transpose[{{Dx+1,Dy,Dz},{zero,ones,zero},{zero,zero,ones}},{4,5,1,2,3}]
			],
		Print["error, unknown direction"]
		]
	];


(* ::Subsubsection::Closed:: *)
(*Drot*)


DRot[D_,F_]:=Module[{val,e1,e2,e3,n1,n2,n3,NN},
{val,{e1,e2,e3}}=Eigensystem[D];
n1=Normalize[F.e1];
n2=Normalize[F.e2-(n1.(F.e2))*n1]//N;
n3=Normalize[Cross[n1,n2]]//N;
NN=Transpose[{n1,n2,n3}];
Chop[NN.(IdentityMatrix[3]val).Transpose[NN]]
];


(* ::Subsection:: *)
(*Deriv*)


(* ::Subsubsection::Closed:: *)
(*Deriv*)


SyntaxInformation[Deriv] = {"ArgumentsPattern" -> {_, _, _}};

Deriv[disp_,vox_]:=
Module[{dim,Dx,Dy,Dz},
	dim=Dimensions[disp];
	Dx=Transpose[DerivFunc[Transpose[disp,{1,3,2}],dim[[2]],vox[[2]]],{1,3,2}];
	Dy=DerivFunc[disp,dim[[3]],vox[[3]]];
	Dz=Transpose[DerivFunc[Transpose[disp,{3,2,1}],dim[[1]],vox[[1]]],{3,2,1}];
	{Dx,Dy,Dz}
	];

Deriv[disp_,vox_,mask_]:=
Module[{dim,Dx,Dy,Dz},
	dim=Dimensions[disp];
	Dx=Transpose[DerivFunc[Transpose[disp,{1,3,2}],Transpose[mask,{1,3,2}],dim[[2]],vox[[2]]],{1,3,2}];
	Dy=DerivFunc[disp,mask,dim[[3]],vox[[3]]];
	Dz=Transpose[DerivFunc[Transpose[disp,{3,2,1}],Transpose[mask,{3,2,1}],dim[[1]],vox[[1]]],{3,2,1}];
	{Dx,Dy,Dz}
	];


(* ::Subsubsection::Closed:: *)
(*DerivFunc*)


DerivFunc[disp_,length_,step_]:=
Module[{coor,f},
	coor=Range[length]*step;
	Map[(
		f=Interpolation[Transpose[{coor,#}],InterpolationOrder->1];
		Table[f'[st],{st,coor}]
		)&,disp,{2}]
	];

DerivFunc[disp_,mask_,length_,step_]:=
Module[{coor,f,fr,df,dfr},
	coor=Range[length]*step;
	MapThread[(
		f=Interpolation[Transpose[{coor,#1}],InterpolationOrder->1];
		fr=Interpolation[Transpose[{coor,Reverse[#1]}],InterpolationOrder->1];
		df=RotateRight[#2,1]*#2*Table[f'[st],{st,coor}];
		dfr=RotateLeft[#2,1]*#2*-Reverse[Table[fr'[st],{st,coor}]];
		MapThread[If[#1==0,#2,#1]&,{df,dfr}]
		)&,{disp,mask},2]
	];


(* ::Subsection::Closed:: *)
(*SplitSets*)


Options[SplitSets] = {ReverseSets -> False, ReverseData -> True, PaddOverlap->0};

SyntaxInformation[SplitSets] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

SplitSets[data_, sets_, overlap_, OptionsPattern[]] := Module[{lengthSet, sels, start, end, dat, over, pad},
  
  dat = If[OptionValue[ReverseData], Reverse[data], data];
  
  pad = OptionValue[PaddOverlap];  
  dat = ArrayPad[dat,{{pad,pad},{0,0},{0,0}}];
  
  over=overlap+2pad;
  
  lengthSet = (Length[dat] + (sets - 1)*over)/sets;
  sels = Table[
    start = (i lengthSet + 1) - i over;
    end = start + lengthSet - 1;
    Range[start, end]
    , {i, 0, sets - 1}];
  
  dat = (dat[[#]] & /@ sels);
  
  dat = If[OptionValue[ReverseData], Reverse[dat, 2], dat];
  dat = If[OptionValue[ReverseSets], Reverse[dat], dat];
  
  dat
  ]


(* ::Subsection:: *)
(*Join sets*)


(* ::Subsubsection::Closed:: *)
(*JoinSets*)


Options[JoinSets]={ReverseSets->True,ReverseData->True, NormalizeSets -> True, MotionCorrectSets -> False, PaddOverlap -> 2};

SyntaxInformation[JoinSets] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

JoinSets[data_?ArrayQ,over_,opts:OptionsPattern[]]:=JoinSets[data,over,{1,1,1},opts]

JoinSets[data_?ArrayQ,over_,vox_,OptionsPattern[]]:=Block[
	{dat, overlap, motion, pad, normalize, depth, meth, target},
	
	(*get the options*)
	motion = OptionValue[MotionCorrectSets];
	pad = OptionValue[PaddOverlap];
	normalize=OptionValue[NormalizeSets];
	depth=ArrayDepth[data];
	overlap = If[ListQ[over],First@over,over];
	
	(*normalize the data*)
	dat=If[normalize,
		PrintTemporary["normalizing data"];
		Switch[ArrayDepth[data],
			5,100 NormalizeDiffData/@data,
			4,100 NormalizeData/@data
		],
		data
	];

	(*reverse the order of the sets if needed*)
	dat=If[OptionValue[ReverseSets],Reverse[dat],dat];
	
	If[motion,
		Switch[depth,
			5,
			motion=False;
			(*define the moving data*)
			Print["motion correct is only for 3D volues"]
			,
			4,
			PrintTemporary["motion correcting data"];
			dat = CorrectJoinSetMotion[dat, vox, over, PaddOverlap->pad];
			overlap = overlap + 2*pad;
		]
	];
	
	(*reverse the order of the slices if needed*)
	dat=N@If[OptionValue[ReverseData],Reverse[dat,2],dat];
	
	PrintTemporary["Joining data"];	
	dat = Switch[depth,
		5,Transpose[(JoinSetsi[dat[[All, All, #]],overlap]) & /@ Range[Length[dat[[1, 1]]]]],
		4,JoinSetsi[dat,overlap],
		_,$Failed
	];
	
	(*give output*)	
	dat = If[motion, ArrayPad[dat, Prepend[ConstantArray[{0, 0}, ArrayDepth[dat] - 1], {-pad, -pad}]],dat];
	
	Return[If[OptionValue[ReverseData],Reverse[dat],dat]]
]


JoinSetsi[data_?ArrayQ,overlap_?IntegerQ]:=
Module[{sets,set1,set2,step,set1over,set2over,joined},
	
	sets=Length[data];
	step=1/(overlap+1);
	
	(*perform the join*)
	For[i=1,i<sets,i++,
		If[i==1,
			set1=Drop[data[[i]],{-overlap,-1}];
			set1over=Take[data[[i]],{-overlap,-1}];
			,
			set1=Drop[joined,{-overlap,-1}];
			set1over=Take[joined,{-overlap,-1}];
			];
		set2=Drop[data[[i+1]],{1,overlap}];
		set2over=Take[data[[i+1]],{1,overlap}];

		joined = Joini[{set1, set2}, {set1over, set2over}, overlap];
		];
	
	joined	

	];


JoinSetsi[data_?ArrayQ,overlap_?ListQ,OptionsPattern[]]:=
Module[{sets,set1,set2,i,step,set1over,set2over,joined,overSet,data1,data2,drop1,drop2,overl},
	
	sets=Length[data];
	
	(*perform the join*)
	For[i=1,i<sets,i++,
		overSet=overlap[[i]];
		If[i==1,
			data1=data[[i]];,
			data1=joined;
			];
		If[Length[overSet]!=3&&!IntegerQ[overSet],
			Return[Message[JoinSets::over,overSet]];
			,
			If[IntegerQ[overSet],
				step=1/(overSet+1);
				data2=data[[i+1]];
				overl=overSet;
				,
				If[Length[overSet]==3,
					overl=overSet[[1]];
					step=1/(overl+1);
					drop1=overSet[[2]];
					drop2=overSet[[3]];
					If[drop1!=0,data1=Drop[data1,{-drop1,-1}]];
					If[drop2!=0,data2=Drop[data[[i+1]],{1,drop2}];,data2=data[[i+1]];];
					]
				]
			];
		set1=Drop[data1,{-overl,-1}];
		set1over=Take[data1,{-overl,-1}];
		set2=Drop[data2,{1,overl}];
		set2over=Take[data2,{1,overl}];
		(*joined=Joini2[{set1,set2},{set1over,set2over},step];*)
		joined=Joini[{set1,set2},{set1over,set2over},overl];
		];
		
	joined

	];
	



(* ::Subsubsection::Closed:: *)
(*Joini*)


Joini[sets_, setover_, step_] := Module[{over,dato,unit,noZero,tot},
  (*define the overlapping voxels*)
  unit = Unitize[setover];
  noZero = Times @@ unit;
  tot = Total[noZero];
  (*prepare the data for listable compliled function*)
  noZero = TransData[noZero, "l"];
  dato = TransData[TransData[setover, "l"], "l"];
  (*merge the overlapping data*)
  over = TransData[JoinFuncC[dato, noZero, tot, step], "r"];
  (*merge the non ovelap with the overlap*)
  Chop[Join[sets[[1]], over, sets[[2]]]]
  ]


JoinFuncC = Compile[{{dat, _Real, 2}, {noZero, _Integer, 1}, {tot, _Integer, 0}, {steps, _Integer, 0}},
	Block[{ran, unit, tot1, out},
    If[tot === 0,
     (*all zeros, no overlap of signals so just the sum of signals*)
     out = Total[(0 dat + 1) dat];
     ,
     (*overlap of signals*)
     (*define the range needed*)
     tot1 = 1./(tot + 1.);
     ran = Reverse@Range[tot1, 1. - tot1, tot1];
     
     (*replace with gradient*)
     If[tot === steps,
      (*full overlap*)
      out = Total[{ran, 1. - ran} dat];
      ,
      (*partial overlap*)
      (*summ all signals*)
      unit = (0 dat + 1);
      (*replace the overlapping signals with a gradient*)
      unit[[All, Flatten[Position[noZero, 1]]]] = {ran, 1 - ran};
      (*sum the signals*)
      out = Total[unit dat]
      ]];
    (*give the output*)
    out],
    {{out, _Real, 1}}, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsection::Closed:: *)
(*CorrectJoinSetMotion*)


Options[CorrectJoinSetMotion] = {JoinSetSplit -> True, PaddOverlap -> 2}

SyntaxInformation[CorrectJoinSetMotion] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

CorrectJoinSetMotion[input_, vox_, over_, OptionsPattern[]] := Module[
	{sets, nmax, dim, d1, d2, maskd1, maskd2, samp, overp, pad, regFunc, depth},
 	
 	(*get the input*)
 	pad = OptionValue[PaddOverlap];
 	depth = ArrayDepth[input];
 	
	(*data which will be joined, make all data sets 4D*)
	sets = Switch[depth,
		5,input,
		4,Transpose[{#}]&/@input
	];
	
	(*add z padding to allow more overlap*)
	sets = ArrayPad[#, Prepend[ConstantArray[{0, 0}, ArrayDepth[#] - 1], {pad, pad}]] & /@ sets;
	
	(*set needed values*)
	nmax = Length[sets];
	overp = over + 2 pad;
	dim = Dimensions[sets[[1,All,1]]];
	
	(*define the registration function*)
	regFunc = If[OptionValue[JoinSetSplit],RegisterDataTransformSplit,RegisterDataTransform];
	
	i=0;
	PrintTemporary[Dynamic[i]];
	
	(*perform the motion correction*)
	Table[
		i=n;
		(*get the seconds overlap stac*)
		d1 = sets[[n, ;; overp,1]];
		maskd1 = Dilation[#, 10] & /@ Mask[d1];
		(*pad to allow motion*)
		d1 = PadLeft[d1, dim];
		maskd1 = PadLeft[maskd1, dim];
		
		(*get the seconds overlap stac*)
		d2 = sets[[n + 1, -overp ;;,1]];
		maskd2 = Dilation[#, 10] & /@ Mask[d2];
		(*pad to allow motion*)
		d2 = PadLeft[d2, dim];
		maskd2 = PadLeft[maskd2, dim];
		
		maskd1 = maskd2 = Dilation[maskd1 maskd2, 1];
		
		(*get the number of samples for the registration*)
		samp = Round[((Total@Flatten@maskd1)+(Total@Flatten@maskd2))/20];
		
		(*perform the registration*)
		sets[[n + 1]] = Last@regFunc[{d1, maskd1, vox}, {d2, maskd2, vox}, {sets[[n + 1]], vox},
				MethodReg -> "translation", Iterations -> 250, NumberSamples -> samp, PrintTempDirectory -> False, InterpolationOrderReg -> 0];
		
		, {n, 1, nmax - 1}
	];
	
	(*output the data, make the 3D data 3D again*)
	Switch[depth,
		5,sets,
		4,sets[[All,All,1]]
		]
  
  ]


(* ::Subsection::Closed:: *)
(*TensMat*)


SyntaxInformation[TensMat] = {"ArgumentsPattern" -> {_}};

TensMat[tens:{_?ArrayQ..}]:=
Transpose[{{tens[[1]],tens[[4]],tens[[5]]},{tens[[4]],tens[[2]],tens[[6]]},{tens[[5]],tens[[6]],tens[[3]]}},{4,5,1,2,3}];

TensMat[tens_?ListQ]:=
	{{tens[[1]],tens[[4]],tens[[5]]},{tens[[4]],tens[[2]],tens[[6]]},{tens[[5]],tens[[6]],tens[[3]]}};


(* ::Subsection::Closed:: *)
(*TensVec*)


SyntaxInformation[TensVec] = {"ArgumentsPattern" -> {_}};

TensVec[tens:{{_,_,_},{_,_,_},{_,_,_}}]:=
	{tens[[1,1]],tens[[2,2]],tens[[3,3]],tens[[1,2]],tens[[1,3]],tens[[2,3]]};

TensVec[tens:{_?ArrayQ..}]:=
Transpose[Map[{#[[1,1]],#[[2,2]],#[[3,3]],#[[1,2]],#[[1,3]],#[[2,3]]}&,tens,{3}],{2,3,4,1}];


(* ::Subsection:: *)
(*Reshape data*)


(* ::Subsubsection::Closed:: *)
(*Data2DToVector*)


SyntaxInformation[Data2DToVector] = {"ArgumentsPattern" -> {_, _.}};

Data2DToVector[datai_]:=Data2DToVector[datai,1]
Data2DToVector[datai_,maski_]:=Block[{data,depth,mask,pos,vecdata,dimd,dimm},

depth=ArrayDepth[datai];
If[!(depth==2||depth==3),Return@Message[Data2DToVector::dim,depth]];
dimd=If[depth==2,Dimensions[datai],Drop[Dimensions[datai],1]];
dimm=Dimensions[maski];

If[!(mask===1)&&dimm!=dimd,Return@Message[Data2DToVector::mask,dimd,dimm]];

data=N@If[depth==3,maski #&/@datai,maski datai];
mask=Unitize@If[depth==3,Total[data],data];
pos=Position[mask,1];

vecdata=If[depth==3,
DeleteCases[Flatten[TransData[data,"l"],1],ConstantArray[0.,Length[datai]]],
DeleteCases[Flatten[data],0.]
];

{vecdata,{dimd,pos}}
];


(* ::Subsubsection::Closed:: *)
(*Data3DToVector*)


SyntaxInformation[Data3DToVector] = {"ArgumentsPattern" -> {_, _.}};

Data3DToVector[datai_]:=Data3DToVector[datai,1]
Data3DToVector[datai_,maski_]:=Module[{data,depth,mask,pos,vecdata,dimd,dimm},

depth=ArrayDepth[datai];
If[!(depth==3||depth==4),Message[Data3DToVector::dim,depth]];
dimd=If[depth==3,Dimensions[datai],Drop[Dimensions[datai],1]];
dimm=Dimensions[maski];

If[!(mask===1)&&dimm!=dimd,Return@Message[Data3DToVector::mask,dimd,dimm]];

data=N@If[depth==4,maski #&/@datai,maski datai];
mask=Unitize@If[depth==4,Total[data],data];
pos=Position[mask,1];

vecdata=If[depth==4,
DeleteCases[Flatten[TransData[data,"l"],2],ConstantArray[0.,Length[datai]]],
DeleteCases[Flatten[data],0.]
];

{vecdata,{dimd,pos}}
];


(* ::Subsubsection::Closed:: *)
(*VectorToData*)


SyntaxInformation[VectorToData] = {"ArgumentsPattern" -> {_, {_, _}}};

VectorToData[vec_,{dim_, pos_}]:=Block[{output,len},
len=Length@First@vec;
output=Switch[len,
0,ConstantArray[0.,dim],
_,ConstantArray[ConstantArray[0.,len],dim]
];
Switch[
Length[dim],
2,MapThread[(output[[#2[[1]],#2[[2]]]]=#1)&,{vec,pos}],
3,MapThread[(output[[#2[[1]],#2[[2]],#2[[3]]]]=#1)&,{vec,pos}]
];
Switch[len,0,output,_,TransData[output,"r"]]
]


(* ::Subsection::Closed:: *)
(*TransData*)


SyntaxInformation[TransData] = {"ArgumentsPattern" -> {_, _}};

TransData[data_,dir_]:=Block[{ran,dep,fun},
	ran=Range[dep=ArrayDepth[data]];
	fun=Switch[dir,"r",RotateLeft[ran],"l",RotateRight[ran]];
	Transpose[data,fun]
]


(* ::Subsection::Closed:: *)
(*DriftCorrect*)


Options[DriftCorrect]={NormalizeSignal->True,UseMask->True}

SyntaxInformation[DriftCorrect] = {"ArgumentsPattern" -> {_, _,_., OptionsPattern[]}}

DriftCorrect[data_, bi_, opts:OptionsPattern[]] := 
 Block[{bval,pos},

  bval = If[ArrayDepth[bi] == 2, BmatrixInv[bi][[1]], bi];
  
  pos = First@UniqueBvalPosition[bval, 6][[2]];

	DriftCorrect[data, bval, pos, opts]

  ];
  
DriftCorrect[data_, bi_,pos_, OptionsPattern[]] := 
 Block[{sig, cor, bval, sol1, sol2, sol3, a, b, c, x,outp},
  bval = If[ArrayDepth[bi] == 2, BmatrixInv[bi][[1]], bi];
  
  sig = MeanSignal[data, pos,UseMask->OptionValue[UseMask]];
  
  {sol1, sol2, sol3} = {a, b, c} /. FindFit[Transpose[{pos, sig}], {c + b x + a x^2}, {a, b, c}, x];
  cor = sol3/Table[sol3 + sol2 x + sol1 x^2, {x, 1, Length[bi]}];
  
  outp = ConstantArray[cor, Length[data]] data;
  
  If[OptionValue[NormalizeSignal], 100 outp / (sig[[1]] cor[[1]]) , outp ]
  ];


(* ::Subsection:: *)
(*Split and merge*)


(* ::Subsubsection::Closed:: *)
(*CutData*)


SyntaxInformation[CutData] = {"ArgumentsPattern" -> {_,_.}}

CutData[data_]:=CutData[data,FindMiddle[data]]

CutData[data_,cut_] := Switch[ArrayDepth[data],
		4,{data[[All, All, All, ;; cut]],data[[All, All, All, (cut + 1) ;;]],cut},
		3,{data[[All, All, ;; cut]], data[[All, All, (cut + 1) ;;]],cut}]

FindMiddle[dati_] := Module[{dat, fdat, len, datf,peaks,mid,peak,center,mask,ran},
  
	(*flatten mean and normalize data*)
	dat=dati;
	fdat = Flatten[dat];
	dat = Clip[dat, {0, Quantile[Pick[fdat, Unitize[fdat], 1], .95]}];
	dat = N@Nest[Mean, dat, ArrayDepth[dat] - 1];
	len = Length[dat];
	dat = len dat/Max[dat];
	mask = UnitStep[dat - .1 len];
	ran = Flatten[Position[mask, 1][[{1, -1}]]];
	
	(*smooth the data a bit*)
	datf = len - GaussianFilter[mask dat, len/20];
	(*find the peaks*)
	peaks = FindPeaks[datf];
	peaks = If[Length[peaks] >= 3, peaks[[2 ;; -2]], peaks];
	peaks = Select[peaks, (ran[[1]] < #[[1]] < ran[[2]]) &];
	
	(*find the most middle peak*)
	mid = Round[Length[dat]/2];
	center = {mid, .75 len};
	peak = Nearest[peaks, center];
	
	Print[Show[
	ListLinePlot[{len-dat,datf}, PlotStyle->{Black,Orange}],
	ListPlot[{peaks,peak,{center}},PlotStyle->(Directive[{PointSize[Large],#}]&/@{Blue,Red,Green})]
	,ImageSize->100]];
	
	(*output*)
	Round[First@First@peak]
  ]



(* ::Subsubsection::Closed:: *)
(*StichData*)


SyntaxInformation[StichData] = {"ArgumentsPattern" -> {_,_}}

StichData[datal_, datar_] := TransData[Join[TransData[datal, "r"], TransData[datar, "r"]], "l"];


(* ::Subsection:: *)
(*TransformData*)


(* ::Subsubsection::Closed:: *)
(*TransformData*)


Options[DataTranformation]={InterpolationOrder->1}

SyntaxInformation[DataTranformation]={"ArgumentsPattern"->{_,_,_,OptionsPattern[]}};

DataTranformation[data_, vox_, wi_,OptionsPattern[]] := 
 Block[{coor, rot, coorR, interFunc, interFuncC, w},
  w = If[Length[wi]==3,Join[wi,{0,0,0,1,1,1,0,0,0}]];
  
  coor = GetCoordinates[data, vox];
  rot = ParametersToTransformFull[w, "Inverse"];
  coorR = ApplyRotC[coor, rot];
  interFunc = 
   Interpolation[
    Transpose[{Flatten[coor, ArrayDepth[coor] - 2], Flatten[data]}], 
    InterpolationOrder -> OptionValue[InterpolationOrder], 
    "ExtrapolationHandler" -> {0. &, "WarningMessage" -> False}];
  interFuncC = Compile[{{coor, _Real, 1}}, interFunc[coor[[1]], coor[[2]], coor[[3]]], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];
  interFuncC[coorR]
  ]


(* ::Subsubsection::Closed:: *)
(*GetCoordinates*)


GetCoordinates[data_, vox_] := Block[{dim, off, coor},
   off = Dimensions[data]/2;
   coor = MapIndexed[#2 &, data, {ArrayDepth[data]}] - 0.5;
   CoordC[coor, off, vox]
   ];


(* ::Subsubsection::Closed:: *)
(*CoordC*)


   
CoordC = Compile[{{coor, _Real, 1}, {off, _Real, 1}, {vox, _Real, 1}},
    vox (coor - off),
   RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*ApplyRotC*)


ApplyRotC = Compile[{{coor, _Real, 1}, {rot, _Real, 2}}, 
  	(rot.Append[coor, 1])[[1 ;; 3]],
   RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*ParametersToTransformFull*)


ParametersToTransformFull[w_] := ParametersToTransformFull[w, "Normal"]

ParametersToTransformFull[w_, opt_] := Block[{
	tx, ty, tz, rx, ry, rz, sx, sy, sz, gx, gy, gz, 
	T, R, G, S, Rx, Ry, Rz, Gx, Gy, Gz, 
	mat, rMat, tMat}, 
	
	{rx, ry, rz, tx, ty, tz, sx, sy, sz, gx, gy, gz} = w;
	rx = -rx Degree; ry = -ry Degree; rz = -rz Degree;
	T = {{1, 0, 0, tx}, {0, 1, 0, ty}, {0, 0, 1, tz}, {0, 0, 0, 1}};
	
	Rx = {{1, 0, 0, 0}, {0, Cos[rx], Sin[rx], 0}, {0, -Sin[rx], Cos[rx], 0}, {0, 0, 0, 1}};
	Ry = {{Cos[ry], 0, -Sin[ry], 0}, {0, 1, 0, 0}, {Sin[ry], 0, Cos[ry], 0}, {0, 0, 0, 1}};
	Rz = {{Cos[rz], Sin[rz], 0, 0}, {-Sin[rz], Cos[rz], 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
	R = Rx.Ry.Rz;
	
	Gx = {{1, 0, gx, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
	Gy = {{1, 0, 0, 0}, {gy, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
	Gz = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, gz, 1, 0}, {0, 0, 0, 1}};
	G = Gx.Gy.Gz;
	
	S = {{sx, 0, 0, 0}, {0, sy, 0, 0}, {0, 0, sz, 0}, {0, 0, 0, 1}};
	
	mat = T.R.G.S;
	
	Switch[opt,
		"Normal",
		mat,
		"Inverse",
		rMat = Inverse[mat[[1 ;; 3, 1 ;; 3]]];
		tMat = -rMat.mat[[1 ;; 3, 4]];
		Append[Flatten /@ Thread[{rMat, tMat}], {0, 0, 0, 1}]
		]
	]


(* ::Subsubsection::Closed:: *)
(*ParametersToTransformFull*)


InvertDataset[data_] := Module[{dep},
  dep = ArrayDepth[data];
  Switch[dep,
   4, Transpose[Inverse3Di /@ Transpose[data]],
   _, Inverse3Di[data]
   ]
  ]

Inverse3Di[data_] := Block[{out},
  out = data;
  (out = Reverse[out, #]) & /@ {1, 2};
  out]



(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
