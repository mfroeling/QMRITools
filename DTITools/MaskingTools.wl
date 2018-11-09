(* ::Package:: *)

(* ::Title:: *)
(*DTITools MaskingTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`MaskingTools`", {"Developer`"}];

$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


Mask::usage =
"Mask[data] creates a mask by automatically finding a threshold.
Mask[data, min]creates a mask which selects only data above the min value.
Mask[data,{min,max}] creates a mask which selects data between the min and max value."

SmoothMask::usage = 
"SmoothMask[mask] generates one clean masked volume form a noisy mask."

SmartMask::usage = 
"SmartMask[input] crates a smart mask of input, which is either the tensor or the tensor parameters calculated using ParameterCalc.
SmartMask[input, mask] crates a smart mask of input and used the mask as a prior selection of the input."

MaskDTIdata::usage =
"MaskDTIdata[data, mask] aplies a mask to a DTI dataset."

MaskTensdata::usage =
"MaskTensdata[data, mask] aplies a mask a tensor."

GetMaskData::usage =
"GetMaskData[data, mask] retruns the data selected by the mask."

SplitSegmentations::usage = 
"SplitSegmentations[segmentation] splits a lable mask from ITKsnap or slicer3D in seperate masks and label numbers.
Output is masks and label numbers, {mask, labs}."

RescaleSegmentation::usage = 
"RescaleSegmentation[data, dim] rescales segmentations to given dimensions.
RescaleSegmentation[data, {vox1, vox2}] rescales segmentations from voxelsize vox1 to voxelsize vox2."

MergeSegmentations::usage = 
"MergeSegmentations[masks, labels] generates an ITKsnap or slices3D compatible segmentation from individual masks and label numbers.
Output is a labled segmentation."

SmoothSegmentation::usage =
"SmoothSegmentation[masks] smooths segmentations and removes the overlaps between multiple segmentations." 

RemoveMaskOverlaps::usage = 
"RemoveMaskOverlaps[mask] removes the overlaps between multiple masks. Mask is a 4D dataset with {z, masks, x, y}"

SegmentMask::usage = 
"SegmentMask[mask, n] devides a mask in n equal segments along the slice direction. n must be an integer."

NormalizeData::usage = 
"NormalizeData[data] normalizes the data to the mean signal of the data.
NormalizeData[data,{min,max}] normalizes the data between min and max."

NormalizeDiffData::usage = 
"NormalizeDiffData[data] normalizes the diffusion data to the mean signal of the first volume."

HomoginizeData::usage = 
"HomoginizeData[data, mask] tries to homoginize the data within the mask by removing intensity gradients."

ROIMask::usage = 
"ROIMask[maskdim, {name->{{{x,y},slice}..}..}] crates mask from coordinates x and y at slice. 
maskdim is the dimensions of the output {zout,xout,yout}."


(* ::Subsection:: *)
(*Options*)


MaskSmoothing::usage = 
"MaskSmoothing is an options for Mask, if set to True it smooths the mask, by closing holse and smoothing the contours."

MaskComponents::usage =
"MaskComponents is an option for Mask and SmoothMask. Determinse the amount of largest clusters used as mask." 

MaskClosing::usage =
"MaskClosing  is an option for Mask and SmoothMask. The size of the holes in the mask that will be closed" 

MaskFiltKernel::usage =
"MaskFiltKernel is an option for Mask, SmoothMask and SmoothSegmentation. How mucht the contours are smoothed." 

Strictness::usage = 
"Strictness is an option for SmartMask value between 0 and 1. Higer values removes more data."

MaskCompartment::usage = 
"MaskCompartment is an option for SmartMask. Can be \"Muscle\" or \"Fat\"."

SmartMethod::usage = 
"SmartMethod is an option for SmartMask. This specifies how the mask is generated. Can be \"Continuous\" or \"Catagorical\""

SmartMaskOutput::usage = 
"SmartMaskOutput is an option for Smartmask. Can be set to \"mask\" to output only the mask or \"full\" to also output the probability mask."

MeanOutput::usage = 
"MeanOutput is an option for NormalizeDiffData. If True it will also output the normalization factor."

GetMaskOutput::usage = 
"GetMaskOutput is an option for GetMaskData. Defaul is \"Slices\" which gives the mask data per slices. Else the entire mask data is given as output."


(* ::Subsection::Closed:: *)
(*Error Messages*)


Mask::tresh = "Given treshhold `1` value is not a vallid input, must be a number for min treshhold only or a vector {min tresh, max tresh}."

GetMaskData::tresh = "The dimensions of the data and the mask must be the same, dataset: `1`, mask: `2`."

ROIMask::war = "there are more slices in the roi set than in the given dimensions"


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*Mask*)


Options[Mask]={MaskSmoothing -> False, MaskComponents -> 1, MaskClosing -> 5, MaskFiltKernel -> 2};

SyntaxInformation[Mask] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

Mask[data_,opts:OptionsPattern[]]:=Mask[data, 0, opts]

Mask[data_?ArrayQ,tr_,opts:OptionsPattern[]]:=
Module[{mask,tresh},
	If[Length[tresh]>2,
		Message[Mask::tresh,tresh]
		,
		mask=If[tr===0,
			(*no threshhold*)
			ImageData[Binarize[Image3D[NormalizeData[data,.95]]]]
			,
			(*threshhold*)
			tresh=If[NumberQ[tr],{tr},tr];		
			
			If[Length[tresh]==1,
				UnitStep[data-tresh[[1]]],
				UnitStep[data-tresh[[1]]]-UnitStep[data-tresh[[2]]]
				]
			];
				
		If[OptionValue[MaskSmoothing],
			SmoothMask[mask,FilterRules[{opts},Options[SmoothMask]]],
			mask
			]
		]
	]


(* ::Subsection::Closed:: *)
(*SmoothMask*)


Options[SmoothMask]={MaskComponents->1, MaskClosing->5, MaskFiltKernel->2}

SyntaxInformation[SmoothMask] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

SmoothMask[mask_,OptionsPattern[]] := Block[{pad, close,obj,filt},
  close = OptionValue[MaskClosing];(*close holes in mask*)
  pad = 3 close;
  obj = OptionValue[MaskComponents];(*number of objects that are maintained*)
  filt = OptionValue[MaskFiltKernel];(*how much smooting*)

  Round[GaussianFilter[ArrayPad[Closing[ImageData[SelectComponents[Image3D[ArrayPad[mask, pad]],"Count", -obj]], close],-pad], filt]]
  ]


(* ::Subsection::Closed:: *)
(*SmartMask*)


Options[SmartMask]={Strictness->0.50, MaskCompartment->"Muscle", SmartMethod->"Continuous", SmartMaskOutput->"mask"};

SyntaxInformation[SmartMask] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

SmartMask[input_,ops:OptionsPattern[]]:=SmartMask[input, 0, ops]

SmartMask[input_,maski_,OptionsPattern[]]:=Module[{
	sol,func,range,map,mask,pmask,pars
	},
	
	(*get the parameter from the tensor else use input parameters*)
	pars = If[Length[input]==6,
		PrintTemporary["Caculating Parameters"];
		ParameterCalc[input],
		input];
	
	pmask = Mask[pars[[4]] , {0.1, 4}];
	
	(*find the histogram solution*)
	sol=If[maski===0,
		Switch[
			OptionValue[MaskCompartment],
			"Muscle",
			ParameterFit2[pars][[All,{3,5,7}]],
			"Fat",
			ParameterFit2[pars][[All,{2,4,6}]]
			]
			,
			ParameterFit[GetMaskData[#,maski pmask]&/@pars,FitOutput->"BestFitParameters"]
		];
	
	
	Switch[OptionValue[SmartMethod],
		"Catagorical",
		range = (func = SkewNormalDistribution[#2[[1]], #2[[2]], #2[[3]]]; Quantile[func, {.02, .98}]) & /@ sol;
		range = Clip[range, {0, Infinity}, {10^-3, 0.}];
		map = Total[MapThread[Mask[#1,#2]&,{pars,range}]]/5;
		mask = pmask * Mask[TotalVariationFilter[map,.15],{OptionValue[Strictness]}];
		,
		"Continuous",
		map = MapThread[PDF[SkewNormalDistribution[#2[[1]], #2[[2]], #2[[3]]], #1] &, {pars, sol}];
		map = Total[{1, 1, 1, 1, 2}*(#/Max[#] & /@ map)]/6;
		mask = pmask * Mask[TotalVariationFilter[map, .25], {OptionValue[Strictness]}];
		];
		
		If[OptionValue[SmartMaskOutput]==="mask",mask,{mask,map}]
	]


(* ::Subsection:: *)
(*Apply Masks*)


(* ::Subsubsection::Closed:: *)
(*MaskDTIdata*)


SyntaxInformation[MaskDTIdata] = {"ArgumentsPattern" -> {_, _}};

MaskDTIdata[data_, mask_] := Transpose[mask # & /@ Transpose[data]]


(* ::Subsubsection::Closed:: *)
(*MaskTensdata*)


SyntaxInformation[MaskTensdata] = {"ArgumentsPattern" -> {_, _}};

MaskTensdata[tens_, mask_] := mask # & /@ tens


(* ::Subsection::Closed:: *)
(*GetMaskData*)


Options[GetMaskData] = {GetMaskOutput -> "All"}

SyntaxInformation[GetMaskData] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

GetMaskData[data_?ArrayQ, mask_?ArrayQ, OptionsPattern[]] := Module[{depth,fdat},
	If[Dimensions[data]!=Dimensions[mask],
		Message[GetMaskData::dim,Dimensions[data],Dimensions[mask]]
		,
		Switch[OptionValue[GetMaskOutput],
			"Slices",
			depth=ArrayDepth[data];
			Map[(
				fdat=Flatten[N@#];
				Pick[fdat, Unitize[fdat], 1]
				)&,	data*mask, {depth-2}],
			_
			,
			fdat=N@Chop[Flatten[N[data mask]]];
			Pick[fdat, Unitize[fdat], 1]
		]
	]
]


(* ::Subsection:: *)
(*Segmentation functions*)


(* ::Subsubsection::Closed:: *)
(*SplitSegmentations*)


SyntaxInformation[SplitSegmentations] = {"ArgumentsPattern" -> {_}};

SplitSegmentations[masksI_] := Block[{vals, masks},
	masks=SparseArray[masksI];
	vals = DeleteCases[Sort@Round[DeleteDuplicates[Flatten[masksI]]],0];
	masks =(1 - Unitize[masks - #]) & /@ vals;
	masks=Normal[Transpose[masks]];
	{masks, vals}
  ]


(* ::Subsubsection::Closed:: *)
(*MergeSegmentations*)


SyntaxInformation[MergeSegmentations] = {"ArgumentsPattern" -> {_,_}};

MergeSegmentations[masks_, vals_] := Total[vals Transpose@masks];


(* ::Subsubsection::Closed:: *)
(*SmoothSegmentation*)


Options[SmoothSegmentation] = {MaskFiltKernel -> 2}

SyntaxInformation[SmoothSegmentation] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

SmoothSegmentation[masks_, OptionsPattern[]] := 
 Block[{maskInp, maskOver, maskOut, smooth, posOver, x, y, z, p},
  
	smooth = OptionValue[MaskFiltKernel];
	(*convert data to sparse Array and transpose*)
	maskInp = Transpose[SparseArray[masks]];
	(*Get smoothed or non smoothed masks*)
	maskInp = SparseArray[If[smooth === False, maskInp, SparseArray[Round[GaussianFilter[#, smooth] + 0.15]] & /@ maskInp]];
	(*find the overlaps*)
	maskOver = Mask[Total[maskInp], 1.5];
	posOver = Position[maskOver, 1, 3];
	(*get smoothed masks without overlap*)
	maskOut = SparseArray[(1 - maskOver) # & /@ maskInp];
	
	maskInp = Normal[maskInp];
	maskOut = Normal[maskOut];
	(
	    (*find values to fill up the overlaps, 
	    always pick the first value*)
	    {z, x, y} = #;
	    p = First[FirstPosition[maskInp[[All, z, x, y]], 1]];
	    maskOut[[p, z, x, y]] = 1;
	    ) & /@ posOver;
	Clear[maskInp, maskOver];
	Normal[Transpose[SparseArray[maskOut]]]
  ]


(* ::Subsubsection::Closed:: *)
(*RemoveMaskOverlaps*)


SyntaxInformation[RemoveMaskOverlaps] = {"ArgumentsPattern" -> {_}};
	
RemoveMaskOverlaps[masks_] := SmoothSegmentation[masks, MaskFiltKernel->False];


(* ::Subsubsection::Closed:: *)
(*RemoveMaskOverlaps*)


SyntaxInformation[RescaleSegmentation] = {"ArgumentsPattern" -> {_}};

RescaleSegmentation[seg_, vox_] := Block[{segs, val},
  If[ArrayDepth[seg] == 3, {segs, val} = SplitSegmentations[seg], 
   segs = seg];
  segs = RemoveMaskOverlaps[
    Transpose[
     Round[RescaleData[#, vox, InterpolationOrder -> 1]] & /@ 
      Transpose[segs]]];
  If[ArrayDepth[seg] == 3, MergeSegmentations[segs, val], segs]
  ]

(* ::Subsubsection::Closed:: *)
(*SegmentMask*)


SyntaxInformation[SegmentMask] = {"ArgumentsPattern" -> {_, _, _.}};

SegmentMask[mask_, seg_?IntegerQ] := Block[{pos, f, l, sel, out},
  pos = Flatten@Position[Unitize[Total[Flatten[#]]] & /@ mask, 1];
  {f, l} = {First[pos], Last[pos]};
  sel = Partition[Round[Range[f, l, ((l - f)/seg)]], 2, 1] + Append[ConstantArray[{0, -1}, seg - 1], {0, 0}];
  out = ConstantArray[0*mask, seg];
  Table[out[[i, sel[[i, 1]] ;; sel[[i, 2]]]] = mask[[sel[[i, 1]] ;; sel[[i, 2]]]], {i, 1, seg}];
  out
  ]


(* ::Subsection:: *)
(*Normalize Functions*)


(* ::Subsubsection::Closed:: *)
(*NormalizeData*)


SyntaxInformation[NormalizeData] = {"ArgumentsPattern" -> {_,_.}};

NormalizeData[data_] := Switch[ArrayDepth[data],
	3, NormalizeData[data, Mask[data]],
	4, NormalizeData[data, Mask[Mean[Transpose[data]]]]
	]

NormalizeData[data_, mask_] := Block[{mn},
	mn = Switch[ArrayDepth[data],
	3, N[MeanNoZero[Flatten[mask data]]/100.],
	4, N[MeanNoZero[Flatten[mask data[[All, 1]]]]/100.]
	];
	data / mn
]


(* ::Subsubsection::Closed:: *)
(*NormalizeDiffData*)


Options[NormalizeDiffData] = {MeanOutput->False}

SyntaxInformation[NormalizeDiffData] = {"ArgumentsPattern" -> {_, _.}};

NormalizeDiffData[data_,ops:OptionsPattern[]] := NormalizeDiffData[data, Mask[Mean[Transpose[data]], MaskSmoothing -> True],ops]
NormalizeDiffData[data_, mask_,OptionsPattern[]] := Block[{mn,dataout},
	mn = N[MeanNoZero[Flatten[mask data[[All, 1]]]]/100.];
	dataout = data/mn;
	If[OptionValue[MeanOutput],{dataout,mn},dataout]
]


(* ::Subsubsection::Closed:: *)
(*HomoginizeData*)


SyntaxInformation[HomoginizeData] = {"ArgumentsPattern" -> {_, _}};

HomoginizeData[datai_, mask_] := 
 Module[{data, mn, fit, datac, maskout},
  data = mask GaussianFilter[datai, 5];
  mn = Mean[Cases[Flatten[N[data]],Except[0.]]];
  fit = FitGradientMap[Erosion[mask, 3] data];
  
  datac = (datai/(fit + 0.001));
  maskout = Mask[datac, 0.1];
  maskout = Dilation[SmoothMask[maskout, 5], 5];
  maskout Clip[mn datac, {0.8, 1.5} MinMax[data], {0, 0}]
  ]
  
  
FitGradientMap[data_] := Module[{func, x, y, z, coor},
  Clear[x, y, z];
  coor = Flatten[
    MapIndexed[
     ReleaseHold@If[#1 == 0, Hold[Sequence[]], Join[#2, {#1}]] &, 
     data, {3}], 2];
  func = Fit[coor, {1, x, y, z, x y, x z, y z, z^2, x^2, y^2(*, z^3, x^3, y^3, z x x, y x x, z y y , x y y , y z z, x z z*)}, {x, y, z}];
  {x, y, z} = TransData[Array[{#1, #2, #3} &, Dimensions[data]], "r"];
  func
  ]


(* ::Subsection::Closed:: *)
(*ROIMask*)


SyntaxInformation[ROIMask] = {"ArgumentsPattern" -> {_, _, _.}};

ROIMask[ROIdim_,maskdim_,ROI:{(_?StringQ->{{{{_?NumberQ,_?NumberQ}..},_?NumberQ}..})..}]:=
Module[{output},
	output=Map[#[[1]]->ROIMask[ROIdim,maskdim,#[[2]]]&,ROI];
	Print["The Folowing masks were Created: ",output[[All,1]]];
	Return[output]
	]

ROIMask[ROIdim_,maskdim_,ROI:{{_?StringQ->{{{{_?NumberQ,_?NumberQ}..},_?NumberQ}..}}..}]:=
Module[{output},
	output=Map[#[[1,1]]->ROIMask[ROIdim,maskdim,#[[1,2]]]&,ROI];
	Print["The Folowing masks were Created: ",output[[All,1]]];
	Return[output]
	]

ROIMask[ROIdim_,maskdim_,ROI:{{{{_?NumberQ,_?NumberQ}..},_?NumberQ}..}]:=
Module[{output,ROIcor,ROIslice,msk},
	output=ConstantArray[0,Join[{ROIdim[[1]]},maskdim]];
	If[ROI[[All,1]]!={{{0,0}}},
		ROIcor=Round[ROI[[All,1]]];
		ROIslice=Clip[ROI[[All,2]],{1,ROIdim[[1]]}];
		msk=1-ImageData[Image[Graphics[Polygon[#],PlotRange->{{0,ROIdim[[3]]},{0,ROIdim[[2]]}}],"Bit",ColorSpace->"Grayscale",ImageSize->maskdim]]&/@ROIcor;
		MapIndexed[output[[#1]]=msk[[First[#2]]];&,ROIslice];
		];
	Return[output];
	]

ROIMask[maskdim_,ROI:{(_?StringQ->{{{{_?NumberQ,_?NumberQ}..},_?NumberQ}..})..}]:=
Module[{output},
	output=Map[#[[1]]->ROIMask[maskdim,#[[2]]]&,ROI];
	Print["The Folowing masks were Created: ",output[[All,1]]];
	Return[output]
	]

ROIMask[maskdim_,ROI:{{_?StringQ->{{{{_?NumberQ,_?NumberQ}..},_?NumberQ}..}}..}]:=
Module[{output},
	output=Map[#[[1,1]]->ROIMask[maskdim,#[[1,2]]]&,ROI];
	Print["The Folowing masks were Created: ",output[[All,1]]];
	Return[output]
	]

ROIMask[maskdim_,ROI:{{{{_?NumberQ,_?NumberQ}..},_?NumberQ}..}]:=
Module[{output, ROIcor, ROIslice, msk},
 output = ConstantArray[0, maskdim];
 If[ROI[[All, 1]] != {{{0, 0}}},
  ROIcor = Round[Map[Reverse[maskdim[[2 ;; 3]]]*# &, ROI[[All, 1]], {2}]];
  If[Max[ROI[[All, 2]]] > maskdim[[1]], Message[ROIMask::war]];
  ROIslice = Clip[ROI[[All, 2]], {1, maskdim[[1]]}];
  msk = 1 - 
      ImageData[
       Image[Graphics[Polygon[#], 
         PlotRange -> {{0, maskdim[[3]]}, {0, maskdim[[2]]}}], "Bit", 
        ColorSpace -> "Grayscale", 
        ImageSize -> maskdim[[2 ;; 3]]]] & /@ ROIcor;
  MapIndexed[output[[#1]] = msk[[First[#2]]]; &, ROIslice];
  ];
 Return[output];]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
