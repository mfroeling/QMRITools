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


(* ::Subsection::Closed:: *)
(*Functions*)


Mask::usage =
"Mask[data,tresh min]creates a mask which selects only data above the treshlow value.
Mask[data,{tresh min,tresh max}] creates a mask which selects data between the tresh min and tresh max value."

SmartMask::usage = 
"SmartMask[input, mask] crates a smart mask of input based on mask."

SmartMask2::usage = 
"SmartMask2[input] crates a smart mask of input based.
SmartMask2[input, mask] crates a smart mask of input based on mask."

GetMaskData::usage =
"GetMaskData[data, mask] retruns the data selected by the mask as a single list for each image."

ROIMask::usage = 
"ROIMask[maskdim, {name->{{{x,y},slice}..}..}] crates mask from coordinates x and y at slice. \
maskdim is the dimensions of the output {zout,xout,yout}."

GofImport::usage = 
"GofImport[file] imports *.gof file to mathematica expresssion" 

ReadGof::usage = 
"ReadGof[file.gof, T1.dcm] imports the gof file to a format that can be used in ROIMask.
ReadGof[{file1.gof, file2.gof, ..}, T1.dcm] imports the gof files to a format that can be used in ROIMask.
ReadGof[{file1.gof, file2.gof, ..}, {T1-1.dcm, T1-2.dcm, ..}] imports the gof files to a format that can be used in ROIMask,\
where each .gof file correstponds to a different T1 file."

ReadROI::usage = 
"ReadROI[file, voxel, dim] imports *.gof file to format that can be used for ROImask (better use ReadGof).
ReadROI[filename, voxel, dim, off] imports *.gof file to format that can be used for ROImask using offset off (imagePosition dicom header)."

MaskDTIdata::usage =
"MaskDTIdata[data, mask] aplies a mask to a DTI dataset."

MaskTensdata::usage =
"MaskTensdata[data, mask] aplies a mask a tensor."

SmoothMask::usage = 
"SmoothMask[mask] generates one clean masked volume form a noisy mask.
SmoothMask[mask, int] higher number of int increases the maks size, default value is 2."

RemoveMaskOverlaps::usage = 
"RemoveMaskOverlaps[mask] removes the overlaps between multiple masks:"

SmoothSegmentation::usage =
"SmoothSegmentation[masks] smooths segmentations and removes the overlaps between multiple masks." 

SplitSegmentations::usage = 
"SplitSegmentations[segmentation] splits a lable mask from ITKsnap or slicer3D in seperate masks and label numbers.

Output is masks and label numbers."

MergeSegmentations::usage = 
"MergeSegmentations[masks,labels] generates an ITKsnap or slices3D compatible segmentation from individual masks and label numbers.

Output is a labled segmentation."

HomoginizeData::usage = 
"HomoginizeData[data, mask] tries to homoginize the data within the mask by removing intensity gradients."

NormalizeData::usage = 
"NormalizeData[data] normalizes the data to the mean signal of the data.
NormalizeData[data,{min,max}] normalizes the data between min and max."

NormalizeDiffData::usage = 
"NormalizeDiffData[data] normalizes the diffusion data to the mean signal of the first volume."

SegmentMask::usage = 
"SegmentMask[mask, n] devides a mask in n equal segments along the slice direction. n must be an integer."


(* ::Subsection::Closed:: *)
(*Options*)


Smoothing::usage = 
"Smoothing is an options for Mask, Maskbin and SmartMask functions, if set to true (default) it smooths (removes holes and smooth edges) the mask"

DataType::usage = 
"DataType is an option for ReadGof"

Strictness::usage = 
"Strictness is an option for SmartMask (value of 1 to 6) and SmartMask2 (value between 0 and 1). Low values selects more."

Compartment::usage = 
"Compartment is an option for SmartMask2. Can be \"Muscle\" or \"Fat\"."

SmoothMaskFactor::usage = 
"SmoothMaskFactor is an option for SmartMask."

OptimizationRuns::usage = 
"OptimizationRuns is an option for SmartMask"

MaskRange::usage = 
"MaskRange is an option for SmartMask"

MaskComponents::usage =
"MaskComponents is an option for SmoothMask. Determinse the amount of largest clusters used as mask." 

MaskPadding::usage =
"MaskPadding is an option for SmoothMask. Prevents the mask merging with the edge." 

MaskClosing::usage =
"MaskClosing  is an option for SmoothMask. The size of the holes in the mask that will be closed" 

MaskFiltKernel::usage =
"MaskFiltKernel is an option for SmoothMask and SmoothSegmentation. How mucht the contours are smoothed." 

MeanOutput::usage = 
"MeanOutput is an option for NormalizeDiffData. If True it will also output the normalization factor."

GetMaskOutput::usage = 
"GetMaskOutput is an option for GetMaskData. Defaul is \"Slices\" which gives the mask data per slices. Else the entire mask data is given as output."


(* ::Subsection::Closed:: *)
(*Error Messages*)


Mask::tresh = "Given treshhold `1` value is not a vallid input, must be a number for min treshhold only or a vector {min tresh, max tresh}."

GetMaskData::tresh = "The dimensions of the data and the mask must be the same, dataset: `1`, mask: `2`."

ReadGof::dat = "The given DataType is not valid: `1`. Can be \"Arm\" or \"Mara\"."

ReadGof::file= "The given file does not exist: `1`."

GofImport::file= "The given file does not exist: `1`."

ROIMask::war = "there are more slices in the roi set than in the given dimensions"


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection:: *)
(*Mask from gof files*)


(* ::Subsubsection::Closed:: *)
(*ReadROI*)


SyntaxInformation[ReadROI] = {"ArgumentsPattern" -> {_, _, _, _.}};

ReadROI[filename_,voxel_,dim_]:=
Module[{data,output,roi,names,off},
	off=("ImagePosition"/. Import["T1\\00001.dcm","MetaInformation"]);
	{names,roi}=GofImport[filename];
	data=Round[Map[(#-off)/voxel&,roi,{3}]];
	output=Map[{Transpose[{#[[All,1]]+1,-(#[[All,2]]-dim[[2]])}],#[[1,3]]+1}&,data,{2}];
	Return[MapThread[{#1->#2}&,{names,output}]]
]

ReadROI[filename_,voxel_,dim_,off_]:=
Module[{data,output,roi,names},
	{names,roi}=GofImport[filename];
	data=Round[Map[(#-off)/voxel&,roi,{3}]];
	output=Map[{Transpose[{#[[All,1]]+1,-(#[[All,2]]-dim[[2]])}],#[[1,3]]+1}&,data,{2}];
	Return[MapThread[{#1->#2}&,{names,output}]]
	]


(* ::Subsubsection::Closed:: *)
(*ReadGof*)


Options[ReadGof]={DataType->"Normal"};

SyntaxInformation[ReadGof] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};

ReadGof[file:{_?StringQ..}, T1:{_?StringQ..},OptionsPattern[]]:= MapThread[ReadGof[#1,#2,DataType->OptionValue[DataType]] &, {file,T1}]

ReadGof[file:{_?StringQ..}, T1_?StringQ,OptionsPattern[] ]:= Map[ReadGof[#,T1,DataType->OptionValue[DataType]] &, file]

ReadGof[file_?StringQ, T1_?StringQ,OptionsPattern[]]:= Module[{slice, pixel, output, off, roifile, T1file, names, roi, data, meta, vox, dim},
	roifile=If[StringTake[file,-4]==".gof",file,file<>".gof"];
	(*T1file=If[StringTake[T1,-4]==".dcm",T1,T1<>".dcm"];*)
	T1file=T1;
	If[FileExistsQ[T1file],
		If[FileExistsQ[roifile],
			{names, roi} = GofImport[roifile], 
			Message[ReadGof::file, roifile]
		];
		meta = Import[T1file,"MetaInformation"];
		slice = "SlicesSpacing" /. meta;
		
		pixel = "PixelSpacing" /. meta;
		If[StringQ[pixel], pixel = "PixelSpacing" /.("(0028,9110)" /.First["(5200,9230)" /. meta])];
		If[OptionValue[DataType] == "Maastricht", pixel = {1, 1}];
		
		vox = Flatten[{slice,pixel}/.meta];
		
		off = "ImagePosition" /. meta;
		off = If[StringQ[off], "ImagePosition" /.First["(0020,9113)" /. First["(5200,9230)" /. meta]],off];
		
		dim = Dimensions[Import[T1file,"Data"]];
		,Message[ReadGof::file, T1file]
	];
	
	Switch[
		OptionValue[DataType]
		,
		"Normal"
		,
		data = Round[Map[(# - off)/{vox[[2]], vox[[3]], vox[[1]]} &, roi, {3}]];
		output = Map[{Transpose[{#[[All, 1]] + 1, -(#[[All, 2]] - dim[[3]])}/{dim[[3]],dim[[2]]}], #[[1, 3]] + 1} & , data, {2}];
		,
		"Maastricht",
		data = Round[Map[#/{vox[[2]], vox[[3]], vox[[1]]} &, roi, {3}]];
		output = Map[{Transpose[{#[[All, 1]] + 1, -(#[[All, 2]] - dim[[3]])}/{dim[[3]], dim[[2]]}], #[[1, 3]] + 1} &, data, {2}];
		,
		"Arm"
		,
		(*off = off - {0, 0, dim[[2]]*vox[[2]]};*)
		data = Round[Map[(# - {off[[1]],off[[3]],off[[2]]})/vox &, roi, {3}]];
		output = Map[{Transpose[{#[[All, 2]] + 1, #[[All, 3]] + 1}/{dim[[3]], dim[[2]]}], -(#[[1, 1]] - dim[[1]])} &, data, {2}];
	];
		
	Return[Thread[names -> output//N]];
];


(* ::Subsubsection::Closed:: *)
(*GofImport*)


SyntaxInformation[GofImport] = {"ArgumentsPattern" -> {_}};

GofImport[file_]:=
Module[{roiData, roiPos, roiSeg, roiNames, roi},
  If[! FileExistsQ[file],
   Message[GofImport::file, file],
   roiData = Import[file, "Lines"];
   roiPos = # - {0, 1} & /@ 
     Partition[
      Append[Flatten[
        Drop[#, -1] & /@ 
         Position[StringPosition[roiData, "ConcaveHull"], {_, _}]], 
       Length[roiData] + 1], 2, 1];
   roiSeg = Take[roiData, #] & /@ roiPos;
   roiNames = 
    StringReplace[
     StringTake[#[[1]], 
        DeleteDuplicates[Flatten[StringPosition[#[[1]], "\""]]]] & /@ 
      roiSeg, "\"" -> ""];
   roi = 1000*
     ToExpression[
      Fold[StringReplace, #, {{") }" -> ")}", "{ (" -> "{{", 
           " " .. -> " ","E"-> "*10^"}, {"( " -> "{", " )" -> "}", "(" -> "{", 
           ")" -> "}"}, {" " -> ","}}] & /@ (Map[
         "{" <> StringTake[#, {StringPosition[#, "("][[1, 1]], 
             StringPosition[#, "}"][[1, 1]]}] &, (StringReplace[#, 
             "}}" -> "}"] & /@ (Drop[#, 1] & /@ roiSeg)), {2}])];
   {roiNames, roi}
   ]
  ]


(* ::Subsubsection::Closed:: *)
(*ROIMask*)


SyntaxInformation[SegmentMask] = {"ArgumentsPattern" -> {_, _, _.}};

SegmentMask[mask_, seg_?IntegerQ] := Block[{pos, f, l, sel, out},
  pos = Flatten@Position[Unitize[Total[Flatten[#]]] & /@ mask, 1];
  {f, l} = {First[pos], Last[pos]};
  sel = Partition[Round[Range[f, l, ((l - f)/seg)]], 2, 1] + 
    Append[ConstantArray[{0, -1}, seg - 1], {0, 0}];
  out = ConstantArray[0*mask, seg];
  Table[out[[i, sel[[i, 1]] ;; sel[[i, 2]]]] = 
    mask[[sel[[i, 1]] ;; sel[[i, 2]]]], {i, 1, seg}];
  out
  ]


(* ::Subsubsection::Closed:: *)
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
 
 


(* ::Subsection::Closed:: *)
(*Mask*)


Options[Mask]={Smoothing -> False, MaskComponents -> 1, MaskPadding -> 40, MaskClosing -> 20, MaskFiltKernel -> 2};

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
				
		If[OptionValue[Smoothing],
			SmoothMask[mask,FilterRules[{opts},Options[SmoothMask]]],
			mask
			]
		]
	]


(* ::Subsection::Closed:: *)
(*GetMaskData*)

Options[GetMaskData] = {GetMaskOutput -> "Slice"}

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
(*SmartMask*)


Options[SmartMask]={Strictness->.50, Compartment->"Muscle", Method->"Continuous", Reject->True, Output->"mask"};

SyntaxInformation[SmartMask] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

SmartMask[input_,ops:OptionsPattern[]]:=SmartMask[input, 0, ops]

SmartMask[input_,maski_,OptionsPattern[]]:=Module[{
	sol,func,range,map,mask,pmask,pars
	},
	
	(*get the parameter from the tensor else use input parameters*)
	pars = If[Length[input]==6,
		PrintTemporary["Caculating Parameters"];
		ParameterCalc[input,Reject->OptionValue[Reject]],
		input];
	
	pmask = Mask[pars[[4]] , {0.1, 4}];
	
	(*find the histogram solution*)
	sol=If[maski===0,
		Switch[
			OptionValue[Compartment],
			"Muscle",
			ParameterFit2[pars][[All,{3,5,7}]],
			"Fat",
			ParameterFit2[pars][[All,{2,4,6}]]
			]
			,
			ParameterFit[GetMaskData[#,maski pmask]&/@pars,FitOutput->"BestFitParameters"]
		];
	
	
	Switch[OptionValue[Method],
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
		
		If[OptionValue[Output]==="mask",mask,{mask,map}]
	]


(* ::Subsection::Closed:: *)
(*MaskDTIdata*)


SyntaxInformation[MaskDTIdata] = {"ArgumentsPattern" -> {_, _}};

MaskDTIdata[data_, mask_] := Transpose[mask # & /@ Transpose[data]]


(* ::Subsection::Closed:: *)
(*MaskTensdata*)


SyntaxInformation[MaskTensdata] = {"ArgumentsPattern" -> {_, _}};

MaskTensdata[tens_, mask_] := mask # & /@ tens


(* ::Subsection::Closed:: *)
(*SmoothMask*)


Options[SmoothMask]={MaskComponents->1,MaskPadding->40,MaskClosing->20, MaskFiltKernel->2}

SyntaxInformation[SmoothMask] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

SmoothMask[mask_,OptionsPattern[]] := Block[{pad, close,obj,filt},
  pad = Clip[OptionValue[MaskPadding],{4,100}];(*to prevent mask joining with edges*)
  close = Clip[OptionValue[MaskClosing],{1,pad-2}];(*close holes in mask*)
  obj = OptionValue[MaskComponents];(*number of objects that are maintained*)
  filt = OptionValue[MaskFiltKernel];(*how much smooting*)

  Round[GaussianFilter[ArrayPad[Closing[ImageData[SelectComponents[Image3D[ArrayPad[mask, pad]],"Count", -obj]], close],-pad], filt]]
  ]


(* ::Subsection:: *)
(*Segmentation functions*)


(* ::Subsubsection::Closed:: *)
(*RemoveMaskOverlaps*)


SyntaxInformation[RemoveMaskOverlaps] = {"ArgumentsPattern" -> {_}};
	
RemoveMaskOverlaps[masks_] := SmoothSegmentation[masks, False];


(* ::Subsubsection::Closed:: *)
(*SmoothSegmentation*)


Options[SmoothSegmentation] = {MaskFiltKernel -> 2}

SyntaxInformation[SmoothSegmentation] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

SmoothSegmentation[masks_, OptionsPattern[]] := 
 Block[{maskInp, maskOver, maskOut, posOver,smooth, x, y, z, p},
  
	smooth=OptionValue[MaskFiltKernel];
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
(*SplitSegmentations*)


SyntaxInformation[SplitSegmentations] = {"ArgumentsPattern" -> {_}};

SplitSegmentations[masksI_] := Block[{vals, masks},
	masks=SparseArray[masksI];
	vals = Sort@DeleteDuplicates[Flatten[masksI]][[2 ;;]];
	masks = Mask[masks, {# - .5, # + .5}] & /@ vals;
	masks=Normal[Transpose[masks]];
	{masks, vals}
  ]


(* ::Subsubsection::Closed:: *)
(*MergeSegmentations*)


SyntaxInformation[MergeSegmentations] = {"ArgumentsPattern" -> {_,_}};

MergeSegmentations[masks_, vals_] := Total[vals Transpose@masks];


(* ::Subsection::Closed:: *)
(*NormalizeData*)


SyntaxInformation[NormalizeData] = {"ArgumentsPattern" -> {_,_.}};

NormalizeData[data_] := NormalizeDatai[data, .75]
NormalizeData[data_, quan_?NumberQ] := NormalizeDatai[data, Clip[quan, {0, 1}]]

NormalizeData[data_,minmax_] := ScaleData[data,minmax]


ScaleData = 
  Compile[{{data, _Real, 0}, {range, _Real, 1}}, 
   (data - range[[1]])/(range[[2]] - range[[1]]), 
   RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


NormalizeDatai = Block[{mn},
   Compile[{{data, _Real, 3}, {quant, _Real, 0}},
    mn = Quantile[Cases[Flatten[data], Except[0.]], quant];
    data/mn, {{mn, _Real, 0}}, RuntimeAttributes -> {Listable}, 
    RuntimeOptions -> "Speed"]
   ];


(* ::Subsection::Closed:: *)
(*NormalizeDiffData*)


Options[NormalizeDiffData] = {MeanOutput->False}

SyntaxInformation[NormalizeDiffData] = {"ArgumentsPattern" -> {_, _.}};

NormalizeDiffData[data_,ops:OptionsPattern[]] := NormalizeDiffData[data, Mask[Mean[Transpose[data]], Smoothing -> True],ops]
NormalizeDiffData[data_, mask_,OptionsPattern[]] := Block[{mn,dataout},
	mn = N[MeanNoZero[Flatten[mask data[[All, 1]]]]/100.];
	dataout = data/mn;
	If[OptionValue[MeanOutput],{dataout,mn},dataout]
]


(* ::Subsection::Closed:: *)
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


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
