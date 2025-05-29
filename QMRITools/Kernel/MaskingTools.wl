(* ::Package:: *)

(* ::Title:: *)
(*QMRITools MaskingTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`MaskingTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`MaskingTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


NormalizeData::usage = 
"NormalizeData[data] normalizes the data to the mean signal of the data. For 4D data it normalizes to the first volume of the 4th dimension.
NormalizeData[data, mask] normalizes the data based on the mean signal only within the mask."

NormalizeMeanData::usage = 
"NormalizeMeanData[data] calculates the mean normalized data from a 4D dataset."

HomogenizeData::usage = 
"HomogenizeData[data, mask] tries to homogenize the data within the mask by removing intensity gradients."

FitGradientMap::usage = 
"FitGradientMap[data, ord] fit of gradient trough all non zero values withing the data."


Mask::usage =
"Mask[data] creates a mask by automatically finding a threshold.
Mask[data, min]creates a mask which selects only data above the min value.
Mask[data,{min,max}] creates a mask which selects data between the min and max value."

SmoothMask::usage = 
"SmoothMask[mask] generates one clean masked volume form a noisy mask."

MaskData::usage = 
"MaskData[data, mask] applies a mask to data. mask can be 2D or 3D, data can be 2D, 3D or 4D."

MaskSegmentation::usage = 
"MaskSegmentation[seg, mask] applies a mask to a splited segmentation seg from SplitSegmentations. 
The mask is 3D, seg is 4D."


SplitSegmentations::usage = 
"SplitSegmentations[segmentation] splits a lable mask from ITKsnap or slicer3D in separate masks and label numbers.
Output is masks and label numbers, {mask, labs}."

GetSegmentationLabels::usage = 
"GetSegmentationLabels[segmentation] gives a list of all labels in the segmentation."

RescaleSegmentation::usage = 
"RescaleSegmentation[data, dim] rescales segmentations to given dimensions.
RescaleSegmentation[data, {vox1, vox2}] rescales segmentations from voxelsize vox1 to voxelsize vox2."

MergeSegmentations::usage = 
"MergeSegmentations[masks, labels] generates an ITKsnap or slices3D compatible segmentation from individual masks and label numbers.
Output is a labled segmentation.
MergeSegmentations[masks] does the same but automatically numbers the segmentations."

JoinSegmentations::usage =
"JoinSegmentations[seg, joinRules] joins the segmentations in seg according to the rules in joinRules.
JoinRules is a list of rules {{join, new}..} where join is a list of labels to be joined and new is the new label number.
For example {{1, 2}, 3} joins the labels 1 and 2 to label 3."

SelectSegmentations::usage =
"SelectSegmentations[seg, labs] selects only the segmentations from seg with label number labs."

ReplaceSegmentations::usage =
"ReplaceSegmentations[seg, labs, new] replaces the labels labs form the segmentation seg for labels new. Both labs and new should
be lists of integers of the same size. If seg contains more labels then given in labs these will be replaced by 0." 

SmoothSegmentation::usage =
"SmoothSegmentation[segmentation] smooths segmentations and removes the overlaps between multiple segmentations.
SmoothSegmentation[segmentation, labs] only smooths the selected label number labs." 

RemoveMaskOverlaps::usage = 
"RemoveMaskOverlaps[mask] removes the overlaps between multiple masks. Mask is a 4D dataset with {z, masks, x, y}."

GetCommonSegmentation::usage = 
"GetCommonSegmentation[dat, seg, vox] For a list of multiple datasets dat the common segmentations from the list seg are determined.
Output is a list of segmentations where for each region only the part present in all datasets is selected."


DilateMask::usage=
"DilateMask[mask,size] if size > 0 the mask is dilated and if size < 0 the mask is eroded."

SelectMaskComponents::usage=
"SelectMaskComponents[mask] selects the largest connected component in the mask.
SelectMaskComponents[mask,n] selects the n largest connected components in the mask."

SegmentMask::usage = 
"SegmentMask[mask, n] divides a mask in n segments along the slice direction, n must be an integer. The mask is divided in n equal parts where each parts has the same number of slices."

SegmentationVolume::usage = 
"SegmentationVolume[seg] calculates the volume of each label in the segmentation."


(* ::Subsection::Closed:: *)
(*Options*)


NormalizeMethod::usage = 
"NormalizeMethod is an option for NormalizeData. Can be \"Set\" or \"Volumes\" wich normalizes to the first volume or normalizes each volume individually, respectively.
If \"Uniform\" normalizes the histogram of the data to have a uniform distribution between 0 and 1 where 0 is treated as background of the data."

FitOrder::usage = 
"FitOrder is an option for HomogenizeData. It specifies the order of harmonics to be used for the homogenization."

MaskSmoothing::usage = 
"MaskSmoothing is an options for Mask, SmoothMask and SmoothSegmentation, if set to True it smooths the mask, by closing holse and smoothing the contours."

MaskComponents::usage =
"MaskComponents is an option for Mask, SmoothMask and SmoothSegmentation. Determines the amount of largest clusters used as mask." 

MaskClosing::usage =
"MaskClosing  is an option for Mask, SmoothMask and SmoothSegmentation. The size of the holes in the mask that will be closed." 

MaskDilation::usage = 
"MaskDilation is an option for Mask, SmoothMask and SmoothSegmentation. If the value is greater than 0 it will dilate the mask, if the value is smaller than 0 it will erode the mask."

MaskFiltKernel::usage =
"MaskFiltKernel is an option for Mask, SmoothMask and SmoothSegmentation. How mucht the contours are smoothed." 

SmoothIterations::usage =
"SmoothIterations is an option for Mask, SmoothMask and SmoothSegmentation and defines how often the smoothing is repeated."


(* ::Subsection:: *)
(*Error Messages*)


Mask::tresh = "Given threshold `1` value is not a valid input, must be a number for min threshold only or a vector {min tresh, max tresh}."

MaskData::dim = "Dimensions are not equal, data: `1`, mask `2`." 

MaskData::dep = "Data dimensions should be 2D, 3D or 4D. Mask dimensions should be 2D or 3D. Data is `1`D and Mask is `2`D."


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection:: *)
(*Normalize Functions*)


(* ::Subsubsection::Closed:: *)
(*NormalizeData*)


SyntaxInformation[NormalizeData] = {"ArgumentsPattern" -> {_,_., OptionsPattern[]}};

Options[NormalizeData] = {NormalizeMethod -> "Set"}

SyntaxInformation[NormalizeData] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

NormalizeData[data_, opts : OptionsPattern[]] := Block[{ndat, mask},
	If[OptionValue[NormalizeMethod]=!="Uniform",
		ndat = Switch[ArrayDepth[data], 3, data, 4, Mean@Transpose@data];
		mask = Mask[NormDat[ndat - Min@ndat, 0.]];
		NormalizeData[data, mask, opts]
		,
		Quiet@NormDatC[data]
	]
]

NormalizeData[data_, msk_, opts : OptionsPattern[]] := Block[{dat, mn, min, dato, mask},
	mask = Normal@msk;
	dato = data - Min[data];
	mn = Switch[ArrayDepth[data],
		3, MeanNoZero[Flatten[mask dato]],
		4, Switch[OptionValue[NormalizeMethod],
			"Volumes", MedianNoZero[Flatten[mask #]] & /@ Transpose[dato],
			_, MedianNoZero[Flatten[mask dato[[All, 1]]]]
		]
	];
	NormDat[dato, mn]
]


(*normalization towards uniform distribution*)
NormDatC = Compile[{{dat, _Real, 3}}, Block[{fl, min, max, bins, cdf, n, tot},
	n = 1024;
	fl = Flatten[dat];
	{min, max} = MinMax[fl];
	If[min === max, dat,
		bins = BinCounts[fl, {min, max, (max - min)/n}];
		bins[[1]] = 0.;
		tot = Total[bins];
		If[tot == 0., dat,
			cdf = Accumulate[bins]/tot;
			Map[cdf[[#]] &, Clip[Floor[(n - 1) (dat - min)/(max - min) + 1], {1, n}, {1, n}], {-2}]]
		]
	],
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", Parallelization -> True]


(*normalization towards median of given value*)
NormDat[dat_, mn_] := ToPackedArray[100. Which[
	mn===0., dat/MedianNoZero[Flatten@dat],
	ListQ[mn],Transpose[Transpose[dat]/mn],
	True, dat/mn]
]


(* ::Subsubsection::Closed:: *)
(*NormalizeMeanData*)


SyntaxInformation[NormalizeMeanData] = {"ArgumentsPattern" -> {_,_., OptionsPattern[]}};

Options[NormalizeMeanData] = Options[NormalizeData]

NormalizeMeanData[data_, opts:OptionsPattern[]] := NormalizeData[Mean@Transpose@data, opts]

NormalizeMeanData[data_, mask_, opts:OptionsPattern[]] := NormalizeData[Mean@Transpose@data, mask, opts]


(* ::Subsubsection::Closed:: *)
(*HomogenizeData*)


Options[HomogenizeData] = {FitOrder->5}

SyntaxInformation[HomogenizeData] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

HomogenizeData[dat_, opts:OptionsPattern[]]:=HomogenizeData[dat, Unitize[dat], opts]

HomogenizeData[dat_, mask_, OptionsPattern[]] := Block[{fit},
	fit = FitGradientMap[Erosion[mask, 1] GaussianFilter[dat, 5], OptionValue[FitOrder]];
	NormalizeData[Clip[mask DivideNoZero[dat, fit], {0, 3}, {0, 0}]]
]


(* ::Subsection:: *)
(*Masking*)


(* ::Subsubsection::Closed:: *)
(*Mask*)


Options[Mask] = {MaskSmoothing -> False, MaskComponents -> 2, MaskClosing -> False, MaskFiltKernel -> 2, MaskDilation -> 0, SmoothIterations->3};

SyntaxInformation[Mask] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

Mask[dat_?ArrayQ, opts : OptionsPattern[]] := Mask[dat, {0, 0}, opts]

Mask[dat_?ArrayQ, tr_?NumberQ, opts : OptionsPattern[]] := Mask[dat, {tr, 1.1 Max[dat]}, opts]

Mask[inp : {{_?ArrayQ, _} ..}, opts : OptionsPattern[]] := Times @@ (Mask[#, opts] & /@ inp)

Mask[{dat_?ArrayQ, tr_?NumberQ}, opts : OptionsPattern[]] := Mask[dat, tr, opts]

Mask[{dat_?ArrayQ, tr_?VectorQ}, opts : OptionsPattern[]] := Mask[dat, tr, opts]

Mask[dat_?ArrayQ, tr_?VectorQ, opts:OptionsPattern[]]:= Block[{mask, tresh, dataD, datN, data, dil},

	(*perform data checks*)
	data = ToPackedArray@N@Normal@dat;
	dataD = ArrayDepth[data];
	If[Length[tr] =!= 2, Return@Message[Mask::tresh, tr]];
	If[ArrayDepth[data] > 3, Return@Message[Mask::dep, dataD]];

	(*perform the masking*)		
	mask = If[tr === {0, 0},
		(*no Threshold*)
		ImageData[Binarize[If[dataD == 2, Image, Image3D][Rescale[data, {1, 0.95} MinMax[data]]]]],
		(*Threshold*)
		UnitStep[data - tr[[1]]] - UnitStep[data - tr[[2]]]
	];

	(*smooth the mask if needed*)
	Round@Normal@If[OptionValue[MaskSmoothing], SmoothMask[mask, FilterRules[{opts, Options[Mask]}, Options[SmoothMask]]], mask]
]


(* ::Subsubsection::Closed:: *)
(*SmoothMask*)


Options[SmoothMask] = {
	MaskComponents->1, 
	MaskClosing->True, 
	MaskFiltKernel->2, 
	MaskDilation -> 0, 
	SmoothIterations -> 3
}

SyntaxInformation[SmoothMask] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

SmoothMask[mask_, OptionsPattern[]] := Block[{dil ,obj, close, itt, ker, maskI, dim, crp},
	(*get the options*)
	dil = OptionValue[MaskDilation];(*how much the mask is eroded or dilated*)
	dil = If[NumberQ[dil],  Round@dil];
	obj = OptionValue[MaskComponents];(*number of objects that are maintained*)
	obj = If[NumberQ[obj], Round[obj]];
	close = OptionValue[MaskClosing];(*close holes in mask*)
	If[IntegerQ[close], If[close>=1, close = True, close = False]];(*for legacy purposes*)
	itt = OptionValue[SmoothIterations];(*how much the smoothing is repeated*)
	ker = OptionValue[MaskFiltKernel];(*how much smoothing*)

	(*prep the segmentation*)
	dim = Dimensions@mask;
	{maskI, crp} = AutoCropData[mask];
	maskI = Image3D[NumericArray[maskI, "Integer8"]];

	(*remove the mask holes*)
	If[close===True, maskI = 1 - TakeObject[1 - maskI, 1]];

	(*select number of mask components*)
	If[obj>=1, maskI = TakeObject[maskI, obj]];

	(*Dilate or Erode the mask*)
	Which[
		dil > 0, maskI = Dilation[maskI, dil],
		dil < 0, maskI = Erosion[maskI, -dil]
	];

	(*Filter the segmentation*)
	If[ker > 0, Which[
		itt === 1, maskI = Erosion[Round[GaussianFilter[Dilation[maskI, 1], ker]], 1],
		itt > 1, maskI = Nest[Erosion[Round[GaussianFilter[Dilation[#, 1], ker] + 0.05], 1] &, maskI, itt]
	]];

	(*reverse cropping and return the mask*)
	ReverseCrop[SparseArray@ImageData@Round@maskI, dim, crp]
]


(* ::Subsubsection::Closed:: *)
(*TakeObject*)


TakeObject[maskI_, obj_] := Block[{morph, keys},
	morph = MorphologicalComponents[maskI, CornerNeighbors -> False];
	If[Max[morph] <= obj, maskI,
		keys = Reverse[k=Keys[SortBy[ComponentMeasurements[morph, "Count", CornerNeighbors -> False], Last]]][[;; obj]];
		Image3D[NumericArray[Total[SparseArray[1 - Unitize[morph - #]] & /@ keys], "Integer8"]]
	]
];

(* ::Subsubsection::Closed:: *)
(*DilateMask*)


SyntaxInformation[DilateMask] = {"ArgumentsPattern" -> {_, _}};

DilateMask[mask_, size_] := If[ArrayDepth[mask]===3,
	SmoothMask[mask, MaskDilation -> size, MaskComponents -> Infinity, MaskClosing -> False, SmoothIterations -> 0],
	Transpose[DilateMask[#,size]&/@Transpose[mask]]
]


(* ::Subsubsection::Closed:: *)
(*DilateMask*)


SyntaxInformation[SelectMaskComponents] = {"ArgumentsPattern" -> {_, _.}};

SelectMaskComponents[mask_] := SelectMaskComponents[mask, 1]

SelectMaskComponents[mask_, n_] := SmoothMask[mask, MaskComponents -> n,
	MaskDilation -> 0, MaskClosing -> False, SmoothIterations -> 0]


(* ::Subsubsection::Closed:: *)
(*MaskData*)


SyntaxInformation[MaskData] = {"ArgumentsPattern" -> {_, _}};

MaskData[data_, mask_]:=Block[{dataD, maskD,dimD,dimM,out},
	(*get the data properties*)
	dataD = ArrayDepth[data];
	maskD = ArrayDepth[mask];
	dimD = Dimensions[data];
	dimM = Dimensions[mask];

	(*determine how to mask the data*)
	out = Switch[{dataD, maskD},
		{2,2}, If[dimD == dimM, mask data, 1],
		{3,3}, If[dimD == dimM, mask data, 1],
		{3,2}, If[dimD[[2;;]] == dimM, mask # &/@ data, 1],
		{4,3}, Which[
			dimD[[2;;]] == dimM, mask # & /@ data, 
			dimD[[{1,3,4}]] == dimM, Transpose[mask # & /@ Transpose[data]], 
			True, Return[Message[MaskData::dim, dimD, dimM];$Failed]],
		_, Return[Message[MaskData::dep, dataD, maskD];,$Failed]
	];

	(*make the output*)
	ToPackedArray@N@Normal@out
]


(* ::Subsubsection::Closed:: *)
(*MaskSegmentation*)


SyntaxInformation[MaskSegmentation] = {"ArgumentsPattern" -> {_, _}};

MaskSegmentation[seg_, mask_] := Block[{msk, dim, sc, cr},
	msk = Normal@mask;
	dim = Dimensions[seg[[All, 1]]];
	Transpose[(
		{sc, cr} = AutoCropData[#];
		SparseArray@Round@ReverseCrop[ApplyCrop[msk, cr] sc, dim, cr]
	) & /@ Transpose[seg]]
]


(* ::Subsection:: *)
(*Segmentation functions*)


(* ::Subsubsection::Closed:: *)
(*GetSegmentationLabels*)


GetSegmentationLabels[segI_]:= Sort[DeleteDuplicates[SparseArray[Round[segI]]["ExplicitValues"]]];


(* ::Subsubsection::Closed:: *)
(*SplitSegmentations*)


SyntaxInformation[SplitSegmentations] = {"ArgumentsPattern" -> {_}};

SplitSegmentations[segI_]:=SplitSegmentations[segI, True]

SplitSegmentations[segI_, sparse_] := Block[{seg, dim, exVals, exPos, vals},
	seg = SparseArray[Round[segI]];
	dim = Dimensions[seg];

	exVals = seg["ExplicitValues"];
	exPos = seg["ExplicitPositions"];
	vals = Sort@DeleteDuplicates@exVals;

	seg = Transpose[SparseArray[Pick[exPos, Unitize[exVals - #], 0] -> 1, dim] & /@ vals];

	{If[sparse, seg, Normal@seg], vals}
]


(* ::Subsubsection::Closed:: *)
(*MergeSegmentations*)


SyntaxInformation[MergeSegmentations] = {"ArgumentsPattern" -> {_, _.}};

MergeSegmentations[seg_]:= MergeSegmentations[seg, Range[Length[First@seg]]]

MergeSegmentations[seg_, lab_] := Block[{mt, nv},
	mt = Transpose[If[!SparseArrayQ@seg, SparseArray@Round@seg, seg]];
	(*mapping over sparse is quicker than direct multiplication and total*)
	Normal[Total[mt[[#]] lab[[#]] & /@ Range[Length@lab]] (1 - UnitStep[Total[# & /@ mt] - 2])]
]


(* ::Subsubsection::Closed:: *)
(*JoinSegmentations*)


JoinSegmentations[segI_, joinRules : {_?ListQ, _?IntegerQ}] := JoinSegmentations[segI, {joinRules}]

JoinSegmentations[segI_, joinRules : {{_?ListQ, _?IntegerQ} ..}] := Block[{seg, lab, keep, keepL, join, new, newL},
	{seg, lab} = SplitSegmentations[segI];

	{keep, keepL} = SelectSegmentations[{seg, lab}, Complement[lab, Flatten[joinRules[[All, 1]]]]];

	{new, newL} = Transpose[(
		{join, newL} = #;
		If[MemberQ[lab, newL] && ! MemberQ[join, newL],	
			Echo[{join, new}, "Skipping, new label is not uniuqe or part of replaced: "]; 0,
			new = Total@Transpose@First@SelectSegmentations[{seg, lab}, #[[1]]];
			{new, newL}
		]
	) & /@ joinRules];

	If[keepL==={},
		MergeSegmentations[Transpose@new,newL],
		MergeSegmentations[Transpose[Join[Transpose@keep, new]], Join[keepL, newL]]
	]
]


(* ::Subsubsection::Closed:: *)
(*SelectSegmentations*)


SyntaxInformation[SelectSegmentations] = {"ArgumentsPattern" -> {_,_}};

SelectSegmentations[seg_, labSel_] := SelectReplaceSegmentations[seg, labSel, labSel]


(* ::Subsubsection::Closed:: *)
(*ReplaceSegmentations*)


SyntaxInformation[ReplaceSegmentations] = {"ArgumentsPattern" -> {_,_,_}};

ReplaceSegmentations[segm_, labSel_, labNew_] := SelectReplaceSegmentations[segm, labSel, labNew]


(* ::Subsubsection::Closed:: *)
(*SelectReplaceSegmentations*)


SelectReplaceSegmentations[segm_, labSel_, labNew_] := Block[{split, seg, lab, sel},
	split = If[Length[segm] == 2, If[VectorQ[segm[[2]]], False, True], True];
	{seg, lab} = If[split,  SplitSegmentations[segm], segm];

	sel = MemberQ[labSel, #] & /@ lab;

	If[AllTrue[sel, # === False &],
		If[split, 0 MergeSegmentations[seg, lab], {0 seg[[All,1]], {}}]
		,
		seg = Transpose[Pick[Transpose[seg], sel, True]];
		lab = Pick[lab /. Thread[labSel->labNew], sel, True];
		or = Ordering[lab];

		If[split, MergeSegmentations[seg, lab], {seg[[All,or]], lab[[or]]}]
	]
]


(* ::Subsubsection::Closed:: *)
(*SmoothSegmentation*)


Options[SmoothSegmentation] = {MaskComponents -> 1, MaskClosing -> 1, MaskFiltKernel -> 2, MaskDilation -> 0, SmoothIterations->1}

SyntaxInformation[SmoothSegmentation] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

SmoothSegmentation[maskIn_,  opts:OptionsPattern[]] :=SmoothSegmentation[maskIn, All, opts]

SmoothSegmentation[maskIn_, what_, opts:OptionsPattern[]] := Block[{smooth, obj, md, masks, labs},
	md = ArrayDepth[maskIn];
	(*split segmentations and make sparse if needed*)
	If[md===3, {masks, labs} = SplitSegmentations[maskIn], masks = maskIn];
	masks = Transpose[If[!SparseArrayQ@masks, SparseArray@Round@masks, masks]];

	(*Get smoothed segmentations of the given labels and merge*)
	masks[[what]] = SmoothMask[#, Sequence@@FilterRules[{opts}, Options[SmoothMask]]]&/@masks[[what]];
	If[md===3, MergeSegmentations[Transpose[masks], labs], masks]
]


(* ::Subsubsection::Closed:: *)
(*RescaleSegmentation*)


SyntaxInformation[RescaleSegmentation] = {"ArgumentsPattern" -> {_, _}};

RescaleSegmentation[seg_, vox_] := Block[{segs, val},
	If[ArrayDepth[seg] == 3, {segs, val} = SplitSegmentations[seg], segs = seg];
	segs = Transpose@RemoveMaskOverlapsI[SparseArray[SparseArray[Round[RescaleData[Normal[#], vox, InterpolationOrder -> 1]]] & /@Transpose[segs]]];
	If[ArrayDepth[seg] == 3, MergeSegmentations[segs, val], segs]
]


(* ::Subsubsection::Closed:: *)
(*RemoveMaskOverlaps*)


SyntaxInformation[RemoveMaskOverlaps] = {"ArgumentsPattern" -> {_}};
	
RemoveMaskOverlaps[masks_] := Transpose@RemoveMaskOverlapsI[Transpose[SparseArray[Round@masks]]];

RemoveMaskOverlapsI[masks_] := Block[{maskOver, posOver, maskInp, maskOut, z, x, y, p},
	maskOver = SparseArray[Mask[Total[masks], 1.5]];
	posOver = maskOver["ExplicitPositions"];
	maskOut = SparseArray[SparseArray[(1 - maskOver) #] & /@ masks];

	maskInp = maskOver Transpose[masks, {4, 1, 2, 3}];
	p = ({z, x, y} = #; {First@First[maskInp[[z, x, y]]["ExplicitPositions"]], z, x, y}) & /@posOver;
	maskOut + SparseArray[p -> 1, Dimensions[maskOut]]
]


(* ::Subsection::Closed:: *)
(*GetCommonSegmentation*)


SyntaxInformation[GetCommonSegmentation] = {"ArgumentsPattern" -> {_, _, _}};

GetCommonSegmentation[dat:{_?ArrayQ ..}, seg:{_?ArrayQ ..}, vox:{_?ListQ ..}] := Block[
	{dims, datC, cr, len, segs, labs, labAll, labSel, l, gr, segR},

	(*auto crop all the datasets*)
	dims = Dimensions /@ dat;
	{datC, cr} = Transpose[AutoCropData /@ dat];
	len = Length[datC];

	(*split and smooth all the segmentations*)
	{segs, labs} = Transpose[SplitSegmentations[ApplyCrop[#[[1]], #[[2]]]] & /@ Thread[{seg, cr}]];
	segs = SmoothSegmentation[#, MaskComponents -> 1, SmoothIterations -> 2] & /@ segs;

	(*find common labels*)
	labAll = Intersection @@ labs;
	labSel = (l = #; Flatten[Position[l, #] & /@ labAll]) & /@ labs;

	gr = 5 Mean@vox;
	(*find the common segmentation for each dataset*)
	Table[
		segR = Times @@ Table[
			If[tar === mov,
				segs[[tar, All, labSel[[tar]]]],
				Round@Last@RegisterDataTransform[
					{datC[[tar]], vox[[tar]]},
					{datC[[mov]], vox[[mov]]},
					{segs[[mov, All, labSel[[mov]]]], vox[[mov]]},
					MethodReg -> {"rigid", "affine", "bspline"}, Resolutions -> 3, NumberSamples -> 5000, BsplineSpacing -> gr
				]
			]
		, {mov, 1, len}];

		MergeSegmentations[ReverseCrop[RemoveMaskOverlaps[segR], dims[[tar]], cr[[tar]]], labAll[[labSel[[tar]]]]]
	, {tar, 1, len}]
]


(* ::Subsection::Closed:: *)
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


(* ::Subsection::Closed:: *)
(*SegmentationVolume*)


SyntaxInformation[SegmentationVolume] = {"ArgumentsPattern" -> {_, _.}};

SegmentationVolume[seg_] := SegmentationVolume[seg, {0, 0, 0}]

SegmentationVolume[seg_, vox : {_?NumberQ, _?NumberQ, _?NumberQ}] := Block[{vol},
	vol = If[vox === {0, 0, 0}, 1, N@((Times @@ vox)/1000)];
	vol  Total[Flatten[#]] & /@ Switch[ArrayDepth[seg],
		3, Transpose[First@SplitSegmentations[seg]],
		4, Transpose[seg],
		_, Return[$Failed]
	]
]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
