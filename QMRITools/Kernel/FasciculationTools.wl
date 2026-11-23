(* ::Package:: *)

(* ::Title:: *)
(*QMRITools FasciculationTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`FasciculationTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`FasciculationTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)

NormalizeFasciculationData::usage = 
"NormalizeFasciculationData[data_, mask, {tens, grad, val}] will normalize DWI data based on the tensor fit. 
For each diffusion direction and b-value the estimated unweighted image is calculated."

FindActivations::usage = 
"FindActivations[data] Finds the activation in MUMRI or DTI data after data normalization. 
FindActivations[data, mask] Finds the activation in MUMRI or DTI data after data normalization within the mask."

EvaluateActivation::usage =
"EvaluateActivation[out] allows to evaluate the activation detection using FindActivations, where out is the output of that function with the option ActivationOutput set to True.
EvaluateActivation[out, actS] The same with the extra analysis of the SelectActivations function output given as actS."

AnalyzeActivations::usage = 
"AnalyzeActivations[actMap, mask] Analysis of the activation map generated from the mask.
AnalyzeActivations[actMap, mask, lab] the same but the output is an association with the name lab, or one name per mask for a 4D mask.
AnalyzeActivations[actMap, mask, vox] the same but with the voxel size vox = {z, y, x} in mm, which adds the columns \"ROI Vol [dm3]\", \"Chance [%/dm3]\", \"Chance Obs [%/dm3]\" and the sizes as area in cm^2 \"Mean Area [cm2]\", \"StDv Area [cm2]\", \"Median Area [cm2]\", \"5% Area [cm2]\" and \"95% Area [cm2]\".
AnalyzeActivations[actMap, mask, lab, vox] the same with the name lab and the voxel size vox.

\"Amount\" is the number of activations and \"Amount Obs\" the number of observations (a slice in one volume) that contain at least one activation."

ListActivations::usage =
"ListActivations[act, vox] lists every activation in the activation map act, one association per activation, with the voxel size vox = {z, y, x} in mm. An activation is one connected component within one slice of one volume, counted the same way as in AnalyzeActivations.
ListActivations[act, vox, seg] the same with the label of the label map seg that covers most of the activation, and the percentage of the voxels within seg in that slice that are active in that volume.
ListActivations[act, vox, seg, {norm, mn}] the same with the mean of the normalized data norm and of the mean mn within the activation, as given by FindActivations with ActivationOutput set to All. seg can be None.

The keys \"slice\", \"volume\" and \"component\" identify each activation. The centroid and bounding box are voxel coordinates {y, x} of the array (1-based), a viewer such as ITK-SNAP shows the y axis flipped."

RemoveActivations::usage =
"RemoveActivations[act, keys] removes the activations given by keys from the activation map act, where keys is a list of {slice, volume, component} as given by ListActivations.
act must be the same map that was given to ListActivations, otherwise the component numbers do not match. The output is the activation map without those activations."

SelectActivations::usage = 
"SelectActivations[act] selects the activations above the given ActivationSize.
SelectActivations[act, vox] selects the activations above the given ActivationSize, the voxel size vox = {z, y, x} in mm is needed when SizeMethod is \"Area\".

SelectActivations[act, mask] selects the activations above the given ActivationSize within the given mask or masks. The mask can be 3D or 4D.
SelectActivations[act, {mask, back}] selects the activations above the given ActivationSize within the given mask or masks. All voxels outside the back are ignored.

Output is {actSelected, actTotal} is mask is 3D.
Output is {{actSelected, Total[actSelected]}, {actTotal, Total[actTotal]}} is mask is 4D where actSelected and actTotal are per mask."


(* ::Subsection::Closed:: *)
(*Options*)


ActivationThreshold::usage =
"ActivationThreshold is an option for FindActivations. Fist value is the number of standard deviations second is the percentage threshold."

ThresholdMethod::usage =
"ThresholdMethod is an option for FindActivations. Values can be \"StandardDeviation\", \"Fraction\" or \"Both\"."

IgnoreSlices::usage =
"IgnoreSlices is an option for FindActivations and SelectActivations. determines how many slices of the start and end of the dataset are ignored."

ActivationBackground::usage = 
"ActivationBackground is an option for FindActivations. If all normalized signals, which range between 0-150, are below this value the algorithm does nothing."

ActivationIterations::usage = 
"ActivationIterations is an option for FindActivations. The maximum number of iteration that can be used for activation detection."

ActivationOutput::usage = 
"ActivationOutput is an option for ActivationOutput. If set to All also the mn and threshold values are returned."

ActivationSize::usage = 
"ActivationSize is an option for SelectActivations. Its the minimal size of the activations selected, in number of connected voxels within a slice if SizeMethod is \"Voxels\" and in mm^2 if SizeMethod is \"Area\".
The default is Automatic which is 4 voxels for \"Voxels\" and 30 mm^2 for \"Area\". For \"Area\" the size is rounded up to whole voxels.";

SizeMethod::usage =
"SizeMethod is an option for SelectActivations. Values can be \"Voxels\" or \"Area\". With \"Voxels\" ActivationSize is the number of connected voxels within a slice, with \"Area\" it is the area in mm^2 for which the voxel size is needed.";


(* ::Subsection::Closed:: *)
(*Error Messages*)


FindActivations::tresh = 
"Given thresholds are not valid. The sd should be >1 and is `1` and the fr should be < 1 and is `2`." ;

SelectActivations::dim = 
"The data and/or mask dimensions are wrong. The activation should be 4D, the mask 3/4D, and the background mask 3D. Current dimensions are `1`D, `2`D, and `3`D respectively.";

SelectActivations::size = 
"The data and/or mask size are wrong. The activation should be the same size as the mask and the background mask. Current sizes are `1`, `2`, and `3` respectively.";

SelectActivations::vox =
"SizeMethod \"Area\" needs the voxel size. Give the voxel size {z, y, x} in mm or set SizeMethod to \"Voxels\".";

AnalyzeActivations::size =
"The mask and activation map must have the same Dimensions.";

ListActivations::dim =
"The activation map should be 4D and the label map and data should match it. Current dimensions are `1` and `2`.";

RemoveActivations::key =
"The activation map should be 4D and all keys {slice, volume, component} within its dimensions `1`. Keys out of range: `2`.";


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*NormalizeFasciculationData*)


SyntaxInformation[NormalizeFasciculationData] = {"ArgumentsPattern" -> {_, _, _.}};

NormalizeFasciculationData[data_, mask_] := NormalizeData[data, mask, NormalizeMethod -> "Volumes"];

NormalizeFasciculationData[data_, mask_, {tens_, grad_, val_}] := Block[{
		ran, tensV, coor, scale, bmat
	},
	ran = {0, 1.5 Max[data]};
	(*vectorize the tensor*)
	{tensV, coor} = DataToVector[Transpose@tens, mask];
	bmat = Bmatrix[{val, grad}];
	(*find diffusion weighting scaling*)
	scale = ExpNoZero[-bmat . Transpose[tensV]];
	(*make the output*)
	NormalizeData[Clip[DivideNoZero[data, VectorToData[Transpose[scale], coor]], ran], mask]
]

(* ::Subsection:: *)
(*FindActivations*)


(* ::Subsubsection::Closed:: *)
(*FindActivations*)


Options[FindActivations] = Options[FindActivationsI] = {
	ActivationThreshold -> {3.0, 0.6},
	ThresholdMethod -> "Both", 
	ActivationOutput -> "Activation",
	MaskDilation -> 0, 
	IgnoreSlices -> {0, 0},
	ActivationBackground -> 10,
	ActivationIterations -> 10
};

SyntaxInformation[FindActivations] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

FindActivations[data_?MatrixQ, ops : OptionsPattern[]]:= FindActivationsI[data, ops]

FindActivations[data_, ops : OptionsPattern[]] := FindActivationsI[NormalizeData[data, NormalizeMethod -> "Volumes"],  ops]

FindActivations[data_, mask_, ops : OptionsPattern[]] := FindActivationsI[
		NormalizeData[data, DilateMask[mask, OptionValue[MaskDilation]], NormalizeMethod -> "Volumes"], ops]

FindActivationsI[data_, OptionsPattern[]] := Block[{met, sc, fr, start, stop, dat, act, mn ,tr, itt, back, output},

	(*Get the options*)
	{met, {sc, fr}, {start, stop}, itt, back, output} = OptionValue[{ThresholdMethod, ActivationThreshold, 
		IgnoreSlices, ActivationIterations, ActivationBackground, ActivationOutput}];
		
	(*set the threshold*)
	If[sc < 0 || fr > 1, Return[Message[FindActivations::tresh, sc, fr]]];
	{sc, fr} = Switch[met,
		"Both", {sc, fr},
		"Fraction", {0, fr},
		"StandardDeviation", {sc, 1}
	];

	If[MatrixQ[data],
		FindActC[data, sc, fr, itt, back],

		(*perform the activation finding in the selected slices*)
		dat = RotateDimensionsLeft[Transpose[data[[start + 1 ;; -stop - 1]]]];	
		act = FindActC[dat, sc, fr, itt, back];

		(*create extra output if needed*)
		If[output=!="Activation",
			{mn, tr, sc, fr} = RotateDimensionsRight[MeanThresh[dat, act, sc, fr, back]];
			mn = ToPackedArray@ArrayPad[mn, {{start, stop}, {0, 0}, {0, 0}}, 0.];
			tr = ToPackedArray@ArrayPad[Transpose[{tr, sc, fr}], {{start, stop}, {0, 0}, {0, 0}, {0, 0}}, 0.];
		];

		(*give output*)
		act = SparseArray[ArrayPad[Round[Transpose[RotateDimensionsRight[act]]], {{start, stop}, {0, 0}, {0, 0}, {0, 0}}]];

		If[output==="Activation", {act, data}, {act, data, mn ,tr}] 
	]
]


(* ::Subsubsection::Closed:: *)
(*FindActC*)


FindActC = Compile[{{t, _Real, 1}, {sc, _Real, 0}, {fr, _Real, 0}, {it, _Real, 0}, {back, _Real, 0}}, Block[{
		tz, ti, ts, mn, tr, i, len, cont, err, nz
	},
	mn = Mean[t];
	(*Check background threshold*)
	If[mn <= back, Return[0.0*t]];

	(*zero samples have no data, never use them*)
	nz = Unitize[t];
	tz = Pick[t, nz, 1];
	ti = tz;
	i = 0;
	cont = True;
	err = False;

	While[cont && i <= it,
		i++;
		len = Length[ti];
		mn = Mean[ti];
		(*samples above the fraction threshold for the sd*)
		ts = Pick[tz, UnitStep[fr mn - tz], 0];
		If[Length[ts] <= 5,
			(*Trigger error state*)
			err = True;
			cont = False;
			,
			(*threshold from the sd and the fraction*)
			tr = Max[{0.1, Min[{1.0 - (sc/mn) StandardDeviation[ts], fr}]}];
			(*new non-active samples*)
			ts = Pick[tz, UnitStep[tr mn - tz], 0];
			If[Length[ts] <= 5,
				(*Trigger error state*)
				err = True;
				cont = False;
				,
				(*converged when the number of non-active samples no longer changes*)
				cont = Length[ts] != len;
				ti = ts;
			];
		];
	];

	(*Final threshold from the sd of the non-active samples*)
	If[err,
		0.0 t,
		mn = Mean[ti];
		tr = Min[1.0 - (sc/mn) StandardDeviation[ti], fr];
		nz UnitStep[tr mn - t]
	]
], RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}]


(* ::Subsubsection::Closed:: *)
(*MeanThresh*)


MeanThresh = Compile[{{t, _Real, 1},{s, _Real, 1}, {sc, _Real, 0}, {fr, _Real, 0}, {back, _Real, 0}}, Block[
	{ti, mn, sd, tr},
	mn = Mean[t];
	If[mn <= back,
		{0., 0., 0., 0.},
		(*non activated samples without zeros*)
		ti = Pick[t, s + UnitStep[-t], 0.0];
		If[Length[ti] < 5,
			{0., 0., 0., 0.},
			mn = Mean[ti];
			sd = Ramp[1.0 - (sc/mn)*StandardDeviation[ti]];
			tr = Min[{sd, fr}];
			mn*{1.0, tr, sd, fr}
		]
	]
], RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}]


(* ::Subsection:: *)
(*SelectActivations*)


(* ::Subsubsection::Closed:: *)
(*SelectActivations*)


Options[SelectActivations] = {
	ActivationSize->Automatic,
	SizeMethod->"Voxels",
	IgnoreSlices->{0,0}
};

SyntaxInformation[SelectActivations]={"ArgumentsPattern"->{_,_.,_.,OptionsPattern[]}};

SelectActivations[act_?ArrayQ, ops:OptionsPattern[]]:=SelectActivations[act,{1,1},None, ops]

SelectActivations[act_?ArrayQ, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, ops:OptionsPattern[]]:=SelectActivations[act,{1,1},vox, ops]

SelectActivations[act_ ?ArrayQ, mask_?ArrayQ, ops:OptionsPattern[]]:=SelectActivations[act,{mask,1},None, ops]

SelectActivations[act_ ?ArrayQ, {mask_?ArrayQ, bMask_?ArrayQ}, ops:OptionsPattern[]]:=SelectActivations[act,{mask,bMask},None, ops]

SelectActivations[act_?ArrayQ, {mask_,bMask_}, vox:(None|{_?NumberQ,_?NumberQ,_?NumberQ}), ops:OptionsPattern[]]:=Block[{
	aDepth,mDepth,bDepth,size, met, aDim,nVol, mDim, bDim, start, stop,back, masks,out,sel
	},
	(*check data dimensions*)
	aDepth = ArrayDepth[act];
	mDepth = If[mask=!=1,ArrayDepth[mask],aDepth-1];
	bDepth = If[bMask=!=1,ArrayDepth[bMask],aDepth-1];
	If[aDepth=!=4||4<mDepth<3||bDepth=!=3, Return[Message[SelectActivations::dim,aDepth,mDepth,bDepth]]];

	(*check data sizes*)
	aDim = Dimensions@act[[All,1]];
	nVol = Length@act[[1]];
	mDim = If[mask=!=1, Dimensions@If[mDepth===4, mask[[All,1]], mask], aDim];
	bDim = If[bMask=!=1, Dimensions[bMask], aDim];
	If[aDim=!=mDim || aDim=!=bDim, Return[Message[SelectActivations::size,aDim,mDim,bDim]]];

	(*get detection size of fasc in voxels, the area in mm^2 is rounded up to whole voxels*)
	{met, size} = OptionValue[{SizeMethod, ActivationSize}];
	If[met==="Area" && vox===None, Return[Message[SelectActivations::vox]; $Failed]];
	size = Switch[met,
		"Area", Ceiling[If[size===Automatic, 30, size] / (vox[[2]] vox[[3]])],
		_, Round[If[size===Automatic, 4, size]]
	];

	(*create the selection mask*)
	sel = If[mask===1, 0act[[All,1]]+1, mask ];

	(*create background mask*)
	{start, stop} = OptionValue[IgnoreSlices];
	back = If[bMask===1,If[mDepth===4, Total@Transpose@sel,sel], bMask];
	back[[start+1;;-stop-1]] = back[[start+1;;-stop-1]]+1;

	(*find the masks for which to perform the segmentation*)
	masks = If[mDepth===4, Transpose[MaskData[sel, back]], {sel back}];

	(*make the activation masks per dataset*)
	out = SelectActivationI[act, masks, size];

	If[mDepth===4,
		{{out, Total@out}, {Transpose[Total[Transpose[#]]&/@out], Total@Transpose[Total[out]]}},
		{First@out, Total@Transpose@First@out}
	]
]


(* ::Subsubsection::Closed:: *)
(*SelectActivationI*)


SelectActivationI[act_?ArrayQ, masks_?ArrayQ, size_]:=Block[{msk, dat},
	SparseArray[(
	(*loop over masks for analysis*)
		msk = #;
		(*normal arrays, per slice work on sparse arrays and images is about 3x slower*)
		dat = Round@Transpose[msk # & /@ Transpose[Normal@act]];
		SparseArray[Map[SelectActivationI[#, size]&, dat, {2}]]
	)& /@ masks]
];


(*actual selection function*)
SelectActivationI[im_?MatrixQ, size_?IntegerQ] := SelectActivationI[im , {size, Infinity}]

SelectActivationI[im_?MatrixQ, {min_, max_}] := Block[{lab, cnt},
	If[Total[im, 2] < min, 0 im,
		(*label the components and keep those with a voxel count in range*)
		lab = MorphologicalComponents[im];
		cnt = BinCounts[Flatten@lab, {1, Max[lab] + 1, 1}];
		ArrayReshape[Prepend[Boole[min <= # <= max] & /@ cnt, 0][[Flatten[lab] + 1]], Dimensions[im]]
	]
];


(* ::Subsection:: *)
(*AnalyzeActivations*)


(* ::Subsubsection::Closed:: *)
(*AnalyzeActivations*)


AnalyzeActivations[act_, msk_]:=AnalyzeActivations[act, msk, "", None]

AnalyzeActivations[act_, msk_, vox:{_?NumberQ, _?NumberQ, _?NumberQ}]:=AnalyzeActivations[act, msk, "", vox]

AnalyzeActivations[act_, msk_, lab_]:=AnalyzeActivations[act, msk, lab, None]

AnalyzeActivations[act_, msk_, lab_, vox_]:=Block[{aDepth, aDim, mDepth, mDim, labs},
	aDepth = ArrayDepth[act];
	aDim = Dimensions[If[aDepth===4,act[[All,1]],act[[1,All,1]]]];
	mDepth = ArrayDepth[msk];
	mDim = Dimensions[If[mDepth==3,msk,msk[[All,1]]]];

	If[aDepth =!= (mDepth+1) || aDim =!= mDim, Return[Message[AnalyzeActivations::size]]];

	If[mDepth === 3,
		AnalyzeActivationsI[act, msk, lab, vox],
		labs=If[lab===""||Length[lab]=!=Length[act],
			"Vol_"<>StringPadLeft[ToString[#],3,"0"]&/@Range[Length[act]],
			lab];
		Association[MapThread[AnalyzeActivationsI[#1,#2,#3,vox]&,{act, Transpose@msk, labs}]]
	]
]


(* ::Subsubsection::Closed:: *)
(*AnalyzeActivationsI*)


AnalyzeActivationsI[act_, msk_, lab_, vox_]:=Block[{sizes,sizesO,nActs,nActsO,mSize,mSizeT,nSlices,nVols,nObs,mSd,quants,chance,chanceO,chanceV,vals,out,vol,area},
	
	(*sizes of the activations per observation (slice in one volume)*)
	sizesO = Map[If[Total[Flatten[#]] == 0, {}, ComponentMeasurements[Image[#, "Bit"], "Count"][[All, 2]]]&, act, {2}];
	sizes = N@Flatten[sizesO];

	nActs = Length@sizes;
	(*number of observations with at least one activation*)
	nActsO = Total[Boole[# =!= {}]& /@ Flatten[sizesO, 1]];

	mSize = Map[Total[Flatten[#]]&,msk];
	
	mSizeT = Total@mSize;

	nSlices = Total@Unitize@mSize;
	nVols = Length@act[[1]];
	nObs = nSlices nVols;

	(*mSd=If[nActs>2,{Mean[sizes],StandardDeviation[sizes]},{0,0}];
	quants=If[nActs>2,Quantile[sizes, {0.5,0.05,0.95}],{0,0,0}];*)

	mSd = Which[
		nActs > 2, {Mean[sizes], StandardDeviation[sizes]},
		0 < nActs <=2, {Mean[sizes], 0.},
		True, {0., 0.}];
	quants = If[0 < nActs, Quantile[sizes, {0.5, 0.05, 0.95}], {0. ,0. ,0.}];

	chance = 100. nActs / nVols;
	chanceO = 100. nActs / nObs (* nVols nSlices*);
	chanceV = 1000. chanceO / mSizeT;

	vals=Flatten@{mSizeT, nActs, nObs, chance, chanceO, chanceV, mSd, quants, nActsO};

	out=Association[Thread[{
		"ROI vol", "Amount", "Observed", "Chance/Vol", "Chance/Obs", 
		"Chance/Vox", "Mean Size", "StDv Size", 
		"Median Size", "5% Size", "95% Size", "Amount Obs"}->vals]];

	(*with the voxel size add the volume in dm3, chance per dm3 and the sizes as area in cm2*)
	If[vox =!= None,
		vol = mSizeT Times@@vox/10.^6;
		area = vox[[2]] vox[[3]]/100.;
		out = Join[out, Association[Thread[{
			"ROI Vol [dm3]", "Chance [%/dm3]", "Chance Obs [%/dm3]", "Mean Area [cm2]", "StDv Area [cm2]",
			"Median Area [cm2]", "5% Area [cm2]", "95% Area [cm2]"}->
			Flatten@{vol, If[vol > 0, {chanceO, 100. nActsO/nObs}/vol, {0., 0.}], area mSd, area quants}]]]
	];

	If[lab==="",out, Association[lab->out]]
]


(* ::Subsubsection::Closed:: *)
(*ListActivations*)


SyntaxInformation[ListActivations] = {"ArgumentsPattern" -> {_, _, _., _.}};

ListActivations[act_?ArrayQ, vox:{_?NumberQ, _?NumberQ, _?NumberQ}] := ListActivations[act, vox, None, None]

ListActivations[act_?ArrayQ, vox:{_?NumberQ, _?NumberQ, _?NumberQ}, seg_?ArrayQ] := ListActivations[act, vox, seg, None]

ListActivations[act_?ArrayQ, vox:{_?NumberQ, _?NumberQ, _?NumberQ}, seg:(None|_?ArrayQ), sig:(None|{_?ArrayQ, _?ArrayQ})] := Block[{
		dim, segD, obs, evs, z, vol, pos, lab, nObs, nSl, nMsk, ev, p, row, lv, com, sv, bv
	},
	(*check dimensions*)
	dim = Dimensions[act];
	If[Length[dim] =!= 4 || (seg =!= None && Dimensions[seg] =!= dim[[{1, 3, 4}]]) || (sig =!= None && Dimensions[sig[[1]]] =!= dim),
		Return[Message[ListActivations::dim, dim, If[seg === None, None, Dimensions[seg]]]; $Failed]
	];
	segD = If[seg === None, None, Round@Normal[seg]];

	(*active voxels grouped per observation {slice, volume}*)
	obs = SortBy[GatherBy[SparseArray[act]["NonzeroPositions"], Take[#, 2] &], Take[First[#], 2] &];

	(*components per observation: {slice, volume, component, positions, active voxels in seg}*)
	evs = Flatten[Table[
		{z, vol} = Take[First[o], 2];
		pos = o[[All, 3 ;; 4]];
		lab = Extract[MorphologicalComponents[Normal[act[[z, vol]]]], pos];
		{z, vol, #[[1, 1]], #[[All, 2]], If[segD === None, Length[pos], Total[Unitize[Extract[segD[[z]], pos]]]]} & /@
			SortBy[GatherBy[Transpose[{lab, pos}], First], #[[1, 1]] &]
	, {o, obs}], 1];

	(*events per observation, slices with events per volume and seg voxels per slice*)
	nObs = Counts[evs[[All, {1, 2}]]];
	nSl = Counts[Keys[nObs][[All, 2]]];
	nMsk = If[segD === None, None, Total[Unitize[segD], {2, 3}]];

	(*one association per event*)
	Table[
		ev = evs[[i]];
		{z, vol, p} = ev[[{1, 2, 4}]];
		row = <|
			"event" -> i,
			"slice" -> z,
			"volume" -> vol,
			"component" -> ev[[3]],
			"size [vox]" -> Length[p],
			"area [mm2]" -> Length[p] vox[[2]] vox[[3]],
			"centroid y" -> Round[N@Mean[p[[All, 1]]], 0.1],
			"centroid x" -> Round[N@Mean[p[[All, 2]]], 0.1],
			"bbox ymin" -> Min[p[[All, 1]]],
			"bbox ymax" -> Max[p[[All, 1]]],
			"bbox xmin" -> Min[p[[All, 2]]],
			"bbox xmax" -> Max[p[[All, 2]]]
		|>;

		(*majority label of the event*)
		If[segD =!= None,
			lv = Extract[segD[[z]], p];
			com = First[Commonest[lv]];
			row = Join[row, <|"label" -> com, "label fraction" -> Round[N[Count[lv, com]/Length[p]], 0.01]|>]
		];

		(*signal and baseline in the event*)
		If[sig =!= None,
			sv = N@Mean[Extract[Normal[sig[[1, z, vol]]], p]];
			bv = N@Mean[Extract[Normal[sig[[2, z]]], p]];
			row = Join[row, <|
				"signal" -> Round[sv, 0.1],
				"baseline" -> Round[bv, 0.1],
				"drop [%]" -> If[bv > 0, Round[100 (1 - sv/bv), 0.1], 0.]
			|>]
		];

		(*context of the observation and volume*)
		row = Join[row, <|"events in obs" -> nObs[{z, vol}], "slices in vol" -> nSl[vol]|>];
		If[segD =!= None,
			row = Join[row, <|"obs active [%]" -> If[nMsk[[z]] > 0, Round[100. ev[[5]]/nMsk[[z]], 0.01], 0.]|>]
		];
		row
	, {i, Length[evs]}]
]


(* ::Subsubsection::Closed:: *)
(*RemoveActivations*)


SyntaxInformation[RemoveActivations] = {"ArgumentsPattern" -> {_, _}};

RemoveActivations[act_?ArrayQ, {}] := SparseArray[act]

RemoveActivations[act_?ArrayQ, keys:{{_Integer, _Integer, _Integer}..}] := Block[{
		dim, bad, out, z, vol, im, lab, tab
	},
	(*check the map and the keys*)
	dim = Dimensions[act];
	bad = If[Length[dim] =!= 4, keys, Select[keys, !(1 <= #[[1]] <= dim[[1]] && 1 <= #[[2]] <= dim[[2]] && #[[3]] >= 1) &]];
	If[bad =!= {}, Return[Message[RemoveActivations::key, dim, bad]; $Failed]];

	out = SparseArray[act];
	(*per observation relabel the components as in ListActivations and remove the given ones*)
	Do[
		{z, vol} = Take[First[k], 2];
		im = Normal[out[[z, vol]]];
		lab = MorphologicalComponents[im];
		tab = ConstantArray[1, Max[lab] + 1];
		tab[[Select[k[[All, 3]], # <= Max[lab] &] + 1]] = 0;
		out[[z, vol]] = im ArrayReshape[tab[[Flatten[lab] + 1]], Dimensions[im]];
	, {k, GatherBy[keys, Take[#, 2] &]}];

	out
]


(* ::Subsection::Closed:: *)
(*EvaluateActivation*)


SyntaxInformation[EvaluateActivation]={"ArgumentsPattern"->{_,_.,_.,_.,_.,OptionsPattern[]}};

EvaluateActivation[act_,dat_,mn_,tr_]:=EvaluateActivation[act,dat,mn,tr,act]

EvaluateActivation[{act_,dat_,mn_,tr_}]:=EvaluateActivation[act,dat,mn,tr,act]

EvaluateActivation[{act_,dat_,mn_,tr_},actS_]:=EvaluateActivation[act,dat,mn,tr,actS]

EvaluateActivation[act_,dat_,mn_,tr_,actS_]:=Module[{datD,actD,actSD,mnD,trD,sc,dim,aim
	(*dDim, zz, dd, yy ,xx, sl, dyn*)},
	NotebookClose[plotWindow];

	actD=Normal[act];
	actSD=Normal[actS];
	datD=dat;mnD=mn;trD=tr;
	(*prep images*)
	sc = Quantile[Flatten[datD],{0.005,0.995}];
	dim = N@Clip[Rescale[datD,sc],{0,1}];
	aim = RotateDimensionsLeft[N@{actD-actSD, actSD, 0.actD, Clip[actD+actSD,{0,1}]}];

	dDim = Dimensions[datD];
	{zz, dd, yy, xx} = dDim;
	{sl,dyn} = Ceiling[{zz,dd}/2];

	PrintTemporary["Prepping Manipulate window"];
	pan = Manipulate[
		(*slice location*)
		sl = Clip[sl,{1,zz},{1,zz}];
		dyn = Clip[dyn,{1,dd},{1,dd}];

		l = Length[pos[[m]]];
		n = Clip[n,{1,l}];

		If[!sel,
			(*use selection bar*)
			{z,y,x} = pos[[m,n]];
			sl = z;
			c = Abs[{x,y}-{0,yy+1}];
			,
			(*use locator pane*)
			cor = Reverse[Abs[Ceiling[c]-{0,yy+1}]];
			{y,x} = cor;
			z = sl;
		];

		(*get the signals*)
		dDat = Range[dd];
		sig = datD[[z,All,y,x]];
		actT = actD[[z,All,y,x]];
		actSt = actSD[[z,All,y,x]];

		(*get the line values*)
		mean = mnD[[z,y,x]];
		tresh = trD[[z,All,y,x]];

		Row[{Show[
			ListLinePlot[Thread[{dDat,sig}],PlotStyle->LightDarkV[],PlotMarkers->Automatic,PlotRange->{0,1.1 Max@sig},ImageSize->500,
				GridLines->If[grid, {
					{{dyn, Directive[LightDarkV[], Thick]}},
					{{mean, Directive[LightDarkV[], Thick]},
						{tresh[[1]], Directive[Red,Thick]},
						{tresh[[2]], Directive[Thick,RGBColor[0.812807, 0.518694, 0.303459], Dashed]},
						{tresh[[3]], Directive[Thick,RGBColor[0.253651, 0.344893, 0.558151], Dashed]}}
					}, None
				], $plotOptions
			],
			Graphics[{
				RGBColor[0.812807, 0.518694, 0.303459], Line[{{0, tresh[[1]]}, {Length@sig, tresh[[1]]}}],
				RGBColor[0.253651, 0.344893, 0.558151], Line[{{0, tresh[[1]]}, {Length@sig, tresh[[1]]}}]
			}],
			ListPlot[Pick[Thread[{dDat,sig}],actT,1],PlotStyle->Red,PlotMarkers->{Automatic,10}],
			ListPlot[Pick[Thread[{dDat,sig}],actSt,1],PlotStyle->Green,PlotMarkers->{Automatic,10}]]
			,
			LocatorPane[Dynamic[c],ImageCompose[
				Image[dim[[sl,dyn]],ColorSpace->"Grayscale",ImageSize->400],{Image[aim[[sl,dyn]],ColorSpace->"RGB"], alpha}
			], Appearance->Style[If[crs, "+", " "], Blue, FontSize->60, Bold]
		]}]

		,{{sel, False, "Use Locator"},{True,False}}
		,Delimiter
		,{{m, 1, "Number of act."},max,ControlType->SetterBar}
		,{{n, 1, "Position number"}, 1, Dynamic[l], 1}
		,Delimiter
		,{{sl, 1, "Slice"}, 1, Dynamic[zz], 1}
		,{{dyn, 1, "Dynamic"}, 1, Dynamic[dd], 1}
		,Delimiter
		,{{alpha, 0.5, "Opacity"},0,1}
		,{{crs, True, "Show cross"},{True,False}}
		,Delimiter
		,{{grid, False, "Show grid"},{True,False}}
		,Row[{Dynamic[{z,y,x}]}]
		,
		{zz,ControlType->None},{dd,ControlType->None},{xx,ControlType->None},{yy,ControlType->None},
		{max,ControlType->None},{pos,ControlType->None},{l,1,ControlType->None},

		{z,ControlType->None},{y,ControlType->None},{x,ControlType->None},
		{c,ControlType->None},{cor,ControlType->None},
		{dDat,ControlType->None},{sig,ControlType->None},
		{actT,ControlType->None},{actSt,ControlType->None},
		{mean,ControlType->None},{tresh,ControlType->None},

		Initialization:>{
			{zz, dd, yy, xx} = dDim,
			{sl,dyn} = Ceiling[{zz,dd}/2],
			max = Select[Sort[DeleteDuplicates[Flatten@Total@Transpose@actSD]],#>0&],
			pos = Position[Round@Total@Transpose@actSD,#]&/@max
		},

		SynchronousInitialization->True,
		SaveDefinitions->False,
		ControlPlacement->Right
	];

	plotWindow = CreateWindow[DialogNotebook[{CancelButton["Close",Clear[datD,actD,actSD,mnD,trD,dim,aim];DialogReturn[]],pan},WindowSize->All,WindowTitle->"Plot data window"]];
]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
