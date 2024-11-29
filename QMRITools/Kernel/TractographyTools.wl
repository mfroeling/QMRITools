(* ::Package:: *)

(* ::Title:: *)
(*QMRITools TractographyTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`TractographyTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`TractographyTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


FiberTractography::usage =
"FiberTractography[tensor, vox] performs fiber tractography on the tensor with voxel dimensions vox.
FiberTractography[tensor, vox, {par, {min, max}}] performs fiber tractography on the tensor with voxel dimensions vox with additional stopping criteria
par, where tracts are only generated between values of par min and max.
FiberTractography[tensor, vox, {{par, {min, max}}, ..}] performs fiber tractography on the tensor with voxel dimensions vox with 
multiple additional stopping criteria."

FitTracts::usage = 
"FitTracts[tract] fits a tract or a list of tracts, each defined as a list of {x, y, z} coordinates with a polynomial function.
FitTracts[tract, vox, dim] does the same but constrains all the tract coordinates to the volume defined by dim."

ResampleTracts::usage = 
"ResampleTracts[tracts, n] resample each Tract to exactly n vertices."

MoveTracts::usage = 
"MoveTracts[tracts, off] moves the tract coordinates by off, which is {x, y, z}."

RescaleTracts::usage = 
"RescaleTracts[tracts, sc] scales the tract coordinates by 1/sc, which is {x, y, z} or single number."

GetTractValues::usage = 
"GetTractValues[tracts, parameter, vox] gets the value of the parameter map at each tract coordinate."


FilterTracts::usage = 
"FilterTracts[tracts, vox, {select.. }] filters the tracts based on the list of select criteria.
Select criteria are defined as {\"logic\",{\"how\", criteria}}.
The \"logic\" parameter can be \"and\", \"or\" and \"not\".
The \"how\" parameter can be:
	- \"x\", \"y\", or \"z\" for slice selection, here criteria is a slice number
	- \"through\" for selecting tract that go through a roi, here criteria is a 3D mask.
	- \"within\" for selecting tract that fit fully within the roi, here criteria is a 3D mask.
	- \"partwithin\" for selecting the part of the tracts that fall within the roi, here criteria is a 3D mask.
Any number of select criteria can be listed."

SegmentTracts::usage = 
"SegmentTracts[tracts, segs, vox, dim] segments the tracts based on segs."


SeedDensityMap::usage = 
"SeedDensityMap[seeds, vox, dim] makes a seed density map based on the seed locations."

TractDensityMap::usage = 
"TractDensityMap[tracts, vox, dim] makes a tract density map based on the tracts vertices."

TractLengthMap::usage = 
"TractLengthMap[tracts, vox, dim] makes a tract length map based on the tracts lengths."

TractLength::usage = 
"TractLength[tracts] calculates the length of each tract."

TractAngleMap::usage = 
"TractAngleMap[tracts, vox, dim] makes a tract angle map based on the tracts angles with the z-plane.
TractAngleMap[tracts, v1, vox, dim] makes a tract angle map based on the tracts angles with the plane normal to v1.
TractAngleMap[tracts, {v1, v2}, vox, dim] makes a tract angle map based on the tracts elevation angles with the plane normal to v1 and the azimuth angle in that plane relative to v2."

TractAngle ::usage = 
"TractAngle[tracts] calculates the angle of each tract segment with the z-plane.
TractAngle[tracts, v1] calculates the angle of each tract segment with the plane normal to v1.
TractAngle[tracts, {v1, v2}] calculates the elevation and azimuth angle of each tract segment with the plane normal to v1 and the azimuth angle in that plane relative to v2."

TractCurvatureMap::usage =
"TractCurvatureMap[tracts, vox, dim] makes a tract curvature map based on the tracts curvature."

TractCurvature::usage =
"TractCurvature[tracts] calculates the curvature of each tract segment."


PlotTracts::usage = 
"PlotTracts[tracts, vox] plots the tracts assuming a BoxRatio based on vox.
PlotTracts[tracts, vox, dim] plots the tracts assuming a BoxRatio based on vox with a PlotRange spanning the full dim."

PlotSegmentedTracts::usage = 
"PlotSegmentedTracts[tracts, segments, dim, vox] plots the tracts after segmenting each segment.
PlotSegmentedTracts[tracts, segments, bones, dim, vox] plots the tracts after segmenting each segment also rendering a bone volume."


FindTensorPermutation::usage = 
"FindTensorPermutation[tensor, vox] performs tractography for all tensor permutations and gives back the one that has the longest tracts.
FindTensorPermutation[tensor, vox, {par, {min, max}}] same but with additional stopping criteria par, where tracts are only generated between values of par min and max.
FindTensorPermutation[tensor, vox, {{par, {min, max}}, ..}] same but with multiple additional stopping criteria.

Output = {permutations, flips, plot}

FindTensorPermutation[] is based on DOI: 10.1016/j.media.2014.05.012."


ImportTracts::usage = 
"ImportTracts[file] imports a *.trk file. It can contain {tracts, vox, dim, seeds}."

ExportTracts::usage = 
"ExportTracts[file, tracts, vox, dim, seeds] exports the tracts, vox, dim and seeds to *.trk file."


(* ::Subsection::Closed:: *)
(*Options*)


FittingOrder::usage = 
"FittingOrder is an option for FitTracts. It specifies the polynomial order of the function to fit the tract."

FitTractSegments::usage =
"FitTractSegments is an option for SegmentTracts. If set True the segmented tracts are fitted with FitTracts."

TractMonitor::usage = 
"TractMonitor is an option for FiberTractography. When set True it prints the progress."

FiberLengthRange::usage = 
"FiberLengthRange is an option for FiberTractography and specifies the allowed tract range."

FiberAngle::usage = 
"FiberAngle is an option for FiberTractography and specifies the allowed angle change per tract step."

TensorFlips::usage =
"TensorFlips is an option for FiberTractography and specifies if the tensor orientation is flipped, see FlipTensorOrientation."

TensorPermutations::usage = 
"TensorPermutations is an option for FiberTractography and specifies if the tensor orientation is permuted, see FlipTensorOrientation."

StopThreshold::usage = 
"StopThreshold is an option for FiberTractography and defines the stop threshold which is a value between 0 and 1."

StepSize::usage = 
"StepSize is an option for FiberTractography and defines the tractography step size."

MaxSeedPoints::usage = 
"MaxSeedPoints is an option for FiberTractography and defines the maximum number of seed points to be used."


MaxTracts::usage = 
"MaxTracts is an option for PlotTracts. It specifies how many tracts are plotted."

TractSize::usage = 
"TractSize is an option for PlotTracts. When tubes are used it specifies the tube width."

TractColoring::usage = 
"TractColoring is an option for PlotTracts and sets how the tracts are colored. Values can be \"Direction\", \"Length\", \"Angle\", {par}, or RGBColor[].
For \"Length\", \"Angle\", {par} it can be defined in the form {..., {min, max}} where the {min, max} specifies the range of the color function."

TractReduction::usage = 
"TractReduction is an option for PlotTracts. Value can be an Integer > 0, which determines with which factor the tract coordinates are subsampled."

TractScaling::usage = 
"TractScaling is an option for PlotTracts. The value can be \"World\" or \"Voxel\", if the value is \"World\" the tracts are in mm else in voxel coordinates."

NormalizeDensity::usage = 
"NormalizeDensity is an option for TractDensityMap. If set True the tract density is normalized, if False then it is the true tract count."


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*Tract Operations*)


(* ::Subsubsection::Closed:: *)
(*FitTracts*)


Options[FitTracts] = {FittingOrder -> 5};

SyntaxInformation[FitTracts] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

FitTracts[tract_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, dim_, OptionsPattern[]]:= SelectValidCoor[FitTractsC[tract, Round@OptionValue[FittingOrder]], vox, dim]

FitTracts[tract_, OptionsPattern[]] := ToPackedArray/@FitTractsC[tract, Round@OptionValue[FittingOrder]]


FitTractsC = Compile[{{trf, _Real, 2}, {ord, _Real, 0}}, Block[{mat, r},
	r = Range[Length[trf]];
	mat = Transpose[r^# & /@ Range[0., ord]];
	mat . (PseudoInverse[mat] . trf)]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsubsection::Closed:: *)
(*SelectValidCoor*)


SyntaxInformation[SelectValidCoor] = {"ArgumentsPattern" -> {_, _, _.}};

SelectValidCoor[trF_, dim_]:=Pick[trF, SelectValidCoorC[trF, dim], 1]

SelectValidCoor[trF_, vox_, dim_]:=Pick[trF, SelectValidCoorV[trF, vox, dim], 1]


SelectValidCoorC = Compile[{{tr, _Real, 2}, {dim, _Integer, 1}}, Block[{x, y, z},
	{x, y, z} = Transpose[tr];
	UnitStep[x - 1] UnitStep[dim[[1]] - x] UnitStep[y - 1] UnitStep[dim[[2]] - y] UnitStep[z - 1] UnitStep[dim[[3]] - z]
], RuntimeOptions -> "Speed", RuntimeAttributes -> {Listable}]


SelectValidCoorV = Compile[{{tr, _Real, 2}, {vox, _Real, 1}, {dim, _Integer, 1}}, Block[{x, y, z},
	{x, y, z} = Ceiling[Transpose[tr]/vox];
	UnitStep[x - 1] UnitStep[dim[[1]] - x] UnitStep[y - 1] UnitStep[dim[[2]] - y] UnitStep[z - 1] UnitStep[dim[[3]] - z]
], RuntimeOptions -> "Speed", RuntimeAttributes -> {Listable}]


(* ::Subsubsection::Closed:: *)
(*ResampleTracts*)


SyntaxInformation[ResampleTracts] = {"ArgumentsPattern" -> {_, _}};

ResampleTracts[tracts_, len_?IntegerQ] := Block[{int, r},
	r = Rescale[N@Range[len]];
	ToPackedArray@Map[(int = ListInterpolation[Transpose[{#}], {{0, 1}}, InterpolationOrder -> 1]; int /@ r) &, tracts]
]

ResampleTracts[tracts_, len_?ListQ] := Block[{int, r},

	ToPackedArray@MapThread[(
		int = ListInterpolation[Transpose[{#1}], {{0, 1}}, InterpolationOrder -> 1];
		r = Rescale[N@Range[#2]];
		ToPackedArray[int /@ r]
	) &, {tracts, len}]
]

(* ::Subsubsection::Closed:: *)
(*MoveTracts*)


SyntaxInformation[MoveTracts] = {"ArgumentsPattern" -> {_, _}};

MoveTracts[tracts_, off:{_?NumberQ,_?NumberQ,_?NumberQ}]:=MoveTractsC[tracts, off]


MoveTractsC = Compile[{{tr, _Real, 2}, {off, _Real, 1}},
	# + off & /@ tr
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsubsection::Closed:: *)
(*RescaleTracts*)


SyntaxInformation[RescaleTracts] = {"ArgumentsPattern" -> {_, _}};

RescaleTracts[tracts_, sc_?NumberQ]:= ToPackedArray[ToPackedArray/@RescaleTractsI[tracts, {sc, sc, sc}]]

RescaleTracts[tracts_, sc:{_?NumberQ,_?NumberQ,_?NumberQ}]:=ToPackedArray[ToPackedArray/@RescaleTractsI[tracts, sc]]

(*rescale to accurate values*)
RescaleTractsI = Compile[{{tr, _Real, 2}, {sc, _Real, 1}},
	Transpose[Transpose[tr]/sc]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]

(*rescale to coordinates, for length*)
RescaleTractsC = Compile[{{tr, _Real, 2}, {sc, _Real, 1}},
	Transpose[Ceiling[Transpose[tr]/sc]]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]

(*rescale to coordinates of midpoint, for angle*)
RescaleTractsCM = Compile[{{tr, _Real, 2}, {vox, _Real, 1}},
	Transpose[Ceiling[Transpose[Mean[{tr[[2 ;; -1]], tr[[1 ;; -2]]}]]/vox]], 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]

(*rescale to coordinates dropping edges, for curvature*)
RescaleTractsCD = Compile[{{tr, _Real, 2}, {sc, _Real, 1}}, 
	Transpose[Ceiling[Transpose[tr]/sc]][[2 ;; -2]], 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsubsection::Closed:: *)
(*FilterTractLength*)


SyntaxInformation[FilterTractLength] = {"ArgumentsPattern" -> {_, _, _.}};

FilterTractLength[trksI_, {lmin_, lmax_}] := FilterTractLength[trksI, {lmin, lmax}, "tracts"]

FilterTractLength[trksI_, {lmin_, lmax_}, what_] := Block[{len, sel, drop, trks},
	len = FLengthC[trksI];
	sel = UnitStep[len - lmin];
	trks = Pick[trksI, sel, 1];
	len = Pick[len, sel, 1];
	trks = MapThread[If[#1 >= lmax, drop = Ceiling[((#1 - lmax)/Mean[Norm /@ Differences[#2]])/2] + 1; #2[[drop ;; -drop]], #2] &, {len, trks}];
	Switch[what, "tracts", trks, "sel", sel, "both", {trks, sel}]
]


(* ::Subsection::Closed:: *)
(*GetTractValues*)


Options[GetTractValues]= {InterpolationOrder->1}

SyntaxInformation[GetTractValues] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

GetTractValues[tracts_, val_, opts : OptionsPattern[]] := GetTractValues[tracts, val, {1., 1., 1.}, opts]

GetTractValues[tracts_, val_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, OptionsPattern[]]:=Block[{fun, tractsSc, int, dim},
	int = OptionValue[InterpolationOrder];

	Which[IntegerQ[int] && int >= 1,
		fun = MakeIntFunction[val, vox, OptionValue[InterpolationOrder]];
		ToPackedArray[fun @@@ #] & /@ tracts,
		True,
		tractsSc = SelectValidCoor[If[vox === {1., 1., 1.}, tracts, RescaleTractsC[tracts, vox]], Dimensions[val]];
		If[IntegerQ[int] && int === 0,
			SelectTractVal[val, tractsSc],
			SelectTractValD[val, tractsSc]
		]
	]
]


SelectTractVal = Compile[{{roi, _Real, 3}, {tract, _Integer, 2}}, 
	Part[roi, #[[1]], #[[2]], #[[3]]] & /@ tract
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];

SelectTractValD = Compile[{{roi, _Real, 3}, {tract, _Integer, 2}}, 
	Part[roi, #[[1]], #[[2]], #[[3]]] & /@ DeleteDuplicates[tract]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsection::Closed:: *)
(*SegmentTracts*)


Options[SegmentTracts] = {
	FiberLengthRange -> {15, 500}, 
	OutputForm -> "Joined", 
	FitTractSegments -> True
}

SyntaxInformation[SegmentTracts] = {"ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}};

SegmentTracts[tracts_, segs_, vox_, dim_, OptionsPattern[]] := Block[{tractsI},
	tractsI = RescaleTractsC[tracts, vox];
	tractsI = FilterTracts[tracts, tractsI, {{"and", {"partwithin", #}}}, FiberLengthRange -> OptionValue[FiberLengthRange]] & /@ Transpose[segs];

	If[OptionValue[FitTractSegments] === True, tractsI = If[Length[#] > 0, FitTracts[#, vox, dim, FittingOrder -> 3], #] & /@ tractsI];

	Switch[OptionValue[OutputForm],
		"Joined", Flatten[tractsI, 1],
		"Individual", tractsI
	]
]


(* ::Subsection:: *)
(*Tract Properties Maps*)


(* ::Subsubsection::Closed:: *)
(*SeedDensityMap*)


SyntaxInformation[SeedDensityMap] = {"ArgumentsPattern" -> {_, _, _}};

SeedDensityMap[seeds_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, dim_] := Normal@SparseArray[Normal@Counts@RescaleTractsC[seeds, vox], dim];


(* ::Subsubsection::Closed:: *)
(*TractDensityMap*)


Options[TractDensityMap] = {NormalizeDensity -> True}

SyntaxInformation[TractDensityMap] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

TractDensityMap[tracts_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, dim_,  OptionsPattern[]] := Block[{dens},
	dens = Normal@SparseArray[Select[
		Normal@Counts[Flatten[DeleteDuplicates /@ RescaleTractsC[tracts, vox], 1]],
		(1 <= #[[1, 1]] <= dim[[1]] && 1 <= #[[1, 2]] <= dim[[2]] && 1 <= #[[1, 3]] <= dim[[3]]) &],
		dim];

	ToPackedArray@N@If[OptionValue[NormalizeDensity], dens/MedianNoZero[Flatten[dens]], dens]
]


(* ::Subsubsection::Closed:: *)
(*TractLengthMap*)


SyntaxInformation[TractLengthMap] = {"ArgumentsPattern" -> {_, _, _}};

TractLengthMap[tracts_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, dim_] := GatherThread[RescaleTractsC[tracts, vox], 
	TractLength[tracts], dim]


(* ::Subsubsection::Closed:: *)
(*TractLength*)


SyntaxInformation[TractLength] = {"ArgumentsPattern" -> {_}};

TractLength[tracts_]:=FLengthC[tracts]


FLengthC = Compile[{{trc, _Real, 2}},
	Total[Norm/@Differences[trc]]
, RuntimeAttributes->{Listable},RuntimeOptions->"Speed"]


(* ::Subsubsection::Closed:: *)
(*TractAngleMap*)


SyntaxInformation[TractAngleMap] = {"ArgumentsPattern" -> {_, _, _, _.}};

TractAngleMap[tracts_, vox : {_?NumberQ, _?NumberQ, _?NumberQ}, dim_] := GatherThread[RescaleTractsCM[tracts, vox], 
	VecAngC[tracts], dim]

TractAngleMap[tracts_, v1_?VectorQ, vox : {_?NumberQ, _?NumberQ, _?NumberQ}, dim_] := GatherThread[RescaleTractsCM[tracts, vox], 
	VecAngCv1[tracts, v1], dim]

TractAngleMap[tracts_, {v1_?VectorQ, v2_?VectorQ}, vox : {_?NumberQ, _?NumberQ, _?NumberQ}, dim_] := TractAngleMap[tracts, {v1, v2, Cross[v1,v2]}, vox, dim]

TractAngleMap[tracts_, {v1_?VectorQ, v2_?VectorQ, v3_?VectorQ}, vox : {_?NumberQ, _?NumberQ, _?NumberQ}, dim_] := Block[{elevation, azimuth},
	{elevation, azimuth} = Transpose[VecAngCv2[tracts, v1, v2]];
	{GatherThread[RescaleTractsCM[tracts, vox], elevation, dim], GatherThread[RescaleTractsCM[tracts, vox], azimuth, dim]}
]


(* ::Subsubsection::Closed:: *)
(*TractAngle*)


SyntaxInformation[TractAngle] = {"ArgumentsPattern" -> {_}};

TractAngle[tracts_]:= VecAngC[tracts]

TractAngle[tracts_, v1_]:= VecAngC[tracts, v1]

TractAngle[tracts_, {v1_, v2_}]:= VecAngCv1[tracts, v1, v2]

TractAngle[tracts_, {v1_, v2_, v3_}]:= VecAngCv2[tracts, v1, v2]


VecAngC = Compile[{{tr, _Real, 2}}, Block[{angs},
	angs = ArcSin[Abs[Dot[#, {1, 0, 0}]/Norm[#]]] & /@ (tr[[2 ;; -1]] - tr[[1 ;; -2]]);
	(180./Pi) angs
], RuntimeOptions -> "Speed", RuntimeAttributes -> {Listable}];

VecAngCv1 = Compile[{{tr, _Real, 2}, {v1, _Real, 1}}, Block[{angs},
	angs = ArcSin[Abs[Dot[#, v1]/Norm[#]]] & /@ (tr[[2 ;; -1]] - tr[[1 ;; -2]]);
	(180./Pi) angs
	], RuntimeOptions -> "Speed", RuntimeAttributes -> {Listable}];

VecAngCv2 = Compile[{{tr, _Real, 2}, {v1, _Real, 1}, {v2, _Real, 1}}, Block[{vec, dotv, proj, dotp, nr, angs},
	angs = (
		(*normalize and align with plane normal*)
		vec = #/Norm[#];
		dotv = Dot[vec, v1];
		(*project, normalize and align with v2*)
		proj = vec - dotv v1;
		nr = Norm[proj];
		dotp = Dot[If[nr <= 0., v2, proj/nr], v2];
		(*calculate angles*)
		{ArcSin[Abs[dotv]], ArcCos[Abs[dotp]]}
	) & /@ (tr[[2 ;; -1]] - tr[[1 ;; -2]]);
	Transpose[(180./Pi) angs]
], RuntimeOptions -> "Speed", RuntimeAttributes -> {Listable}];


(* ::Subsubsection::Closed:: *)
(*TractCurvatureMap*)


SyntaxInformation[TractCurvatureMap] = {"ArgumentsPattern" -> {_, _, _, _.}};

TractCurvatureMap[tracts_, vox : {_?NumberQ, _?NumberQ, _?NumberQ}, dim_] := GatherThread[RescaleTractsCD[tracts, vox],  
	CurvC[tracts], dim]


(* ::Subsubsection::Closed:: *)
(*TractCurvature*)


SyntaxInformation[TractCurvature] = {"ArgumentsPattern" -> {_}};

TractCurvature[tracts_]:= CurvC[tracts]


CurvC = Compile[{{tract, _Real, 2}}, Block[{diff, ds, d1, d2},
	diff = tract[[2 ;; -1]] - tract[[1 ;; -2]];
	ds = Norm /@ diff;
	d1 = diff / ds;
	d2 = (d1[[2 ;; -1]] - d1[[1 ;; -2]]) / Mean[{ds[[2 ;; -1]], ds[[1 ;; -2]]}];
	1000. Norm[Cross[#[[1]], #[[2]]]] & /@ Transpose[{d1[[;; -2]], d2}]
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsubsection::Closed:: *)
(*GatherThread*)


GatherThread[coor_?ListQ, val_?ListQ, dim_?VectorQ] := Block[{len, ran, list, out},
	(*Split in sub lists to prevent memory crash for large tracts sets*)
	len = Length@coor;
	out = If[len <= 50000,
		GatherThread[coor, val],
		ran = DeleteDuplicates[Append[Range[0, len, 50000], len]];
		list = Thread[{ran[[;; -2]] + 1, ran[[2 ;;]]}];
		GatherThread[Table[GatherThread[coor[[i[[1]] ;; i[[2]]]], val[[i[[1]] ;; i[[2]]]]], {i, list}]]
	];

	(*clip coordinates and calculate median*)
	out[[All, 1]] = Transpose[{
		Clip[out[[All, 1, 1]], {1, dim[[1]]}], 
		Clip[out[[All, 1, 2]], {1, dim[[2]]}], 
		Clip[out[[All, 1, 3]], {1, dim[[3]]}]
	}];
	out[[All, 2]] = N[Median /@ out[[All, 2]]];

	(*matrix from coordinate rule*)
	ToPackedArray@Normal@SparseArray[out, dim, 0.]
]

GatherThread[coor_, val_] := GatherThread[Thread /@ Thread[coor -> val]]

GatherThread[rule_] := Block[{out},
	out = GatherBy[Flatten[rule, 1], First];
	Thread[out[[All, 1, 1]] -> (Flatten /@ out[[All, All, 2]])]
]


(* ::Subsection:: *)
(*FiberTractography*)


(* ::Subsubsection::Closed:: *)
(*FiberTractography*)


Options[FiberTractography] = {
	FiberLengthRange -> {20, 500},
	FiberAngle -> 30,
	TensorFlips -> {1, 1, 1},
	TensorPermutations -> {"x", "y", "z"},
	InterpolationOrder -> 0,
	StopThreshold -> 0.5,
	StepSize -> Automatic,
	Method -> "RK4",
	MaxSeedPoints -> Automatic,
	TractMonitor -> True,
	Parallelization -> False
};

SyntaxInformation[FiberTractography] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

FiberTractography[tensor_, voxi_, opts : OptionsPattern[]] := FiberTractography[tensor, voxi, {0. First@tensor, {0., 1.}}, opts]

FiberTractography[tensor_, voxi_, {par_?ArrayQ, {min_?NumberQ, max_?NumberQ}}, opts : OptionsPattern[]] := FiberTractography[tensor, voxi, {{par, {min, max}}}, opts]

FiberTractography[tensor_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, inp : {{_, {_, _}} ...}, OptionsPattern[]] := Block[{
		lmin, lmax, amax, maxSeed, flip, per, int, stopT, step, tracF, vecF, trFunc, ran,
		tens, tensMask, inpTr, treshs, stop, coors, vecInt, stopInt, ones, dim, crp, t2,
		seedN, seedI, seedT, seeds, t1, tracts, iii, drop, smax, len, sel, mon
	},

	(*get the options*)
	{{lmin, lmax}, amax, maxSeed, flip, per, int, stopT, step, mon} = OptionValue[{
		FiberLengthRange, FiberAngle, MaxSeedPoints, TensorFlips, TensorPermutations, 
		InterpolationOrder, StopThreshold, StepSize, TractMonitor}];

	step = N@If[NumberQ[step],step, Min[0.5 vox]];
	smax = Ceiling[(lmax/step)];
	tracF = Switch[OptionValue[Method], "RungeKutta" | "RK" | "RK2", RK2, "RungeKutta4" | "RK4", RK4, _, Euler];
	
	(*prepare tensor and stop data, remove background for performance and make int functions*)
	dim = Rest@Dimensions@tensor;
	{stop, crp} = AutoCropData[Unitize[First@tensor] Mask[inp], CropPadding -> 0];
	stop = ToPackedArray@N@Normal@stop;

	(*for 0th order interpolation tensor interpolation is not needed, vector field can be precalculated*)
	tens = ApplyCrop[#, crp]& /@ FlipTensorOrientation[tensor, per, flip];
	{vecF, tens} = If[int>=1, 
		{EigVec, ToPackedArray@N@RotateDimensionsLeft[{tens}, 2]},
		{Vec, ToPackedArray@N@RotateDimensionsLeft[{RotateDimensionsRight[EigVec[RotateDimensionsLeft@tens, {1, 0, 0}]]}, 2]}
	];

	(*make the random seed points*)
	SeedRandom[1234];
	seedN = Total@Flatten@stop;
	maxSeed = Round@Which[
		NumberQ[maxSeed], maxSeed, 
		maxSeed === Automatic, 0.3 seedN,
		maxSeed === All, seedN,
		Head[maxSeed]===Scaled, First[maxSeed] seedN
	];

	(*random seedpoint selection until maxSeed is found*)
	seeds = {};
	seedI = SparseArray[stop]["NonzeroPositions"];
	stopInt = MakeIntFunction[stop, vox, int];
	While[Length[seeds] < maxSeed,
		seedT = Flatten[RandomSample[seedI, #] & /@ Block[{i = maxSeed},
			First@Last@Reap[While[i > 0, Sow[If[i > seedN, seedN, i]];
			i = i - seedN;]]], 1];
		seedT = # vox & /@ (seedT + RandomReal[{-0.999, 0}, {maxSeed, 3}]);
		seedT = Pick[seedT, stopInt @@@ seedT, 1.];		
		seeds = Join[seeds, seedT];
		seeds = seeds[[;; Min[{Length@seeds, maxSeed}]]]
	];

	(*finalize seed points*)
	seeds = ToPackedArray@N@seeds;
	seedN = Length@seeds;

	(*start the tractography*)
	If[mon,
		Echo["Starting with "<>ToString[seedN]<>" seed points and stepsize "<>ToString[step]<>" mm"];
	];

	(*check if parallel or normal computing is needed*)
	If[OptionValue[Parallelization],
		(*parallel tractography preparation*)
		If[mon, Echo["Starting parallel preparation"]];
		ran = Thread[{vox/2, vox Dimensions[stop] - vox/2}];
		t2 = First@AbsoluteTiming[
			DistributeDefinitions[step, amax, smax, stopT, tracF, vecF, seeds, tens, stop, vox, int, ran,
				TractFunc, TractFuncI, EigVec, Vec, VecAng, Euler, RK2, RK4];
			ParallelEvaluate[
				vecInt = ListInterpolation[tens, ran, InterpolationOrder -> int, 
					"ExtrapolationHandler" -> {(0. tens[[1, 1, 1, 1]] &), "WarningMessage" -> False}];
				stopInt = ListInterpolation[stop, ran, InterpolationOrder -> int, 
					"ExtrapolationHandler" -> {(0. &), "WarningMessage" -> False}];
				trFunc = TractFunc[#, step, {amax, smax, stopT}, {vecInt, stopInt, tracF, vecF}]&;
			, DistributedContexts -> None]];
		If[mon, Echo["Parallel preparation time: "<>ToString[Round[t2, .1]]<>" seconds"]];
		(*actual tracto for parallel with memory clear*)
		{t1, tracts} = AbsoluteTiming@ParallelMap[trFunc, seeds, 
			Method -> "EvaluationsPerKernel" -> 10, ProgressReporting -> mon];
		ParallelEvaluate[Clear[vecInt, stopInt, trFunc, tens, stop]];
		ClearDistributedDefinitions[];
		,
		(*normal tractography with monitoring*)
		vecInt = MakeIntFunction[tens, vox, int];
		stopInt = MakeIntFunction[stop, vox, int];
		trFunc = TractFunc[#, step, {amax, smax, stopT}, {vecInt, stopInt, tracF, vecF}]&;
		{t1, tracts} = AbsoluteTiming@If[mon,
			Monitor[Table[trFunc@seeds[[iii]], {iii, 1, seedN}], ProgressIndicator[iii, {0, seedN}]],
			Table[trFunc@seeds[[iii]], {iii, 1, seedN}]
		];
	];

	(*select only tracts within correct range and clip tracts that are longer*)
	If[mon, Echo["Checking tract lengths"]];
	{tracts, sel} = FilterTractLength[tracts, {lmin, lmax}, "both"];
	seeds = Pick[seeds, sel, 1];

	(*report timing an results*)
	If[mon,
		Echo["Tractography took "<>ToString[Round[t1,.1]]<>" seconds ("<>ToString[Round[seedN/t1]]<>" tracts/s)"];
		Echo[ToString[Length[tracts]]<>" valid tracts with length "<>ToString[Round[step Mean[Length /@ tracts], 0.1]]<>"\[PlusMinus]"<>ToString[Round[step StandardDeviation[Length /@ tracts], 0.1]]<>" mm"];	
	];

	(*output tracts move them to coordinates before cropping*)
	MoveTracts[#, vox (crp[[;; ;; 2]] - 1)]&/@{tracts, seeds}
]


(* ::Subsubsection::Closed:: *)
(*TractFunc*)


TractFunc[loc0_?VectorQ, h_, stp_, fun_] := Block[{dir0, mh = 0.75 h, out},
	(*startpoint*)
	dir0 = h Last[fun][First[fun]@@loc0, {0,0,1}];
	(*only tract if step is actuall step size*)
	If[Norm[dir0] <= mh, {loc0},
		(*tracts go bi-directional from start with half step offset*)
		out = ToPackedArray@Join[
			Reverse@TractFuncI[{loc0 + dir0/2, dir0}, {h, mh}, stp, fun], 
			TractFuncI[{loc0 - dir0/2, -dir0}, {h, mh}, stp, fun]];
		If[out==={}, {loc0}, out]
	]	
]


TractFuncI[{loci_?VectorQ, stepi_?VectorQ}, {h_, mh_}, {amax_, smax_, stoptr_}, {vecInt_, stopInt_, tracF_, vecF_}] := Block[{loc1, step0, step1, loc2},
	loc1 = loci;(*current point*)
	step0 = stepi;(*incomming direction*)
	(*break if start point is not valid*)
	If[stopInt @@ loc1 < stoptr, Return[{}]];

	(*perform the tractography while conditions are valid*)
	Flatten[Last[Reap[Do[
		(*sow the location and update steps*)
		Sow[loc1];
		step1 = tracF[loc1, step0, h, vecInt, vecF];
		loc2 = loc1 + step1;
		(*check for stop*)
		If[Norm[step1] < mh || VecAng[step0, step1] > amax || (stopInt @@ loc2) < stoptr, Break[]];
		(*update points for next step*)
		loc1 = loc2;
		step0 = step1;
	, smax]]], 1]
];


(* ::Subsubsection::Closed:: *)
(*Euler*)


Euler[y_, v_, h_, int_, vec_] := Block[{k1},
	k1 = h vec[int@@(y), v];
	k1
]


(* ::Subsubsection::Closed:: *)
(*RK2*)


RK2[y_, v_, h_, int_, vec_] := Block[{k1, k2},
	k1 = h vec[int@@(y), v];
	k2 = h vec[int@@(y + k1/2), v];
	k2	
]


(* ::Subsubsection::Closed:: *)
(*RK4*)


RK4[y_, v_, h_, int_, vec_] := Block[{k1, k2, k3, k4},
	k1 = h vec[int@@(y), v];
	k2 = h vec[int@@(y + k1/2), v];
	k3 = h vec[int@@(y + k2/2), v];
	k4 = h vec[int@@(y + k3), v];
	k1/6 + k2/3 + k3/3 + k4/6
]


(* ::Subsubsection::Closed:: *)
(*VecAng*)


VecAng = Compile[{{v1, _Real, 1}, {v2, _Real, 1}}, Block[{v, n1, n2},
	(*norm and angles*)
	n1 = Norm[v1];
	n2 = Norm[v2];
	(*if one of the two vectors has norm 0. output 90 degrees*)
	If[n1 === 0. || n2 === 0.,90.,
		(*normalize and constrain before calculating the angle*)
		180./Pi ArcCos[Min[1., Max[-1., Dot[v1, v2]/(n1 n2)]]]
	]
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*EigVec*)


(*Solves the cubic equation analytically using the characteristic polynomial. This is an exact method.*)
(*other possible methods could have been power iterations or newtons methods, both are itterative and are not faster than the analytical*)
EigVec = Compile[{{tens, _Real, 1}, {vdir, _Real, 1}}, Block[{
		dxx, dyy, dzz, dxy, dxz, dyz, dxy2, dxz2, dyz2, 
		i1, i2, i3, i, v, s, v2, vv2, l1, a, b, c, norm, vec
	},
	(*tens = 1000 t;*)
	(*method https://doi.org/10.1016/j.mri.2009.10.001*)
	vec = If[Total[Abs[tens]]<10.^-15, {0., 0., 0.},

		(* Extract tensor components *)
		{dxx, dyy, dzz, dxy, dxz, dyz} = tens;
		{dxy2, dxz2, dyz2} = {dxy, dxz, dyz}^2;

		(* First and second invariants *)
		i1 = dxx + dyy + dzz;
		i2 = dxx dyy + dxx dzz + dyy dzz - dxy2 - dxz2 - dyz2;
		i3 = dxx dyy dzz + 2 dxy dxz dyz - dzz dxy2 - dyy dxz2 - dxx dyz2;

		(* Calculate the first eigenvalue using an approximation *)
		i = i1/3;
		v = i^2 - i2/3;
		(*for v<=0 l1 does not make sence for DTI*)
		If[v <= 0., {0., 0., 0.},
			(* Use trigonometric solution for the largest eigenvalue *)
			v2 = Sqrt[v];
			s = i^3 - (i1 i2)/6 + i3/2;
			l1 = i + 2 v2 Cos[ArcCos[Min[1., Max[-1., s/(v v2)]]]/3];

			(* Calculate the corresponding eigenvector components*)
			{a, b, c} = {dxz dxy, dxy dyz, dxz dyz} - {dyz, dxz, dxy} ({dxx, dyy, dzz} - l1);
			vec = {b c, a c, a b};

			(*normalize the vector*)
			norm = Norm[vec];
			If[norm < 10.^-30 || l1 < 0., {0., 0., 0.}, 
				vec / norm
			]
		]
	];

	(*align the vector with the tracts*)
	Sign[Sign[Dot[vdir, vec]] + 0.1] vec
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsubsection::Closed:: *)
(*Vec*)


Vec = Compile[{{vec, _Real, 1}, {vdir, _Real, 1}},
	Sign[Sign[Dot[vdir, vec]] + 0.1] vec
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsection::Closed:: *)
(*FindTensorPermutation*)


Options[FindTensorPermutation] = {
	FiberLengthRange -> {20, 500},
	FiberAngle -> 30,
	InterpolationOrder -> 0,
	StopThreshold -> 0.5,
	StepSize -> 2,
	Method -> "Euler",
	MaxSeedPoints -> 500
};

SyntaxInformation[FindTensorPermutation] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

FindTensorPermutation[tensor_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, opts : OptionsPattern[]] := FindTensorPermutation[tensor, vox, {0. First@tensor, {0., 1.}}, opts]

FindTensorPermutation[tensor_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, {par_?ArrayQ, {min_?NumberQ, max_?NumberQ}}, opts : OptionsPattern[]] := FindTensorPermutation[tensor, vox, {{par, {min, max}}}, opts]

FindTensorPermutation[tens_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, stop_, opts : OptionsPattern[]] := Block[{
	step, tracts, perms, flips, lengs, i, j, pl, ind
	},

	step = OptionValue[StepSize];
	step = N@If[NumberQ[step],step, Min[0.5 vox]];

	perms = {{"x", "y", "z"}, {"x", "z", "y"}, {"y", "x", "z"}, {"y","z", "x"}, {"z", "x", "y"}, {"z", "y", "x"}};
	flips = {{1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1}};
	ind = {1, 1};

	PrintTemporary[Dynamic[{ind, flips[[ind[[1]]]], perms[[ind[[2]]]]}]];
	lengs = Table[
		ind = {i, j};
		tracts = First@FiberTractography[tens, vox, stop, TractMonitor -> False,TensorFlips -> flips[[i]], TensorPermutations -> perms[[j]],
			(Sequence[# -> OptionValue[#] & /@ FilterRules[Options[FiberTractography], Options[FindTensorPermutation]][[All, 1]]])
		];
		N@Mean[Length /@ tracts] step
	, {i, 1, 4}, {j, 1, 6}];

	pl = ArrayPlot[lengs, ColorFunction -> "Rainbow", ImageSize -> 150, Frame -> True, 
		FrameTicks -> {{Thread[{Range[4], flips}], None},{None, Thread[{Range[6], Rotate[#, 90 Degree] & /@ perms}]}}];

	{i, j} = FirstPosition[lengs, Max[lengs]];
	{perms[[j]], flips[[i]], pl}
]


(* ::Subsection:: *)
(*ROI Filter Tracts *)


(* ::Subsubsection::Closed:: *)
(*FilterTracts*)


Options[FilterTracts] = {FiberLengthRange -> {20, 500}}

SyntaxInformation[FilterTracts] = {"ArgumentsPattern" -> {_, _, _, _., OptionsPattern[]}};

FilterTracts[tracts_, tractsI_, select_, opts:OptionsPattern[]]:=FilterTracts[tracts, tractsI, None, select, opts]

FilterTracts[tracts_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, select_, opts:OptionsPattern[]]:= FilterTracts[tracts, tracts, vox, select, opts]

FilterTracts[tracts_, tractsI_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}|None, select_, opts:OptionsPattern[]] := FilterTractLength[
	SelectTracts[tracts,
		CombineROIs[{#[[1]], Switch[ToLowerCase[#[[2, 1]]],
			"through", SelectTractTroughVol,
			"within", SelectTractInVol,
			"partwithin", SelectTractPartInVol,
			"x" | "y" | "z", SelectTractTroughPlane
		][tractsI, If[StringLength[#[[2, 1]]] == 1, #[[2]], #[[2, 2]]], vox]} & /@ select]
	]
	, OptionValue[FiberLengthRange]
]


(* ::Subsubsection::Closed:: *)
(*SelectTracts*)


SelectTracts[tracts_, sel_] := ToPackedArray/@PartTracts[DeleteCases[Pick[tracts, Clip[Round@sel, {0, 1}], 1], {}]]


(* ::Subsubsection::Closed:: *)
(*PartTracts*)


PartTracts[tracts_] := Block[{tmp, sts, st}, 
	Select[Flatten[(
		tmp = #;
		sts = 1.1 Min[st = Norm /@ Differences[tmp]];
		tmp[[#[[1]] + 1 ;; #[[2]]]] & /@ Partition[Flatten[{0, Position[UnitStep[st - sts], 1], Length[tmp]}], 2, 1]
	) & /@ tracts, 1], Length[#] > 3 &]
]


(* ::Subsubsection::Closed:: *)
(*CombineROIs*)


CombineROIs[rois : {{_?StringQ, _?ListQ} ..}] := Block[{or, and, not, test},
	or = Select[rois, (ToLowerCase[#[[1]]]==="or") &][[All, 2]];
	or = If[or === {}, 1, Clip[Total[or], {0, 1}]];

	not = Select[rois, (ToLowerCase[#[[1]]]==="not") &][[All, 2]];
	not = 1 - If[not === {}, 0, Clip[Total[not], {0, 1}]];

	and = Select[rois, (ToLowerCase[#[[1]]]==="and") &][[All, 2]];
	and = If[and === {}, 1, Clip[Times @@ and, {0, 1}]];

	Times @@ {or, and, not}
]


(* ::Subsubsection::Closed:: *)
(*SelectTractTroughPlane*)


SelectTractTroughPlane[tracts_, {dir_, slice_}, None]:= SelectTractTroughPlane[tracts, {dir, slice}, {1, 1, 1}]

SelectTractTroughPlane[tracts_, {dir_, slice_}, vox:{_?NumberQ,_?NumberQ,_?NumberQ}]:= SelectTractTroughPlaneV@@Switch[ToLowerCase[dir],
	"x", {tracts[[All, All, 3]], slice vox[[3]] - 0.5},
	"y", {tracts[[All, All, 2]], slice vox[[2]] - 0.5},
	"z", {tracts[[All, All, 1]], slice vox[[1]] - 0.5}
];


SelectTractTroughPlaneV = Compile[{{tract, _Real, 1}, {val, _Real, 0}},
	Boole[0 < Total[UnitStep[tract - val]] < Length[tract]], 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*SelectTractTroughVol*)


SelectTractTroughVol[tracts_, roi_, None]:= SelectTractTroughVolC[Normal@roi, tracts]

SelectTractTroughVol[tracts_, roi_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}]:= SelectTractTroughVolV[Normal@roi, tracts, vox]


SelectTractTroughVolV = Compile[{{roi, _Integer, 3}, {tract, _Integer, 2}},
	Max[Part[roi, #[[1]], #[[2]], #[[3]]] & /@ tract], 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


SelectTractTroughVolV = Compile[{{roi, _Integer, 3}, {tract, _Real, 2}, {vox, _Real, 1}},
	Max[Part[roi, #[[1]], #[[2]], #[[3]]] & /@ Transpose[Ceiling[Transpose[tract]/vox]]], 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*SelectTractInVol*)


SelectTractInVol[tracts_, roi_, None]:= SelectTractInVolC[Normal@roi, tracts]

SelectTractInVol[tracts_, roi_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}]:= SelectTractInVolV[Normal@roi, tracts, vox]


SelectTractInVolC = Compile[{{roi, _Integer, 3}, {tract, _Integer, 2}},
	Floor[Total[Part[roi, #[[1]], #[[2]], #[[3]]] & /@ tract] / Length[tract]],
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


SelectTractInVolV = Compile[{{roi, _Integer, 3}, {tract, _Real, 2}, {vox, _Real, 1}},
	Floor[Total[Part[roi, #[[1]], #[[2]], #[[3]]] & /@ Transpose[Ceiling[Transpose[tract]/vox]]]/Length[tract]],
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*SelectTractPartInVol*)


SelectTractPartInVol[tracts_, roi_, None]:= SelectTractPartInVolC[Normal@roi, tracts]

SelectTractPartInVol[tracts_, roi_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}]:= SelectTractPartInVolV[Normal@roi, tracts, vox]


SelectTractPartInVolC = Compile[{{roi, _Integer, 3}, {tract, _Integer, 2}},
	Part[roi, #[[1]], #[[2]], #[[3]]] & /@ tract
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


SelectTractPartInVolV = Compile[{{roi, _Integer, 3}, {tract, _Real, 2}, {vox, _Real, 1}},
	Part[roi, #[[1]], #[[2]], #[[3]]] & /@ Transpose[Ceiling[Transpose[tract]/vox]]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsection:: *)
(*PlotTracts*)


(* ::Subsubsection::Closed:: *)
(*PlotTracts*)


Options[PlotTracts] = {
	MaxTracts -> 2000, 
	ImageSize -> 600, 
	Method -> "line", 
	TractColoring -> "Direction",
	ColorFunction -> "SouthwestColors",
	Boxed -> True,
	TractSize -> 1,
	TractReduction -> 3,
	TractScaling -> "World",
	PerformanceGoal -> "Quality"
}

SyntaxInformation[PlotTracts] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

PlotTracts[tracts_, voxi_, opts : OptionsPattern[]] := PlotTracts[tracts, voxi, 0, opts]

PlotTracts[tracts_, voxi_, dimi_, OptionsPattern[]] := Block[{
		range, vox, size, select, opts, col, tube, line, plot, colOpt, 
		met, pran, max, colf, colArray, sc, n, scale, red, qual
	},

	(*reduce points along tracts*)
	(*tracts = Select[tractsI[[All,1;;-1;;OptionValue[TractReduction]]],Length[#]>3&];*)
	(*Graphics`RenderTiming*)

	(*select correct number of tracts*)
	max = OptionValue[MaxTracts];
	(*If[OptionValue[Method]==="tube", max =  Min[4000, max]];*)
	SeedRandom[1234];
	select = ToPackedArray/@RandomSample[tracts, Min[max, Length[tracts]]];

	red = OptionValue[TractReduction];
	n = If[IntegerQ[red],
		Max[{#, 3}]&/@Ceiling[Length/@select / red],
		Round[Quantile[Length/@select, 0.75] / First@red]
	];
	select = ResampleTracts[select, n];

	(*calculated needed sizes ranges and scales*)
	scale = OptionValue[TractScaling];
	vox = Switch[scale, "World", {1,1,1}, _, Reverse@voxi];
	range = If[dimi === 0,
		Round[Reverse[MinMax /@ Transpose@Flatten[tracts, 1]]/vox],
		Reverse@Thread[{{0, 0, 0}, Switch[scale, "World", voxi dimi, _, dimi]}]];
	size = vox Flatten[Differences /@ range];
	sc = OptionValue[TractSize];

	(*get the tract vertex colors*)
	colf = OptionValue[ColorFunction];
	colOpt = OptionValue[TractColoring];

	{met, pran} = Which[
		Head[colOpt]===MaterialShading, {"Material", colOpt},
		StringQ[colOpt]||ArrayQ[colOpt], {colOpt, Automatic}, 
		Length[colOpt]===2, colOpt, 
		MemberQ[{Hue, RGBColor, GrayLevel, CMYKColor, LABColor, LCHColor, LUVColor, XYZColor}, Head[colOpt]], {"Color", RGBColor@colOpt},
		True, {"Direction", Automatic}
	];
	If[ArrayQ[met], colArray=met; met="Array"];

	col = ToPackedArray@Switch[met,
		"Direction",
		MakeDirectionColor[select],
		"Length",
		MakeLengthColor[select, {pran, colf}],
		"Angle",
		MakeAngleColor[select, {pran, colf}],
		"Array",
		MakeArrayColor[select, {pran, colf}, colArray, voxi],
		"Color",
		MakeConstantColor[select, pran],
		"Material",
		None,
		_,
		White
	];

	(*make the plot*)
	qual = OptionValue[PerformanceGoal];
	opts = Sequence[{
		Method -> {"TubePoints" -> {Which[qual==="Quality", 7, qual==="Speed", 5, IntegerQ[qual], qual, True, 5], 2}}, 
		ImageSize -> OptionValue[ImageSize], 
		Boxed -> OptionValue[Boxed],
		Axes -> OptionValue[Boxed],
		BoxRatios -> size,
		PlotRange -> range,
		SphericalRegion -> True, Lighting -> "ThreePoint", 
		Background -> Lighter@Gray,  LabelStyle -> Directive[{Bold, 16, White}]
	}];

	select = Reverse[select, 3];
	plot = Graphics3D[Switch[OptionValue[Method],
		"tube", 
		{If[met==="Material", pran, Nothing], CapForm["Butt"], JoinForm["Miter"], 
			Switch[scale, 
				"World", Tube[select, sc, VertexColors -> col],
				_, Scale[Tube[select, sc, VertexColors -> col], 1/vox, {0,0,0}]
			]
		},
		"line", 
		Switch[scale, 
			"World", Line[select, VertexColors -> col], 
			_, Line[RescaleTracts[select, vox], VertexColors -> col]
		],
		_, $Failed
	], opts]
]


(* ::Subsubsection::Closed:: *)
(*MakeDirectionColor*)


MakeDirectionColor[tract : {{_?NumberQ, _?NumberQ, _?NumberQ} ..}] := Block[{dirs},
	dirs = Abs[Differences[tract]];
	ToPackedArray[Reverse[Normalize[#]]] & /@ Mean[{Prepend[dirs, dirs[[1]]], Append[dirs, dirs[[-1]]]}]
]

MakeDirectionColor[tracts : {_?ListQ ..}] := ToPackedArray[MakeDirectionColor /@ tracts]


(* ::Subsubsection::Closed:: *)
(*MakeLengthColor*)


MakeLengthColor[tracts_, {pran_, colf_}]:=Block[{len, col},
	len = FLengthC[tracts];
	len = Rescale[len, If[pran === Automatic, Quantile[len, {.05, 0.95}], pran]];
	col = ColorData[colf];

	ToPackedArray@MapThread[ToPackedArray[ConstantArray[#2 /. RGBColor -> List, Length@#1]]&, {tracts, col /@ len}]
];


(* ::Subsubsection::Closed:: *)
(*MakeAngleColor*)


MakeAngleColor[tracts_, {pran_, colf_}] := Block[{ang, col},
	ang = (Mean[{Prepend[#, #[[1]]], Append[#, #[[-1]]]}]) & /@VecAngC[tracts];
	ang = Rescale[ang, If[pran === Automatic, {0, 90}, pran]];
	col = ColorData[colf];

	ToPackedArray[ToPackedArray[(col /@ #) /. RGBColor -> List] & /@ ang]
];


(* ::Subsubsection::Closed:: *)
(*MakeArrayColor*)


MakeArrayColor[tract_, {pran_, colf_}, dat_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}] := Block[{vals, col},
	vals = GetTractValues[tract, dat, vox, InterpolationOrder -> 0];
	vals = Rescale[vals, If[pran === Automatic, Quantile[Flatten[vals], {.05, 0.95}], pran]];
	col = ColorData[colf];

	ToPackedArray[ToPackedArray[(col /@ #) /. RGBColor -> List] & /@ vals]
]


(* ::Subsubsection::Closed:: *)
(*MakeConstantColor*)


MakeConstantColor[tract_, pran_] := MakeConstantColor[tract, pran, True]

MakeConstantColor[tract_, pran_, rand_] := Block[{vals, col},
	ToPackedArray[If[rand,
		col = Table[Blend[{Darker[pran,.2], pran, Lighter[pran,.2]}, x], {x, 0., 1., 1./10}] /. RGBColor -> List;
		ToPackedArray[ConstantArray[RandomChoice[col], Length@#1]] & /@ tract
		,
		col = pran /. RGBColor -> List;
		ToPackedArray[ConstantArray[col, Length@#1]] & /@ tract
	]]
]


(* ::Subsection::Closed:: *)
(*PlotSegmentedTracts*)


Options[PlotSegmentedTracts] := {
	MaxTracts -> 5000,
	FiberLengthRange -> {20, 500},
	Method -> "line",
	OutputForm -> "All",
	ImageSize->400,
	Monitor -> False,
	TractSize -> 1,
	ColorFunction -> "RomaO",
	PerformanceGoal -> "Quality",
	ContourOpacity -> 0.3
}

SyntaxInformation[PlotSegmentedTracts] = {"ArgumentsPattern" -> {_, _, _, _, _., OptionsPattern[]}};

PlotSegmentedTracts[tracts_, segments_, dim_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, opts : OptionsPattern[]] := PlotSegmentedTracts[tracts, segments, None, dim, vox, opts]

PlotSegmentedTracts[tracts_, segments_, bones_, dim_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, opts : OptionsPattern[]] := Block[{
		ntr, fran, type, segs, tractsF, tractsFI, ran, rand, colListT, colListC, showF,
		ref, bon, musc, tract, tracksSel, lengs, nTracts, sel, mon, size, output, colFunc
	},

	(*get options*)
	{ntr, fran, type, mon, size, output, colFunc, sizeT, qual, opa} = OptionValue[{
		MaxTracts, FiberLengthRange, Method, Monitor, ImageSize, OutputForm, 
		ColorFunction, TractSize, PerformanceGoal, ContourOpacity}];

	(*prepare data*)
	segs = Transpose@segments;
	If[mon, Echo["Fitting tracts"]];
	SeedRandom[1234];
	tractsF = RandomSample[tracts, Min[{5 ntr, Length@tracts}]];
	tractsFI = RescaleTractsC[tractsF, vox];

	(*make colors*)
	SeedRandom[1234];
	ran = Range[Length@segs];
	If[Head[colFunc] === MaterialShading,
		colListC = White;
		colListT = If[ListQ[colFunc[[1]]]&&AllTrue[colFunc[[1]], StringQ], 
			MaterialShading/@RandomChoice[colFunc[[1]], Length@segs], 
			ConstantArray[colFunc, Length@segs]
		];
		qual = "Speed";
		,
		colListT = colListC = Reverse[ColorData[colFunc] /@ Rescale[ran]][[RandomSample[ran]]];
	];

	(*reference environement*)
	If[mon, Echo["Making muscle iso volumes"]];
	ref = PlotTracts[tractsF, vox, dim, MaxTracts -> 1, Method -> "line", TractColoring -> RGBColor[0, 0, 0, 0]];
	ref[[1]] = {};

	(*make the muscle contours*)
	musc = If[colListC =!= None, Table[PlotContour[segs[[i]], vox, ContourOpacity -> opa, ContourColor -> If[ColorQ[colListC],colListC,colListC[[i]]], 
		ContourSmoothRadius -> 2, ContourResolution -> 2], {i, ran}], Graphics3D[]];
	bon = If[bones =!= None, PlotContour[bones, vox, ContourOpacity -> 1, ContourColor -> Lighter@Gray, 
		ContourSmoothRadius -> 2, ContourResolution -> 2], Graphics3D[]];

	If[mon, Echo["Making per muscle tracts"]];
	(*select the tracts per muscle and make fiber plots*)
	tracksSel = FilterTracts[tractsF, tractsFI, {{"and", {"partwithin", #}}}, FiberLengthRange -> fran] & /@ segs;
	tracksSel = If[#=!={}, FitTracts[#, vox, dim, FittingOrder -> 3], {}]& /@ tracksSel;

	lengs = Length /@ tracksSel;
	nTracts = Round[ntr lengs/Total[lengs]];
	sel = UnitStep[nTracts - 11];
	tract = MapThread[If[#4 =!= 1, Graphics3D[],
		PlotTracts[#1, vox, dim, MaxTracts -> #2, Method -> type, TractSize -> sizeT, TractColoring -> #3, 
		TractReduction -> 5, PerformanceGoal -> "Speed"]
	] &, {tracksSel, nTracts, colListT, sel}];

	If[mon, Echo["Finalizing scenes"]];
	showF = Show[ref, ##, ImageSize -> size, Axes -> False, Boxed -> False, ViewPoint -> {0., -1.5, 0.5}, 
		BaseStyle -> RenderingOptions -> {"3DRenderingMethod" -> If[qual==="Speed",  "HardwareDepthBuffer", Automatic]}
	] & @@ # &;

	Switch[output,
		"All", showF[{bon, tract, musc}],
		"Groups", showF[{#}] & /@ {bon, tract, musc},
		"Joined", showF[{bon, #}] & /@ Thread[{tract, musc}],
		"Individual", {showF[{bon}], showF[{#}] & /@ tract, showF[{#}] & /@ musc}
	]
]


(* ::Subsection::Closed:: *)
(*ImportTracts*)


(* ::Subsubsection::Closed:: *)
(*RegisterImport*)


ImportExport`RegisterImport["trk", 
	ImportTractsDefault, {}, 
	"AvailableElements" -> {"Tracts", "VoxelSize", "Dimensions", "Seeds"}, 
	"OriginalChannel" -> True
];


(* ::Subsubsection::Closed:: *)
(*ImportTracts*)


SyntaxInformation[ImportTracts] = {"ArgumentsPattern" -> {_}};

ImportTracts[] := ImportTracts[""]

ImportTracts[file_] := Block[{fileI},
	fileI = If[file == "",
		FileSelect["FileOpen", {"*.trk"}, "trk", WindowTitle -> "Select the trk file to import"],
		If[FileExistsQ[file], file, $Failed]
	];
	Import[fileI, "trk"]
]


(* ::Subsubsection::Closed:: *)
(*ImportTractsDefault*)


ImportTractsDefault[file_, ___] := Block[{strm, all, nTr, nTrLeng, dim, vox, seeds, tracts},

	strm = OpenRead[file, BinaryFormat -> True];

	(*nDim, nVox, nSeed, nTrCoor*)
	all = BinaryReadList[strm, "Integer32", 4];

	(*the number of tracts and each tracts length*)
	nTr = First@BinaryReadList[strm, "Integer32", 1];
	nTrLeng = BinaryReadList[strm, "Integer32", nTr];

	(*read the data*)
	{dim, vox, seeds, tracts} = DynamicPartition[BinaryReadList[strm, "Real32", Total[all]], all];

	(*partition the seeds*)
	seeds = Partition[seeds, 3];
	tracts = DynamicPartition[Partition[tracts, 3], nTrLeng];
	Close[strm];

	(*give the output*)
	{tracts, vox, Round[dim], seeds}
]


(* ::Subsection::Closed:: *)
(*ExportTracts*)


(* ::Subsubsection::Closed:: *)
(*RegisterExport*)


ImportExport`RegisterExport["trk",
	ExportTractsDefault,
	"DefaultElement" -> {"Tracts", "VoxelSize", "Dimensions", "Seeds"},
	"AvailableElements" -> {"Tracts", "VoxelSize", "Dimensions", "Seeds"},
	"OriginalChannel" -> True
];


(* ::Subsubsection::Closed:: *)
(*ExportTracts*)


SyntaxInformation[ExportTracts] = {"ArgumentsPattern" -> {_, _., _., _., _.}};

ExportTracts[tracts : {_?ListQ ..}] := ExportTracts["", tracts, {0, 0, 0}, {0, 0, 0}, {}]

ExportTracts[tracts : {_?ListQ ..}] := ExportTracts["", tracts, {0, 0, 0}, {0, 0, 0}, {}]

ExportTracts[tracts : {_?ListQ ..}, vox : {_?NumberQ, _?NumberQ, _?NumberQ}] := ExportTracts["", tracts, vox, {0, 0, 0}, {}]

ExportTracts[tracts : {_?ListQ ..}, vox : {_?NumberQ, _?NumberQ, _?NumberQ}, dim : {_?NumberQ, _?NumberQ, _?NumberQ}] := ExportTracts["", tracts, vox, dim, {}]

ExportTracts[tracts : {_?ListQ ..}, vox : {_?NumberQ, _?NumberQ, _?NumberQ}, dim : {_?NumberQ, _?NumberQ, _?NumberQ}, seeds_?ListQ] := ExportTracts["", tracts, vox, dim, seeds]

ExportTracts[file_?StringQ, tracts : {_?ListQ ..}] := ExportTracts[file, tracts, {0, 0, 0}, {0, 0, 0}, {}]

ExportTracts[file_?StringQ, tracts : {_?ListQ ..}, vox : {_?NumberQ, _?NumberQ, _?NumberQ}] := ExportTracts[file, tracts, vox, {0, 0, 0}, {}]

ExportTracts[file_?StringQ, tracts : {_?ListQ ..}, vox : {_?NumberQ, _?NumberQ, _?NumberQ}, dim : {_?NumberQ, _?NumberQ, _?NumberQ}] := ExportTracts[file, tracts, vox, dim, {}]

ExportTracts[file_?StringQ, tracts : {_?ListQ ..}, vox : {_?NumberQ, _?NumberQ, _?NumberQ}, dim : {_?NumberQ, _?NumberQ, _?NumberQ}, seeds_?ListQ] := Block[{fileo},
	fileo = If[file == "",
		FileSelect["FileSave", {"*.trk"}, "trk", WindowTitle -> "Select the destination file"],
		ConvertExtension[file, ".trk"]
	];

	Export[fileo, {tracts, vox, dim, seeds}, {"trk", {"Tracts", "VoxelSize", "Dimensions", "Seeds"}}]
]


(* ::Subsubsection::Closed:: *)
(*ExportTractsDefault*)


ExportTractsDefault[file_, rule_, ___] := Block[{
		tracts, vox, dim, seeds, strm, nDim, nVox, nSeed, nTr, nTrCoor, nTrLeng
	},

	{tracts, vox, dim, seeds} = {"Tracts", "VoxelSize", "Dimensions", "Seeds"} /. rule;

	(*voxel size and dimensions*)
	nDim = If[dim =!= {0, 0, 0}, 3, 0];
	nVox = If[vox =!= {0, 0, 0}, 3, 0];
	(*Number of seed coordinates*)
	nSeed = 3 Length[seeds];
	(*number of tracts tracts coordinates and tract lenghts*)
	nTr = Length[tracts];
	nTrLeng = Length /@ tracts;
	nTrCoor = 3 Total[nTrLeng];

	(*open the stream*)
	strm = OpenWrite[file, BinaryFormat -> True];
	(*how to partition the stream*)
	BinaryWrite[strm, Flatten@{nDim, nVox, nSeed, nTrCoor, nTr, nTrLeng}, "Integer32"];
	(*Write the data*)
	BinaryWrite[strm, Flatten@{dim, vox, seeds, tracts}, "Real32"];
	(*close the stream*)
	Close[strm];
]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
