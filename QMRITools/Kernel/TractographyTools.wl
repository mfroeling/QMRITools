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
"FiberTractography[tensor, vox] performs fibertractography on the tensor with voxels dimensions vox.
FiberTractography[tensor, vox, {par, {min, max}}] performs fibertractography on the tensor with voxels dimensions vox with additional stoppin criteria
par, where tracts are only generated between values of par min and max.
FiberTractography[tensor, vox, {{par, {min, max}}, ..}] performs fibertractography on the tensor with voxels dimensions vox with 
multiple additional stopping criteria."


FitTracts::usage = 
"FitTracts[tract] fits a tract or a list of tracts, each defined as a list of {x, y, z} coordinates with a polinomial function.
FitTracts[tract, vox, dim] does the same but contrains all the tract coordinates to the volume difined by dim."

ResampleTracts::usage = 
"ResampleTracts[tracts, n] resample each Tract to exactly n vertices."

MoveTracts::usage = 
"MoveTracts[tracts, off] moves the tract coordicantes by off, which is {x, y, z}."

RescaleTracts::usage = 
"RescaleTracts[tracts, sc] scales the tract coordinates by 1/sc, which is {x, y, z} or single number."

FiberLength::usage = 
"FiberLength[tracts] calculates the length of each tract."

GetTractValues::usage = 
"GetTractValues[tracts, parameter, vox] gets the value of the parameter map at each tract coordinate."


FilterTracts::usage = 
"FilterTracts[tracts, vox, {select.. }] filters the tracts based on the list of select criteria.
Select criteria are defined as {\"logic\",{\"how\", criteria}}.
The \"logic\" parameter can be \"and\", \"or\" and \"not\".
The \"how\" parameter can be:
	- \"x\", \"y\", or \"z\" for slice selection, here criteria is a slice number
	- \"thourgh\" for selecting tract that go through a roi, here criteria is a 3D mask.
	- \"within\" for selecting tract that fit fully within the roi, here criteria is a 3D mask.
	- \"partwithin\" for selecting the part of the tracts that fall within the roi, here criteria is a 3D mask.
Any number of select criteria can be listed."

SegmentTracts::usage = 
"SegmentTracts[tracts, segs, vox, dim] segments the tracts based on segs."


SeedDensityMap::usage = 
"SeedDensityMap[seeds, vox, dim] makes a seed density map based on the seed loactions."

TractDensityMap::usage = 
"TractDensityMap[tracts, vox, dim] makes a tract density map based on the tracts vertices."

TractLengthMap::usage = 
"TractLengthMap[tracts, vox, dim] makes a tract length map based on the tracts lengths."

TractAngleMap::usage = 
"TractAngleMap[tracts, vox, dim] makes a tract angle map based on the tracts angles with the z-plane."


PlotTracts::usage = 
"PlotTracts[tracts, vox] plots the tracts assuming an Boxratio based on vox.
PlotTracts[tracts, vox, dim] plots the tracts assuming an Boxratio based on vox with a PlotRange spanning the full dim."

PlotSegmentedTracts::usage = 
"PlotSegmentedTracts[tracts, segments, dim, vox] plots the tracts after segmenting each segments.
PlotSegmentedTracts[tracts, segments, bones, dim, vox] plots the tracts after segmenting each segments also rendering a bone volume."


FindTensorPermutation::usage = 
"FindTensorPermutation[tensor, vox] performs tractography for all tensor permutations and gives back the one that has the longest tracts.
FindTensorPermutation[tensor, vox, {par, {min, max}}] same but with additional stoppin criteria par, where tracts are only generated between values of par min and max.
FindTensorPermutation[tensor, vox, {{par, {min, max}}, ..}] same but with with multiple additional stopping criteria.

Ouput = {permutations, flips, plot}

FindTensorPermutation[] is based on DOI: 10.1016/j.media.2014.05.012."


ImportTracts::usage = 
"ImportTracts[file] imports a *.trk file. It can contain {tracts, vox, dim, seeds}."

ExportTracts::usage = 
"ExportTracts[file, tracts, vox, dim, seeds] exports the tracts, vox, dim and seeds to *.trk file."


(* ::Subsection::Closed:: *)
(*Options*)


FittingOrder::usage = 
"FittingOrder is an option for FitTracts. It specifies the polinominal order of the function to fit the tract."

FitTractSegments::usage =
"FitTractSegments is an option for SegmentTracts. If set True the segmented tracts are fitted with FitTracts."

TracMonitor::usage = 
"TracMonitor is an option for FiberTractography. When set True it prints the progress."

FiberLengthRange::usage = 
"FiberLengthRange is an option for FiberTractography and specifies the allowed tract range."

FiberAngle::usage = 
"FiberAngle is an option for FiberTractography and specifies the allowed angle change per tract step."

TensorFilps::usage =
"TensorFilps is an option for FiberTractography and speciefies if the tensor orientation is fliped, see FlipTensorOrientation."

TensorPermutations::usage = 
"TensorPermutations is an option for FiberTractography and speciefies if the tensor orientation is permuted, see FlipTensorOrientation."

StopThreshhold::usage = 
"StopThreshhold is an option for FiberTractography and defines the stop threshhold which is a value between 0 and 1."

StepSize::usage = 
"StepSize is an option for FiberTractography and defines the tractography step size."

MaxSeedPoints::usage = 
"MaxSeedPoints is an option for FiberTractography and defines the maximum number of seedspoints to be used."


MaxTracts::usage = 
"MaxTracts is an option for PlotTracts. It specifies how many tracts are plotted."

TractSize::usage = 
"TractSize is an option for PlotTracts. When tubes are used it specifies the tube withd."

TractColoring::usage = 
"TractColoring is an option for PlotTracts and sets how the tracts are colored. Values can be \"Direction\", \"Length\", \"Angle\", {par}, or RGBColor[].
For \"Length\", \"Angle\", {par} it can be defined in the form {..., {min, max}} where the {min, max} specifies the range of the color function."

TractReduction::usage = 
"TractReduction is an option for PlotTracts. Value can be an Integer > 0, which determines with which facter the tract coordinates are subsampled."

TractScaling::usage = 
"TractScaling is an option for PlotTracts. The value can be \"World\" or \"Voxel\", if the value is \"Wold\" the tracts are in mm else in voxel coordinates."

NormalizeDensity::usage = 
"NormalizeDensity is an option for TractDensityMap. If set True the tractdensity is normalized, if False then it is the true tract count."



(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*Tract Operations*)


(* ::Subsubsection::Closed:: *)
(*FitTracts*)


Options[FitTracts] = {FittingOrder -> 3};

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

ResampleTracts[tracts_, len_] := Block[{int, r},
	r = Rescale[N@Range[len]];
	ToPackedArray@Map[(int = ListInterpolation[Transpose[{#}], {{0, 1}}, InterpolationOrder -> 1]; int /@ r) &, tracts]
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


RescaleTractsI = Compile[{{tr, _Real, 2}, {sc, _Real, 1}},
	Transpose[Transpose[tr]/sc]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


RescaleTractsC = Compile[{{tr, _Real, 2}, {sc, _Real, 1}},
	Transpose[Ceiling[Transpose[tr]/sc]]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


RescaleTractsCM = Compile[{{tr, _Real, 2}, {vox, _Real, 1}},
	Transpose[Ceiling[Transpose[Mean[{tr[[2 ;; -1]], tr[[1 ;; -2]]}]]/vox]], 
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


(* ::Subsubsection::Closed:: *)
(*FiberLength*)


SyntaxInformation[FiberLength] = {"ArgumentsPattern" -> {_}};

FiberLength[tracts_]:=FLengthC[tracts]

FLengthC = Compile[{{trc,_Real,2}},
	Total[Norm/@Differences[trc]]
,RuntimeAttributes->{Listable},RuntimeOptions->"Speed"]


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
(*GetTractValues*)


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

TractLengthMap[tracts_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, dim_] := GatherThread[RescaleTractsC[tracts, vox], FiberLength[tracts], dim]


(* ::Subsubsection::Closed:: *)
(*TractAngleMap*)


SyntaxInformation[TractAngleMap] = {"ArgumentsPattern" -> {_, _, _}};

TractAngleMap[tracts_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, dim_] := GatherThread[RescaleTractsCM[tracts, vox], VecAngC[tracts], dim]


(* ::Subsubsection::Closed:: *)
(*VecAngC*)


VecAngC = Compile[{{tr, _Real, 2}},
	Abs[(180./Pi) ArcCos[(#[[1]]/Sqrt[Total[Abs[#]^2]]) & /@ (tr[[2 ;; -1]] - tr[[1 ;; -2]])] - 90]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*GatherThread*)


GatherThread[coor_, val_, dim_] := Block[{len, ran, list, out},
	(*Split in sublists to prevent memory crash for large tracts sets*)
	len = Length[coor];
	ran = Range[0, len, 10000];
	list = {1, 0} + # & /@ Partition[If[Last@ran =!= len, Append[ran, len], ran], 2, 1];
	(*Median of Medians is a good approximation*)
	out = Table[GatherThread[coor[[i[[1]] ;; i[[2]]]], val[[i[[1]] ;; i[[2]]]]], {i, list}];
	out = Select[GatherThread[out], (1 <= #[[1, 1]] <= dim[[1]] && 1 <= #[[1, 2]] <= dim[[2]] && 1 <= #[[1, 3]] <= dim[[3]]) &];
	ToPackedArray@N@Normal@SparseArray[out, dim]
]

GatherThread[coor_, val_] := GatherThread[Thread /@ Thread[coor -> val]]

GatherThread[rule_] := Block[{out},
	out = GatherBy[Flatten[rule, 1], First];
	Thread[out[[All, 1, 1]] -> (Median /@ out[[All, All, 2]])]
]


(* ::Subsection:: *)
(*FiberTractography*)


(* ::Subsubsection::Closed:: *)
(*FiberTractography*)


Options[FiberTractography] = {
	FiberLengthRange -> {20, 500},
	FiberAngle -> 30,
	TensorFilps -> {1, 1, 1},
	TensorPermutations -> {"x", "y", "z"},
	InterpolationOrder -> 1,
	StopThreshhold -> 0.5,
	StepSize -> Automatic,
	Method -> "Euler",
	MaxSeedPoints -> Automatic,
	TracMonitor->True
};

SyntaxInformation[FiberTractography] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

FiberTractography[tensor_, voxi_, opts : OptionsPattern[]] := FiberTractography[tensor, voxi, {0. First@tensor, {0., 1.}}, opts]

FiberTractography[tensor_, voxi_, {par_?ArrayQ, {min_?NumberQ, max_?NumberQ}}, opts : OptionsPattern[]] := FiberTractography[tensor, voxi, {{par, {min, max}}}, opts]

FiberTractography[tensor_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, inp : {{_, {_, _}} ...}, OptionsPattern[]] := Block[{
		lmin, lmax, amax, maxSeed, flip, per, int, stopT, step, tracF, trFunc,
		tens, tensMask, inpTr, treshs, stop, coors, intE, intS, ones, dim, crp,
		seedN, seedI, seedT, seeds, t1, tracts, iii, drop, smax, len ,sel
	},

	(*get the options*)
	{lmin, lmax} = OptionValue[FiberLengthRange];
	amax = OptionValue[FiberAngle];
	maxSeed = OptionValue[MaxSeedPoints];
	
	flip = OptionValue[TensorFilps];
	per = OptionValue[TensorPermutations];
	
	int = OptionValue[InterpolationOrder];
	stopT = OptionValue[StopThreshhold];
	
	step = OptionValue[StepSize];
	step = N@If[NumberQ[step],step, Min[0.5 vox]];
	
	tracF = Switch[OptionValue[Method], "RungeKutta" | "RK" | "RK2", RK2, "RungeKutta4" | "RK4", RK4, _, Euler];
	smax = Ceiling[(lmax/step)];

	(*prepare tensor and stop data, remove background for performance and make int functions*)
	dim = Rest@Dimensions@tensor;
	{stop, crp} = AutoCropData[Unitize[First@tensor] Mask[inp], CropPadding -> 0];
	tens = ApplyCrop[#,crp]&/@FlipTensorOrientation[tensor, per, flip];

	intE = MakeIntFunction[RotateDimensionsLeft[{tens}, 2], vox, int];
	intS = MakeIntFunction[N@Normal@stop, vox, int];

	(*make the random seed points*)
	SeedRandom[1234];
	seedN = Total@Flatten@stop;

	maxSeed = Round@Which[
		NumberQ[maxSeed], maxSeed, 
		maxSeed === Automatic, seedN,
		Head[maxSeed]===Scaled, First[maxSeed] seedN
	];

	seedI = SparseArray[stop]["NonzeroPositions"];
	seeds = {};

	While[Length[seeds] < maxSeed,
		seedT = Flatten[RandomSample[seedI, #] & /@ Block[{i = maxSeed},
			First@Last@Reap[While[i > 0, Sow[If[i > seedN, seedN, i]];
			i = i - seedN;]]], 1];
		seedT = # vox & /@ (seedT + RandomReal[{-0.99, 0}, {maxSeed, 3}]);
		seedT = Pick[seedT, intS @@@ seedT, 1.];
		
		seeds = Join[seeds, seedT];
		seeds = seeds[[;; Min[{Length@seeds, maxSeed}]]]
	];
	seedN = Length@seeds;

	If[OptionValue[TracMonitor],
		Echo["Starting tractography for "<>ToString[seedN]<>" seed points with stepsize "<>ToString[step]<>" mm"];
	];

	(*perform tractography for each seed point*)	
	t1 = First[AbsoluteTiming[
		trFunc = TractFunc[#, step, {amax, smax, stopT}, {intE, intS, tracF}]&;
		tracts = If[OptionValue[TracMonitor],
			Monitor[Table[trFunc@seeds[[iii]], {iii, 1, seedN}], ProgressIndicator[iii, {0, seedN}]],
			Table[trFunc@seeds[[iii]], {iii, 1, seedN}]
		];
	]];

	If[OptionValue[TracMonitor],
		Echo["Checking Lengths"];
	];

	(*select only tracts within correct range and clip tracts that are longer*)
	{tracts, sel} = FilterTractLength[tracts, {lmin, lmax}, "both"];
	seeds = Pick[seeds, sel, 1];
	
	(*report timing*)
	If[OptionValue[TracMonitor],
		Echo["Tractography took "<>ToString[Round[t1,.1]]<>" seconds ("<>ToString[Round[seedN/t1]]<>" tracts/s)"];
		Echo[ToString[Length[tracts]]<>" valid tracts with length "<>ToString[Round[step Mean[Length /@ tracts], 0.1]]<>"\[PlusMinus]"<>ToString[Round[step StandardDeviation[Length /@ tracts], 0.1]]<>" mm"];	
	];

	(*output tracts*)
	MoveTracts[#, vox (crp[[;; ;; 2]] - 1)]&/@{tracts, seeds}
]


(* ::Subsubsection::Closed:: *)
(*TractFunc*)


TractFunc[loc0_?VectorQ, h_, stp_, fun_] := Block[{dir0},
	(*startpoint*)
	dir0 = h EigVec[First[fun]@@loc0, {0,0,1}];
	(*tracts bi-directional*)
	ToPackedArray@Join[
		Reverse@TractFuncI[{loc0 + dir0/2, dir0}, h, stp, fun], 
		TractFuncI[{loc0 - dir0/2, -dir0}, h, stp, fun]
	]
]


TractFuncI[{loc0_?VectorQ, step0_?VectorQ}, h_, {amax_, smax_, stoptr_}, {intE_, intS_, tracF_}] := Block[{loc, locn, step, stepn},
	(*perform tractography*)
	NestWhileList[(
			(*get new location and the step at new location*)
			{loc, step} = #[[1;;2]];
			locn = loc + step;
			stepn = tracF[locn, step, h, intE];
			
			(*output, new location and its direction, angle and stop*)
			{locn, stepn, VecAng[step, stepn], intS@@(locn+stepn)}
		)&, {loc0, step0, 0, intS@@(loc0+step0)}, 
	(#[[3]] < amax && #[[4]] > stoptr)&, 1, smax][[All,1]]
];


(* ::Subsubsection::Closed:: *)
(*VecAng*)


VecAng = Compile[{{v1, _Real, 1}, {v2, _Real, 1}}, Block[{v, n1 ,n2},
	n1 = Sqrt[v1 . v1];
	n2 = Sqrt[v2 . v2];
	v = If[n1 > 0. && n2 > 0., (v1 . v2)/(n1 n2), v1 . v2];
	Re[180. ArcCos[v] / Pi]
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*EigVec*)


EigVec = Compile[{{t, _Real, 1}, {vdir, _Real, 1}}, Block[{
		dxx, dyy, dzz, dxy, dxz, dyz, dxy2, dxz2, dyz2, tens,
		i1, i2, i3, i, v, s, v2, l1, a, b, c, norm, vec
	},
	tens = 1000 t;
	(*method https://doi.org/10.1016/j.mri.2009.10.001*)
	vec = If[Total[Abs[tens]]<10.^-15,
		{0., 0., 0.},

		{dxx, dyy, dzz, dxy, dxz, dyz} = tens;
		{dxy2, dxz2, dyz2} = {dxy, dxz, dyz}^2;

		i1 = dxx + dyy + dzz;
		i2 = dxx dyy + dxx dzz + dyy dzz - dxy2 - dxz2 - dyz2;
		i3 = dxx dyy dzz + 2 dxy dxz dyz - dzz dxy2 - dyy dxz2 - dxx dyz2;
		
		i = i1/3;
		v = i^2 - i2/3;
		s = i^3 - (i1 i2)/6 + i3/2;

		v2 = Sqrt[v];
		l1 = i + 2 v2 Cos[Re[ArcCos[If[(v v2)===0, 0., s/(v v2)]]/3]];
		
		{a, b, c} = {dxz dxy, dxy dyz, dxz dyz} - {dyz, dxz, dxy} ({dxx, dyy, dzz} - l1);
		{a, b, c} = {b c, a c, a b};
		norm = Sqrt[a^2 + b^2 + c^2];
		
		If[norm < 10.^-30, {0., 0., 0.}, {a, b, c} / norm]
	];

	 Sign[Sign[vdir . vec] + 0.1] vec
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsubsection::Closed:: *)
(*Euler*)


Euler[y_, v_, h_, func_] := Block[{k1},
	k1 = h EigVec[func@@(y), v];
	k1
]


(* ::Subsubsection::Closed:: *)
(*RK2*)


RK2[y_, v_, h_, func_] := Block[{k1, k2},
	k1 = h EigVec[func@@(y), v];
	k2 = h EigVec[func@@(y + k1/2), v];
	k2	
]


(* ::Subsubsection::Closed:: *)
(*RK4*)


RK4[y_, v_, h_, func_] := Block[{k1, k2, k3, k4},
	k1 = h EigVec[func@@(y), v];
	k2 = h EigVec[func@@(y + k1/2), v];
	k3 = h EigVec[func@@(y + k2/2), v];
	k4 = h EigVec[func@@(y + k3), v];
	k1/6 + k2/3 + k3/3 + k4/6
]


(* ::Subsection::Closed:: *)
(*FindTensorPermutation*)


Options[FindTensorPermutation] = {
	FiberLengthRange -> {20, 500},
	FiberAngle -> 30,
	InterpolationOrder -> 0,
	StopThreshhold -> 0.5,
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
		tracts = First@FiberTractography[tens, vox, stop, TracMonitor -> False,TensorFilps -> flips[[i]], TensorPermutations -> perms[[j]],
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
	Method->"line", 
	TractColoring->"Direction",
	ColorFunction->"SouthwestColors",
	Boxed->True,
	TractSize-> 1,
	TractReduction->1,
	TractScaling -> "World",
	PerformanceGoal->"Quality"
}

SyntaxInformation[PlotTracts] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

PlotTracts[tracts_, voxi_, opts : OptionsPattern[]] := PlotTracts[tracts, voxi, 0, opts]

PlotTracts[tracts_, voxi_, dimi_, OptionsPattern[]] := Block[{
		range, vox, size, select, opts, col, tube, line, plot, colOpt, 
		met, pran, max, colf, colArray, sc, n
	},
	
	(*reduce points along tracts*)
	(*tracts = Select[tractsI[[All,1;;-1;;OptionValue[TractReduction]]],Length[#]>3&];*)

	(*select correct number of tracts*)
	max = OptionValue[MaxTracts];
	(*If[OptionValue[Method]==="tube", max =  Min[4000, max]];*)
	SeedRandom[1234];
	select = ToPackedArray/@RandomSample[tracts, Min[max, Length[tracts]]];
	n = Round[Quantile[Length/@select, 0.75] / OptionValue[TractReduction]];
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
		StringQ[colOpt]||ArrayQ[colOpt], {colOpt, Automatic}, 
		Length[colOpt]===2, colOpt, 
		MemberQ[{Hue, RGBColor, GrayLevel, CMYKColor, LABColor, LCHColor, LUVColor, XYZColor},Head[colOpt]], {"Color", RGBColor@colOpt},
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
		_,
		White
	];
		
	(*make the plot*)
	opts = Sequence[{
		Method -> {"TubePoints" -> {If[OptionValue[PerformanceGoal]==="Quality",5,3], 2}}, 
		Lighting -> "Accent", 
		ImageSize -> OptionValue[ImageSize], 
		SphericalRegion -> True, Boxed -> OptionValue[Boxed],
		Background -> Lighter@Gray, BoxRatios -> size, PlotRange -> range, 
		Axes -> OptionValue[Boxed], LabelStyle -> Directive[{Bold, 16, White}]
	}];

	
	select = Reverse[select, 3];
	plot = Graphics3D[Switch[OptionValue[Method],
		"tube", 
		{CapForm["Butt"], JoinForm["Miter"], 
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
	Monitor -> False
}

SyntaxInformation[PlotSegmentedTracts] = {"ArgumentsPattern" -> {_, _, _, _, _., OptionsPattern[]}};

PlotSegmentedTracts[tracts_, segments_, dim_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, opts : OptionsPattern[]] := PlotSegmentedTracts[tracts, segments, None, dim, vox, opts]

PlotSegmentedTracts[tracts_, segments_, bones_, dim_, vox:{_?NumberQ,_?NumberQ,_?NumberQ}, opts : OptionsPattern[]] := Block[{
		ntr, fran, type, segs, tractsF, tractsFI, ran, rand, colList, showF,
		ref, bon, musc, trac, tracksSel, lengs, nTracts, sel, mon
	},
	
	(*get options*)
	ntr = OptionValue[MaxTracts];
	fran = OptionValue[FiberLengthRange];
	type = OptionValue[Method];
	mon = OptionValue[Monitor];

	(*prepare data*)
	segs = Transpose@segments;
	If[mon, Echo["Fitting tracts"]];
	SeedRandom[1234];
	tractsF = RandomSample[tracts, Min[{5 ntr, Length@tracts}]];
	tractsF = FitTracts[tractsF, vox, dim, FittingOrder -> 3];
	tractsFI = RescaleTractsC[tractsF, vox];

	(*make colors*)
	SeedRandom[1234];
	ran = Range[Length@segs];
	colList = Reverse[ColorData["DarkRainbow"] /@ Rescale[ran]][[RandomSample[ran]]];

	(*rand = RandomSample[ran];
	colList = ColorData["DarkRainbow"] /@ N[Rescale[rand]];*)

	(*reference environement*)
	If[mon, Echo["Making muscle iso volumes"]];
	ref = PlotTracts[tractsF, vox, dim, MaxTracts -> 1, Method -> "line", TractColoring -> RGBColor[0, 0, 0, 0]];

	(*make the muscle contours*)
	musc = Table[PlotContour[segs[[i]], vox, ContourOpacity -> 0.3, ContourColor -> colList[[i]], ContourSmoothRadius -> 2, ContourResolution -> {12, 3, 3}], {i, ran}];
	bon = If[bones =!= None, PlotContour[bones, vox, ContourOpacity -> 1, ContourColor -> Gray, ContourSmoothRadius -> 2, ContourResolution -> {12, 3, 3}], Graphics3D[]];

	If[mon, Echo["Making per muscle tracts"]];
	(*select the tracts per muscle and make fiber plots*)
	tracksSel = FilterTracts[tractsF, tractsFI, {{"and", {"partwithin", #}}}, FiberLengthRange -> fran] & /@ segs;
	lengs = Length /@ tracksSel;
	nTracts = Round[ntr lengs/Total[lengs]];
	sel = UnitStep[nTracts - 11];
	trac = MapThread[If[#4 =!= 1, Graphics3D[],
		PlotTracts[#1, vox, dim, MaxTracts -> #2, Method -> type, TractSize -> 1, TractColoring -> #3, TractReduction -> 5, PerformanceGoal -> "Speed"]
	] &, {tracksSel, nTracts, colList, sel}];

	If[mon, Echo["Finalizing scenes"]];
	showF = Show[ref, ##, ImageSize -> OptionValue[ImageSize], Axes -> False, Boxed -> False, ViewPoint -> {0., -2., 1.}] & @@ # &;
	
	Switch[OptionValue[OutputForm],
		"All", showF[{bon, trac, musc}],
		"Groups", showF[{#}] & /@ {bon, trac, musc},
		"Joined", showF[{bon, #}] & /@ Thread[{trac, musc}],
		"Individual", {showF[{bon}], showF[{#}] & /@ trac, showF[{#}] & /@ musc}
	]
]


(* ::Subsection::Closed:: *)
(*ImportTracts*)


(* ::Subsubsection::Closed:: *)
(*RegisterImport*)


ImportExport`RegisterImport["trk", ImportTractsDefault, {}, 
  "AvailableElements" -> {"Tracts", "VoxelSize", "Dimensions", 
    "Seeds"}, "OriginalChannel" -> True];


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
	(*number of tracts tracts coordinates and trac lenghts*)
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
