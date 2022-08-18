(* ::Package:: *)

(* ::Title:: *)
(*QMRITools VisteTools*)


 


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


FitTract::usage = "
FitTract[tract] fits a tract defined as a list of {x,y,z} coordinates with a polinomial function.
FitTract[{tract, ...}] fits a list of tracts defined as a list of {x,y,z} coordinates with a polinomial function."

FiberTractography::usage =
"FiberTractography[tensor, vox] performs fibertractography on the tensor with voxels dimensions vox.
FiberTractography[tensor, vox, {par, {min, max}}] performs fibertractography on the tensor with voxels dimensions vox with additional stoppin criteria
par, where tracts are only generated between values of par min and max.
FiberTractography[tensor, vox, {{par, {min, max}}, ..}] performs fibertractography on the tensor with voxels dimensions vox with 
multiple additional stopping criteria."

FindTensorPermutation::usage = 
"FindTensorPermutation[tensor, vox] performs tractography for all tensor permutations and gives back the one that has the longest tracts.
FindTensorPermutation[tensor, vox, {par, {min, max}}] same but with additional stoppin criteria par, where tracts are only generated between values of par min and max.
FindTensorPermutation[tensor, vox, {{par, {min, max}}, ..}] same but with with multiple additional stopping criteria.

Ouput = {permutations, flips, plot}

FindTensorPermutation[] is based on DOI: 10.1016/j.media.2014.05.012."


SelectTractTroughPlane::usage = ""
SelectTractTroughVol::usage = ""
SelectTractInVol::usage = ""
SelectTractPartInVol::usage = "" 
PartTracts::usage = ""
SelectTracts::usage = ""
CombineROIs::usage = ""
TractDensityMap::usage = "TractDensityMap[tracts_, vox_, dim_]"
FilterTracts::usage = "FilterTracts[tracts_, vox_, select_]"
MakeColor::usage = ""

(* ::Subsection::Closed:: *)
(*Options*)


FittingOrder::usage = 
"FittingOrder is an option for FitTract. It specifies the polinominal order of the function to fit the tract."

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


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*FiberTractography*)


(* ::Subsubsection::Closed:: *)
(*FiberTractography*)


Options[FiberTractography] = {
	FiberLengthRange -> {10, 200},
	FiberAngle -> 30,
	TensorFilps -> {1, 1, 1},
	TensorPermutations -> {"x", "y", "z"},
	InterpolationOrder -> 1,
	StopThreshhold -> 0.5,
	StepSize -> 2,
	Method -> "Euler",
	MaxSeedPoints -> Infinity,
	TracMonitor->True
};

SyntaxInformation[FiberTractography] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

FiberTractography[tensor_, voxi_, opts : OptionsPattern[]] := FiberTractography[tensor, voxi, {0. First@tensor, {0., 1.}}, opts]

FiberTractography[tensor_, voxi_, {par_?ArrayQ, {min_?NumberQ, max_?NumberQ}}, opts : OptionsPattern[]] := FiberTractography[tensor, voxi, {{par, {min, max}}}, opts]

FiberTractography[tensor_, vox_, inp : {{_, {_, _}} ...}, opts : OptionsPattern[]] := Block[{
		lmin, lmax, amax, maxSeed, flip, per, int, stopT, step, tracF, trFunc,
		tens, tensMask, inpTr, treshs, stop, coors, intE, intS,
		ones, n, seeds, t1, tracts, iii, drop, smax		 
	},
	
	(*get the options*)
	{lmin, lmax} = OptionValue[FiberLengthRange];
	amax = OptionValue[FiberAngle];
	maxSeed = OptionValue[MaxSeedPoints];
	
	flip = OptionValue[TensorFilps];
	per = OptionValue[TensorPermutations];
	
	int = OptionValue[InterpolationOrder];
	stopT = OptionValue[StopThreshhold];
	
	step = N@OptionValue[StepSize];
	tracF = Switch[OptionValue[Method], "RungeKutta" | "RK" | "RK2", RK2, "RungeKutta4" | "RK4", RK4, _, Euler];
 	smax = Ceiling[(lmax/step)];
 	
 	
	(*prepare tensor and stop data*)
	tens = FlipTensorOrientation[tensor, per /. Thread[{"x","y","z"} -> {"z","y","x"}], Reverse[flip]];
	stop = Unitize[Total@Abs@tens] Mask[inp];
			
	(*make the trac and stop interpolation functions*)
	intE = ToIntFunction[RotateDimensionsLeft[{tens}, 2], vox, int];
	intS = ToIntFunction[stop, vox, int];
	
	(*make the random seed points*)
	SeedRandom[1234];
	seeds = Threaded[vox] (RandomChoice[SparseArray[stop]["NonzeroPositions"], Round[1.2 maxSeed]] + RandomReal[{-0.5, 0.5}, {Round[1.2 maxSeed], 3}]);
    seeds = Pick[seeds, intS @@@ seeds, 1.];
    seeds = seeds[[;;Min[{Length@seeds,maxSeed}]]];

	(*perform tractography for each seed point*)	
	t1 = First[AbsoluteTiming[
		trFunc = TractFunc[#, step, {amax, smax, stopT}, {intE, intS, tracF}]&;
		tracts = Monitor[
			Table[trFunc@seeds[[iii]], {iii, 1, Length[seeds]}]
			, ProgressIndicator[iii, {0, Length[seeds]}]];
	]];
	
	len = Length /@ tracts;
	sel = UnitStep[len - lmin/step];
	tracts = Pick[tracts, sel, 1];
	seeds = Pick[seeds, sel, 1];
	
	(*select only tracts within correct range and clip tracts that are longer*)
	tracts = If[step Length[#] >= lmax, 
		drop = Ceiling[(step Length[#] - lmax)/(2 step)] + 1; #[[drop ;; -drop]], 
		#] & /@ tracts;
	
	(*report timing*)
	If[OptionValue[TracMonitor],
		Print["number of valid seedvoxels: ", Total@Flatten@stop];
		Print["number of selected seedpoints: ", Length[seeds]];
		Print["Tractography took ", t1, " seconds"];
		Print["number valid tracts: ", Length[tracts]];	
		Print[Length[tracts], " tracts within length range were selecte with length ", Round[step Mean[Length /@ tracts], 0.1],"\[PlusMinus]", Round[step StandardDeviation[Length /@ tracts], 0.1] , " mm"];
	];

	(*output tracts*)
	{tracts, seeds}
]


(* ::Subsubsection::Closed:: *)
(*ToIntFunction*)


ToIntFunction[dat_, vox_, int_] := Block[{def, range},
	range = Thread[{vox, vox Dimensions[dat][[1;;3]]}]-0.5 vox;
	def = 0. dat[[1,1,1]];
	def =If[ListQ[def], Flatten@def, def];
	With[{ex=def},ListInterpolation[dat, range, InterpolationOrder -> int, 
		"ExtrapolationHandler" -> {ex &, "WarningMessage" -> False}, Method -> "Hermite"]]
]


(* ::Subsubsection::Closed:: *)
(*Euler*)


Euler[y_, v_, h_, func_] := Block[{k1},
	k1 = VecAlign[v, h func[y]];
	k1
]


(* ::Subsubsection::Closed:: *)
(*RK2*)


RK2[y_, v_, h_, func_] := Block[{k1, k2},
	k1 = VecAlign[v, h func[y]];
	k2 = VecAlign[v, h*func[y + k1/2]];
	k2	
]


(* ::Subsubsection::Closed:: *)
(*RK4*)


RK4[y_, v_, h_, func_] := Block[{k1, k2, k3, k4},
	k1 = VecAlign[v, h func[y]];
	k2 = VecAlign[v, h func[y + k1/2]];
	k3 = VecAlign[v, h func[y + k2/2]];
	k4 = VecAlign[v, h func[y + k3]];
	k1/6 + k2/3 + k3/3 + k4/6
]


(* ::Subsubsection::Closed:: *)
(*VecAlign*)
  
  
VecAlign = Compile[{{v1, _Real, 1}, {v2, _Real, 1}}, Sign[Sign[v1 . v2] + 0.1] v2,
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsubsection::Closed:: *)
(*VecAng*)


VecAng = Compile[{{v1, _Real, 1}, {v2, _Real, 1}}, Block[{v,n},
	n = Norm[v1] Norm[v2];
	If[n > 0.,
		v = (v1 . v2) / n;
		If[-1. < v < 1., 180 ArcCos[v] / Pi, 0.],
		0.
	]
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsubsection::Closed:: *)
(*EigV*)


(*fast eigenvector calculation for tractography*)
EigV[{0., 0., 0., 0., 0., 0.}] := {0., 0., 1.}
EigV[{xx_, yy_, zz_, xy_, xz_, yz_}] := First@Eigenvectors@{{xx, xy, xz}, {xy, yy, yz}, {xz, yz, zz}};


(* ::Subsubsection::Closed:: *)
(*TractFunc*)


TractFunc[loc0_, h_, {amax_, smax_, stoptr_}, {intE_, intS_, tracF_}] := Block[{st0, dir0, step1, in1, step2, in2},
	(*initial direction*)
	dir0 = Quiet[EigV[intE@@loc0]];
	
	(*place two seed points arount initial points*)
	step1 = tracF[loc0, dir0, h, Quiet[EigV[intE@@#]]&];
	in1 = {loc0 - 0.5 step1, step1};
	step2 = tracF[loc0, -dir0, h, Quiet[EigV[intE@@#]]&];
	in2 = {loc0 - 0.5 step2, step2};
	
	(*tractography in the two directions*)
	ToPackedArray@Join[
		Reverse@TractFunci[in1, h, {amax, smax, stoptr}, {intE, intS, tracF}],
		TractFunci[in2, h, {amax, smax, stoptr}, {intE, intS, tracF}]
	]
]


(* ::Subsubsection::Closed:: *)
(*TractFunci*)


TractFunci[{loc0_, step0_}, h_, {amax_, smax_, stoptr_}, {intE_, intS_, tracF_}] := Block[
	{ang, stop, loc, locn, step, stepn, angt,out},
	
	ang = 0;
	stop = 1;
	
	If[(intS@@(loc0 + step0)) < stoptr || (intS@@loc0) < stoptr,
		{},
		out = NestWhileList[(
			{loc,step} = #;
						
			(*get new location and the step at new location*)
			locn = loc + step;
			stepn = tracF[locn, step, h, Quiet[EigV[intE@@#]]&];
			
			(*get stop criteria*)
			ang = VecAng[step, stepn];
			stop = intS@@(locn + stepn);
						
			(*output*)
			{locn, stepn})&,
			
			{loc0, step0}, (ang < amax && stop > stoptr)&, All, smax
		];
		If[Length[out] > 1, out[[2 ;;, 1]]]
	]
];


(* ::Subsection::Closed:: *)
(*FitTract*)


Options[FitTract] = {FittingOrder -> 4};

SyntaxInformation[FitTract] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

FitTract[tract_, OptionsPattern[]] := ToPackedArray/@FitTractC[tract, Round@OptionValue[FittingOrder]]

FitTractC = Compile[{{trf, _Real, 2}, {ord, _Real, 0}}, Block[{mat = #^Range[0., ord] & /@ Range[Length[trf]]},
	mat . (PseudoInverse[mat] . trf)
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsection:: *)
(*ROIs*)


(* ::Subsubsection::Closed:: *)
(*SelectTractTroughPlane*)


SelectTractTroughPlane[tracts_, {dir_, slice_}, vox_]:=	SelectTractTroughPlaneC@@Switch[ToLowerCase[dir],
		"x", {tracts[[All,All,3]], slice vox[[3]] - 0.5},
		"y", {tracts[[All,All,2]], slice vox[[2]] - 0.5},
		"z", {tracts[[All,All,1]], slice vox[[1]] - 0.5}
	];

SelectTractTroughPlaneC = Compile[{{tract, _Real, 1}, {val, _Real, 0}},
	Boole[0 < Total[UnitStep[tract - val]] < Length[tract]], 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*SelectTractTroughVol*)


SelectTractTroughVol[tracts_, roi_, vox_]:= SelectTractTroughVolC[Normal@roi, tracts, vox]

SelectTractTroughVolC = Compile[{{roi, _Integer, 3}, {tract, _Real, 2}, {vox, _Real, 1}},
	Max[roi[[#[[1]], #[[2]], #[[3]]]] &[Ceiling[#/vox]] & /@ tract], 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*SelectTractInVol*)


SelectTractInVol[tracts_, roi_, vox_]:= SelectTractInVolC[Normal@roi, tracts, vox]

SelectTractInVolC = Compile[{{roi, _Integer, 3}, {tract, _Real, 2}, {vox, _Real, 1}},
	Floor[Total[roi[[#[[1]], #[[2]], #[[3]]]] &[Ceiling[#/vox]] & /@ tract]/Length[tract]],
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*SelectTractPartInVol*)


SelectTractPartInVol[tracts_, roi_, vox_]:= SelectTractPartInVolC[Normal@roi, tracts, vox]

SelectTractPartInVolC = Compile[{{roi, _Integer, 3}, {tract, _Real, 2}, {vox, _Real, 1}},
	(roi[[#[[1]], #[[2]], #[[3]]]] &[Ceiling[#/vox]]) & /@ tract, 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*FilterTracts*)


FilterTracts[tracts_, vox_, select_] := SelectTracts[tracts,
	CombineROIs[{#[[1]], Switch[ToLowerCase[#[[2, 1]]],
		"through", SelectTractTroughVol,
		"within", SelectTractInVol,
		"partwithin", SelectTractPartInVol,
		"x" | "y" | "z", SelectTractTroughPlane
	][tracts, If[StringLength[#[[2, 1]]] == 1, #[[2]], #[[2, 2]]], vox]} & /@ select]
]


(* ::Subsubsection::Closed:: *)
(*SelectTracts*)


SelectTracts[tracts_, sel_] := ToPackedArray/@PartTracts[Pick[tracts, Clip[Round@sel, {0,1}], 1]]


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


(* ::Subsection::Closed:: *)
(*TractDensityMap*)


TractDensityMap[tracts_, vox_, dim_] := Normal@SparseArray[Normal@Counts@ThreadedDiv[Flatten[tracts, 1], vox], dim];

ThreadedDiv = Compile[{{coor, _Real, 1}, {vox, _Real, 1}}, 
	Ceiling[coor/vox],
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsection:: *)
(*plotTracts*)


MakeColor[tract : {{_?NumberQ, _?NumberQ, _?NumberQ} ..}] := 
 Block[{dirs},
  dirs = Abs[Differences[tract]];
  ToPackedArray[
   Normalize[#] & /@ 
    Mean[{Prepend[dirs, dirs[[1]]], Append[dirs, dirs[[-1]]]}]]
  ]
MakeColor[tracts : {_?ListQ ..}] := MakeColor /@ tracts

MakeColor[tract_, fun_, ran_] := Block[{prange, colFunc, cval, fore},
  prange = ran;
  colFunc = ColorData["ThermometerColors"];
  cval = fun @@ # & /@ tract;
  fore = Unitize[cval];
  (fore (colFunc /@ Rescale[cval, prange])) + Gray (1 - fore)
  ]
  

(* ::Subsection::Closed:: *)
(*FindTensorPermutation*)


Options[FindTensorPermutation] = {
	FiberLengthRange -> {10, 200},
	FiberAngle -> 30,
	InterpolationOrder -> 0,
	StopThreshhold -> 0.5,
	StepSize -> 2,
	Method -> "Euler",
	MaxSeedPoints -> 500
};

SyntaxInformation[FindTensorPermutation] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

FindTensorPermutation[tensor_, vox_, opts : OptionsPattern[]] := FindTensorPermutation[tensor, vox, {1, {0, 1}}, opts]

FindTensorPermutation[tensor_, vox_, {par_?ArrayQ, {min_?NumberQ, max_?NumberQ}}, opts : OptionsPattern[]] := FindTensorPermutation[tensor, vox, {{par, {min, max}}}, opts]

FindTensorPermutation[tens_, vox_, stop_, OptionsPattern[]] := Block[{
	lmin, lmax, amax, maxSeed, int, stopT, step, intF, tracts, perms, flips, lengs, i, j, pl, ind
	},
	
	{lmin, lmax} = OptionValue[FiberLengthRange];
	amax = OptionValue[FiberAngle];
	maxSeed = OptionValue[MaxSeedPoints];
	int = OptionValue[InterpolationOrder];
	stopT = OptionValue[StopThreshhold];
	step = OptionValue[StepSize];
	intF = OptionValue[Method];
	
	perms = {{"x", "y", "z"}, {"x", "z", "y"}, {"y", "x", "z"}, {"y","z", "x"}, {"z", "x", "y"}, {"z", "y", "x"}};
	flips = {{1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1}};
	ind = {1, 1};
	PrintTemporary[Dynamic[{ind, flips[[ind[[1]]]], perms[[ind[[2]]]]}]];
	
	lengs = Table[
		ind = {i, j};
		tracts = FiberTractography[tens, vox, stop,
			FiberLengthRange -> {lmin, lmax}, FiberAngle -> amax, StopThreshhold -> stopT,
			InterpolationOrder -> int, StepSize -> step, Method -> intF, MaxSeedPoints -> maxSeed,
			TensorFilps -> flips[[i]], TensorPermutations -> perms[[j]], TracMonitor -> False
		];
		N@Mean[Length /@ tracts] step
	, {i, 1, 4}, {j, 1, 6}];
	
	pl = ArrayPlot[lengs, ColorFunction -> "Rainbow", ImageSize -> 150, Frame -> True, 
		FrameTicks -> {{Thread[{Range[4], flips}], None},{None, Thread[{Range[6], Rotate[#, 90 Degree] & /@ perms}]}}];
	
	{i, j} = FirstPosition[lengs, Max[lengs]];
	{perms[[j]], flips[[i]], pl}
]


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
