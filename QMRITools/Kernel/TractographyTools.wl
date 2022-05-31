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


(* ::Subsection::Closed:: *)
(*Error Messages*)



(* ::Section:: *)
(*Functions*)


Begin["`Private`"] 


(* ::Subsection:: *)
(*FitTract*)


Options[FitTract] = {FittingOrder -> 4};

SyntaxInformation[FitTract] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

FitTract[tract_, OptionsPattern[]] := Block[{x, y, z, xx, rule, ord, tr},
	ord = xx^# & /@ Range[0, OptionValue[FittingOrder]];
	If[ArrayDepth[tract] === 2,
		rule = xx -> Range[Length[tract]];
		Transpose[(Fit[#, ord, xx] /. rule) & /@ Transpose[tract]]
		,
		(
			tr = #;
			rule = xx -> Range[Length[tr]];
			Transpose[(Fit[#, ord, xx] /. rule) & /@ Transpose[tr]]
		) & /@ tract
	]
];


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

FiberTractography[tensor_, vox_, opts : OptionsPattern[]] := FiberTractography[tensor, vox, {1, {0, 1}}, opts]

FiberTractography[tensor_, vox_, {par_?ArrayQ, {min_?NumberQ, max_?NumberQ}}, opts : OptionsPattern[]] := FiberTractography[tensor, vox, {{par, {min, max}}}, opts]

FiberTractography[tensor_, vox_, inp : {{_, {_, _}} ...}, opts : OptionsPattern[]] := Block[{
		lmin, lmax, amax, maxSeed, flip, per, int, stopT, step, intF,
		tensTr, tensMask, dim, voxTr, inpTr, treshs, stop, range, coors, intE, intS,
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
	intF = Switch[OptionValue[Method], "RungeKutta", RK2, "RungeKutta4", RK4, _, Euler];
	
	(*prepare data*)
	tensTr = FlipTensorOrientation[tensor, per, flip];
	tensTr = FlipTensorOrientation[RotateDimensionsLeft[#] & /@ tensTr, {"y", "x", "z"}];
	tensMask = Unitize[Mean@Abs@tensTr];
	tensTr = RotateDimensionsLeft[tensTr];
	dim = Dimensions[tensTr][[;; -2]];
	voxTr = RotateLeft[vox, 1];
	
	(*convert dimensions of stoppings*)
	{inpTr, treshs} = Transpose[inp];
	inpTr = Transpose[{RotateDimensionsLeft[#] & /@ inpTr, treshs}];
	
	(*create stopping mask*)
	stop = tensMask (Times @@ (Mask[#[[1]], #[[2]]] & /@ inpTr));
	
	(*get the data coordinates on center of voxels*)
	range = Thread[{voxTr, voxTr dim}] - 0.5 voxTr;
	coors = Flatten[CoordinateBoundsArray[# + {0., 0.001} & /@ range, voxTr], 2];
	
	(*make the trac and stop interpolation functions*)
	intE = Interpolation[Thread[{coors, Flatten[tensTr, 2]}], InterpolationOrder -> int,
		"ExtrapolationHandler" -> {{0., 0., 0., 0., 0., 0.} &, "WarningMessage" -> False}, Method -> "Hermite"];
	intS = ListInterpolation[stop, range, InterpolationOrder -> Clip[int,{1,Infinity}], 
		"ExtrapolationHandler" -> {0. &, "WarningMessage" -> False}, Method -> "Hermite"];
	
	(*make the seed point*)
	ones = voxTr # & /@ SparseArray[stop]["NonzeroPositions"];
	(*n = Min[{maxSeed, Length[ones]}];*)
	seeds = RandomChoice[ones, maxSeed] + Transpose[voxTr RandomReal[{-0.5, 0.5} , {3, maxSeed}]];
	
	If[OptionValue[TracMonitor],Print["number of selected seedpoints: ", Length[seeds]]];
	
	t1 = First[AbsoluteTiming[
		(*perform tractography for each seed point*)
		iii=1;
		smax = Ceiling[(lmax/step)];
		tracts = Monitor[Table[
			TractFunc[seeds[[iii]], step, {amax, smax, stopT}, {intE, intS, intF}]
			, {iii, 1, Length[seeds]}]
		, ProgressIndicator[iii, {0, Length[seeds]}]];
	
	(*select only tracts within correct range and clip tracts that are longer*)
	tracts = If[step Length[#] > lmax, 
			drop = Ceiling[(step Length[#] - lmax)/(2 step)] + 1; #[[drop ;; -drop]], #
		] & /@ Select[tracts, lmin < step Length[#] &];
	]];
	
	(*report timing*)
	If[OptionValue[TracMonitor],
		Print["Tractography took ", t1, " seconds"];
		Print[Length[tracts], " tracts were generated with average length ", Round[step Mean[Length /@ tracts], 0.1],"\[PlusMinus]", Round[step StandardDeviation[Length /@ tracts], 0.1] , " mm"];
	];

	(*output tracts*)
	tracts
]


(* ::Subsubsection::Closed:: *)
(*Euler*)


Euler[y_, v_, h_, func_] := VecAlign[v, h func[y]]


(* ::Subsubsection::Closed:: *)
(*RK2*)


RK2[y_, v_, h_, func_] := Block[{k1},
	k1 = VecAlign[v, h func[y]];
	VecAlign[v, h*func[y + k1/2]]
]


(* ::Subsubsection::Closed:: *)
(*RK4*)


RK4[y_, v_, h_, func_] := Block[{k1, k2, k3, k4},
	k1 = VecAlign[v, h func[y]];
	k2 = VecAlign[v, h func[y + k1/2]];
	k3 = VecAlign[v, h func[y + k2/2]];
	k4 = VecAlign[v, h func[y + k3]];
	VecAlign[v, (k1/6 + k2/3 + k3/3 + k4/6)]
  ]


(* ::Subsubsection::Closed:: *)
(*VecAng*)


VecAng = Compile[{{v1, _Real, 1}, {v2, _Real, 1}}, Block[{v,n},
	n = Norm[v1] Norm[v2];
	If[n > 0.,
		v = v1 . v2/n;
		If[-1. < v < 1., 180 ArcCos[v]/Pi, 0.],
		0.
	]
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]
  

(* ::Subsubsection::Closed:: *)
(*VecAlign*)
  
  
VecAlign = Compile[{{v1, _Real, 1}, {v2, _Real, 1}}, Sign[Sign[v1 . v2]+0.1] v2,
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsubsection::Closed:: *)
(*EigV*)


(*fast eigenvector calculation for tractography*)
EigV[{0., 0., 0., 0., 0., 0.}] := {0., 0., 1.}
EigV[{xx_, yy_, zz_, xy_, xz_, yz_}] := First[Eigenvectors[{{xx, xy, xz}, {xy, yy, yz}, {xz, yz, zz}}]];


(* ::Subsubsection::Closed:: *)
(*TractFunc*)


TractFunc[loc0_, h_, {amax_, smax_, str_}, {trF_, stF_, intF_}] := Block[{step0, in1, in2},
	
	(*initial direction*)
	step0 = intF[loc0, {1,0,0}, h, Quiet[EigV[trF @@ #]] &];
	
	(*place two seed points arount initial points*)
	in1 = {loc0 + 0.5 step0, -step0};
	in2 = {loc0 - 0.5 step0, step0};
	
	(*tractography in two directions and *)
	Join[
		Reverse@TractFunci[in1, -h, {amax, smax, str}, {trF, stF, intF}][[2 ;;, 1]],
		TractFunci[in2, h, {amax, smax, str}, {trF, stF, intF}][[2 ;;, 1]]
	]
]


(* ::Subsubsection::Closed:: *)
(*TractFunci*)


TractFunci[{loc0_, step0_}, h_, {amax_, smax_, str_}, {trF_, stF_, intF_}] := Block[
	{i, ang, stop, loc, step, stepn, angt},

	i = ang = 0;
	stop = 1;
	
	NestWhileList[(
			i++;
			(*location from itteration n*)
			step = #[[2]];
			loc = #[[1]] + step;
			(*get the step itteration n*)
			stepn = intF[loc, step, h, EigV[trF @@ #] &];
			(*find direction change with previous step direction*)
			ang = VecAng[step, stepn];
			(*get location of itteration n+1 to see if outside int area*)
			stop = stF @@ (loc + stepn);
			(*output*)
			{loc, stepn}
		) &,
		{loc0, step0},
		(i <= smax && ang < amax && stop > str) &
	]
];


(* ::Subsubsection::Closed:: *)
(*FiberTractography*)


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
