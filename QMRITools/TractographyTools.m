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
FitTract[{tract, ...}] fits a list of tracts defined as a list of {x,y,z} coordinates with a polinomial function.
"

FiberTractography::usage =
"FiberTractography[tensor, vox] performs fibertractography on the tensor with voxels dimensions vox.
FiberTractography[tensor, vox, {par, {min, max}}] performs fibertractography on the tensor with voxels dimensions vox with additional stoppin criteria
par, where tracts are only generated between values of par min and max.
FiberTractography[tensor, vox, {{par, {min, max}}, ..}] performs fibertractography on the tensor with voxels dimensions vox with 
multiple additional stopping criteria
"

(* ::Subsection::Closed:: *)
(*Options*)


FittingOrder::usage = 
"FittingOrder is an option for FitTract. It specifies the polinominal order of the function to fit the tract."


FiberLengthRange::usage = 
"FiberLengthRange is an option for FiberTractography and specifies the allowed tract range."

FiberAngle::usage = 
"FiberAngle is an option for FiberTractography and specifies the allowed angle change per tract step."

SeedSpacing::usage =
"SeedSpacing is an option for FiberTractography and specifies the distance between the seedpoints in mm."

SeedRandomization::usage =
"SeedRandomization is an option for FiberTractography and specifies if there is a random translation of the seedpoints."

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
	SeedSpacing -> 2,
	SeedRandomization -> 0.5,
	TensorFilps -> {1, 1, 1},
	TensorPermutations -> {"x", "y", "z"},
	InterpolationOrder -> 1,
	StopThreshhold -> 0.5,
	StepSize -> 2,
	Method -> "Euler",
	MaxSeedPoints -> Infinity
};

SyntaxInformation[FiberTractography] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

FiberTractography[tensor_, vox_, opts : OptionsPattern[]] := FiberTractography[tensor, vox, {1, {0, 1}}, opts]

FiberTractography[tensor_, vox_, {par_?ArrayQ, {min_?NumberQ, max_?NumberQ}}, opts : OptionsPattern[]] := FiberTractography[tensor, vox, {{par, {min, max}}}, opts]

FiberTractography[tensor_, vox_, inp : {{_, {_, _}} ...}, opts : OptionsPattern[]] := Block[{
		lmin, lmax, amax, seedSp, maxSeed, seedRand, flip, per, int, stopT, step,
		tensTr, voxTr, inpTr, treshs, coors, ran, stop, intE, intS, intF,
		drop, sx, sy, sz, ex, ey, ez, seeds, tracts, iii, t1
	},
	(*get the options*)
	{lmin, lmax} = OptionValue[FiberLengthRange];
	amax = OptionValue[FiberAngle];
	seedSp = OptionValue[SeedSpacing];
	seedRand = Clip[OptionValue[SeedRandomization], {0, 1}];
	maxSeed = OptionValue[MaxSeedPoints];
	
	flip = OptionValue[TensorFilps];
	per = OptionValue[TensorPermutations];
	
	int = OptionValue[InterpolationOrder];
	stopT = OptionValue[StopThreshhold];
	
	step = OptionValue[StepSize];
	intF = Switch[OptionValue[Method], "RungeKutta", RK2, "RungeKutta4", RK4, _, Euler];
	
	(*prepare data*)
	tensTr = FlipTensorOrientation[TransData[#, "l"] & /@ tensor, per(*{"y","x","z"}*), flip];
	voxTr = RotateLeft[vox, 1];
	
	(*convert dimensions of stoppings*)
	{inpTr, treshs} = Transpose[inp];
	inpTr = Transpose[{TransData[#, "l"] & /@ inpTr, treshs}];
	
	(*get the data coordinates on center of voxels*)
	coors = Flatten[Array[voxTr {#1, #2, #3} - 0.5 voxTr &, Dimensions[tensTr][[2 ;;]]], 2];
	
	(*create stopping mask*)
	stop = Unitize[Mean[Abs[tensTr[[1 ;; 3]]]]] (Times @@ (Mask[#[[1]], #[[2]]] & /@ inpTr));
	
	
	(*make the trac and stop interpolation functions*)
	intE = Interpolation[Thread[{coors, Flatten[TransData[tensTr, "l"], 2]}], InterpolationOrder -> int];
	intS = Interpolation[Thread[{coors, Flatten[stop, 2]}], InterpolationOrder -> int];
	
	(*make the seed point*)
	{sx, sy, sz} = 0.5 voxTr;
	{ex, ey, ez} = voxTr Dimensions[tensTr][[2 ;;]] - 0.5 vox;
	
	seeds = Flatten[Table[{x, y, z}, {x, sx, ex, seedSp}, {y, sy, ey, seedSp}, {z, sz, ez, seedSp}], 2];
	seeds = seeds + RandomReal[NormalDistribution[0, seedRand], Dimensions[seeds]];
	seeds = Quiet@Pick[seeds, UnitStep[(intS @@ Transpose[seeds]) - 0.5], 1];
	seeds = RandomSample[seeds, Min[{maxSeed, Length[seeds]}]];
	
	Print["number of selected seedpoints: ", Length[seeds]];
	
	t1 = First[AbsoluteTiming[
		(*perform tractography for each seed point*)
		tracts = Monitor[Table[
			TractFunc[seeds[[iii]], step, {amax, lmax, stopT}, {intE, intS, intF}]
			, {iii, 1, Length[seeds]}]
		, ProgressIndicator[iii, {0, Length[seeds]}]];
	
	(*select only tracts within correct range and clip tracts that are longer*)
	tracts = If[step Length[#] > lmax, 
			drop = Ceiling[(step Length[#] - lmax)/(2 step)] + 1; #[[drop ;; -drop]], #
		] & /@ Select[tracts, lmin < step Length[#] &];
	]];
	
	(*report timing*)
	Print["Tractography took ", t1, " seconds"];
	Print[Length[tracts], " tracts were generated with average length ", step Round[Mean[Length /@ tracts], .1], "mm"];

	(*output tracts*)
	tracts
]


(* ::Subsubsection::Closed:: *)
(*Euler*)


Euler[y_, h_, func_] := h Normalize[func[y]]


(* ::Subsubsection::Closed:: *)
(*RK2*)


RK2[y_, h_, func_] := h Normalize[func[y + h func[y]/2]]


(* ::Subsubsection::Closed:: *)
(*RK4*)


RK4[y_, h_, func_] := Module[{k1, k2, k3, k4},
	k1 = func[y];
	k2 = func[y + h k1/2];
	k3 = func[y + h k2/2];
	k4 = func[y + h k3];
	h Normalize[(k1/6 + k2/3 + k3/3 + k4/6)]
]


(* ::Subsubsection::Closed:: *)
(*EigV*)


(*fast eigenvector calculation for tractography*)
EigV[{0., 0., 0., 0., 0., 0.}] := {0., 0., 1.}
EigV[t_] := Eigenvectors[{{t[[1]], t[[4]], t[[5]]}, {t[[4]], t[[2]], t[[6]]}, {t[[5]], t[[6]], t[[3]]}}][[1]]


(* ::Subsubsection::Closed:: *)
(*TractFunc*)


TractFunc[loc0_, h_, {amax_, smax_, str_}, {trF_, stF_, intF_}] := Block[{step0, in1, in2},
	
	(*initial direction*)
	step0 = intF[loc0, h, Quiet[EigV[trF @@ #]] &];
	
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


TractFunci[{loc0_, step0_}, h_, {amax_, smax_, str_}, {trF_, stF_, intF_}] := Block[{
		i, met, step02, loc, locn, aMax, step, stepn, ang, angt, max, func, stop,
		xmin, xmax, ymin, ymax, zmin, zmax, stopFunc, trFunc
	},
	
	(*make tract and stop function*)
	{{xmin, xmax}, {ymin, ymax}, {zmin, zmax}} = trF[[1]];
	trFunc = Quiet[EigV[trF @@ #]] &;
	stopFunc = Quiet[stF @@ #] &;
	
	(*get input and constraints *)
	max = Ceiling[Abs[(smax/h)]];
	aMax = Abs[h amax];(*angle per mm*)
	i = ang = 0;
	stop = 1;
	
	(*inital values*)
	loc = locn = loc0;
	
	NestWhileList[(
			i++;
			(*location and step from itteration n-1*)
			loc = #[[1]];
			step = #[[2]];
		
			(*location from itteration n*)
			loc = loc + step;
		
			(*get the step and location from itteration n*)
			stepn = intF[loc, h, trFunc];
		
			(*find direction change with previous step direction*)
			angt = N@VectorAngle[step, stepn]/Degree;
		
			(*prevent step from going backwards*)
			stepn = If[angt < 90, stepn, -stepn];
			ang = If[angt < 90, angt, 180 - angt];
		
			(*get location of itteration n+1 to see if outside int area*)
			locn = loc + stepn;
			stop = stopFunc[locn];
		
			(*output*)
			{loc, stepn, ang}
		) &, {loc0, step0, ang}, 
		(xmin <= locn[[1]] <= xmax && ymin <= locn[[2]] <= ymax && zmin <= locn[[3]] <= zmax && i <= max && ang < amax && stop > str) &
	]
];



(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
