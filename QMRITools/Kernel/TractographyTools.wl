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

FitTract::usage = 
"FitTract[tract] fits a tract defined as a list of {x,y,z} coordinates with a polinomial function.
FitTract[{tract, ...}] fits a list of tracts defined as a list of {x,y,z} coordinates with a polinomial function."

PlotTracts::usage = 
"PlotTracts[tracts, vox] ...
PlotTracts[tracts, vox, dim] ..."


FindTensorPermutation::usage = 
"FindTensorPermutation[tensor, vox] performs tractography for all tensor permutations and gives back the one that has the longest tracts.
FindTensorPermutation[tensor, vox, {par, {min, max}}] same but with additional stoppin criteria par, where tracts are only generated between values of par min and max.
FindTensorPermutation[tensor, vox, {{par, {min, max}}, ..}] same but with with multiple additional stopping criteria.

Ouput = {permutations, flips, plot}

FindTensorPermutation[] is based on DOI: 10.1016/j.media.2014.05.012."


TractDensityMap::usage = 
"TractDensityMap[tracts, vox, dim] ..."

SeedDensityMap::usage = 
"SeedDensityMap[seeds, vox, dim] ..."

TractLengthMap::usage = 
"TractLengthMap[tracts, vox, dim] ..."

TractAngleMap::usage = 
"TractAngleMap[tracts, vox, dim] ..."


FilterTracts::usage = 
"FilterTracts[tracts, vox, select] ..."

SelectTractTroughPlane::usage = 
"SelectTractTroughPlane ..."

SelectTractTroughVol::usage = 
"SelectTractTroughVol ..."

SelectTractInVol::usage = 
"SelectTractInVol ..."

SelectTractPartInVol::usage = 
"SelectTractPartInVol ..." 

PartTracts::usage = 
"PartTracts ..."

SelectTracts::usage = 
"SelectTracts ..."

CombineROIs::usage = 
"CombineROIs ..."


FiberLength::usage = 
"FiberLength[tracts] ..."

GetTractValues::usage = 
"GetTractValues[tracts, val, vox, int] ..."


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

MaxTracts::usage = 
"MaxTracts ..."

TractColoring::usage = 
"TractColoring ..."


NormalizeDensity::usage = 
"NormalizeDensity is an option for TractDensityMap ..."



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
	StepSize -> Automatic,
	Method -> "Euler",
	MaxSeedPoints -> Automatic,
	TracMonitor->True
};

SyntaxInformation[FiberTractography] = {"ArgumentsPattern" -> {_, _, _., OptionsPattern[]}};

FiberTractography[tensor_, voxi_, opts : OptionsPattern[]] := FiberTractography[tensor, voxi, {0. First@tensor, {0., 1.}}, opts]

FiberTractography[tensor_, voxi_, {par_?ArrayQ, {min_?NumberQ, max_?NumberQ}}, opts : OptionsPattern[]] := FiberTractography[tensor, voxi, {{par, {min, max}}}, opts]

FiberTractography[tensor_, vox_, inp : {{_, {_, _}} ...}, OptionsPattern[]] := Block[{
		lmin, lmax, amax, maxSeed, flip, per, int, stopT, step, tracF, trFunc,
		tens, tensMask, inpTr, treshs, stop, coors, intE, intS, ones, n, 
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
 	 	
	(*prepare tensor and stop data*)
	(*tens = FlipTensorOrientation[tensor, per /. Thread[{"x","y","z"} -> {"z","y","x"}], Reverse[flip]];*)
	tens = FlipTensorOrientation[tensor, per, flip];
	stop = SparseArray@Unitize[Total@Abs@tens] Mask[inp];
			
	(*make the trac and stop interpolation functions*)
	intE = MakeIntFunction[RotateDimensionsLeft[{tens}, 2], vox, int];
	intS = MakeIntFunction[Normal@stop, vox, int];
	
	(*make the random seed points*)
	SeedRandom[1234];
	seedN = Total@Flatten@stop;
	seedI = stop["NonzeroPositions"];
	seeds = {};
	maxSeed = If[NumberQ[maxSeed],maxSeed,seedN];
	
	While[Length[seeds] < maxSeed,
		seedT = Flatten[RandomSample[seedI, #] & /@ Block[{i = maxSeed},
			First@Last@Reap[Until[i < 0, Sow[If[i > seedN, seedN, i]];
			i = i - seedN;]]], 1];
		seedT = Threaded[vox] (seedT + RandomReal[{-0.99, 0}, {maxSeed, 3}]);
		seedT = Pick[seedT, intS @@@ seedT, 1.];
		
		seeds = Join[seeds, seedT];
		seeds = seeds[[;; Min[{Length@seeds, maxSeed}]]]
	];
	seedN = Length@seeds;

	If[OptionValue[TracMonitor],
		Print["Starting tractography for ", seedN, " seed points with stepsize ", step, " mm"];
	];

	(*perform tractography for each seed point*)	
	t1 = First[AbsoluteTiming[
		trFunc = TractFunc[#, step, {amax, smax, stopT}, {intE, intS, tracF}]&;
		tracts = If[OptionValue[TracMonitor],
			Monitor[Table[trFunc@seeds[[iii]], {iii, 1, seedN}], ProgressIndicator[iii, {0, seedN}]],
			Table[trFunc@seeds[[iii]], {iii, 1, seedN}]];
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
		Print["Tractography took ", t1, " seconds"];
		Print[Length[tracts], " valid tracts with length ", Round[step Mean[Length /@ tracts], 0.1],"\[PlusMinus]", Round[step StandardDeviation[Length /@ tracts], 0.1] , " mm"];	
	];

	(*output tracts*)
	{tracts, seeds}
]


(* ::Subsubsection::Closed:: *)
(*TractFunc*)


TractFunc[loc0_?VectorQ, h_, stp_, fun_] := Block[{dir0},
	dir0 = h EigVec[First[fun]@@loc0];
	(*tractography in the two directions around seed point*)
	ToPackedArray@Join[
		Reverse@TractFunc[{loc0 + dir0/2, -dir0}, h, stp, fun], 
		TractFunc[{loc0 - dir0/2, dir0}, h, stp, fun]
	]
]

TractFunc[{loc0_?VectorQ, step0_?VectorQ}, h_, {amax_, smax_, stoptr_}, {intE_, intS_, tracF_}] := Block[
	{ang, stop, loc, locn, step, stepn, angt,out},
	
	ang = 0;
	stop = intS@@(loc0 + step0);
	(*perform tractography*)
	out = NestWhileList[(
		(*get new location and the step at new location*)
		{loc, step} = #;
		locn = loc + step;
		stepn = tracF[locn, step, h, EigVec[intE@@#]&];
		(*get stop criteria*)
		ang = VecAng[step, stepn];
		stop = (intS@@(locn + stepn));			
		(*output*)
		{locn, stepn})&, {loc0, step0}, 
		(ang < amax && stop > stoptr)&, 1, smax
	];
	If[Length[out] > 1, out[[2;;, 1]], {}]
];


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
	k2 = VecAlign[v, h func[y + k1/2]];
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
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


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


EigVec = Compile[{{tens, _Real, 1}}, Block[{
	dxx, dyy, dzz, dxy, dxz, dyz, i1, i2, i3, i13, phi, v1, v, s, l1, a, b, c},
	(*method from 10.1016/j.mri.2009.10.001*)
	If[Total[tens] === 0.,
		{0., 0., 1.},
		
		{dxx, dyy, dzz, dxy, dxz, dyz} = tens;
		
		i1 = dxx + dyy + dzz;
		i2 = (dxx dyy + dxx dzz + dyy dzz) - (dxy^2 + dxz^2 + dyz^2);
		i3 = dxx dyy dzz + 2 dxy dxz dyz - (dzz (dxy)^2 + dyy (dxz)^2 + dxx (dyz)^2);
		
		i13 = i1/3;
		v = i13^2 - i2/3 + .1^30;
		v1 = 1./v;
		s = (i13)^3 - i1 i2/6 + i3/2;
		phi = ArcCos[(v1 s) Sqrt[v1]]/3;
		
		l1 = i13 + 2 Sqrt[v] Cos[phi];
		
		{a, b, c} = {dxz dxy, dxy dyz, dxz dyz} - {dyz, dxz, dxy} ({dxx, dyy, dzz} - l1);
		
		Normalize[{b c, a c, a b}]
	]
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


EigVec2 = Compile[{{tens, _Complex, 1}}, Block[{
	def, xx, yy, zz, xy, xz, yz, tr, xx2, yy2, zz2, xy2, xz2, yz2, tr2, 
	p, b, b2, c, d, e, f, e1, u1, u11, u12, u13, u2, u21, u22, u23, u3, u1n, u3n},
	
	def = {0., 0., 1.};
	
	(*tens and tensor squared indices*)
	{xx, yy, zz, xy, xz, yz} = tens;
	tr = xx + yy + zz;
	
	Re@If[tr <= 0., def,
		(*root finding parameters*)
		{xx2, yy2, zz2, xy2, xz2, yz2} = tens^2;
		tr2 = xx2 + yy2 + zz2 + 2 xy2 + 2 xz2 + 2 yz2;
				
		p = 1/3;
		b = tr; b2 = b^2;
		c = 0.5 (tr2 - tr^2);
		d = xx yy zz - xx yz2 - yy xz2 - zz xy2 + 2 xy xz yz;
		
		e = Sqrt[-3 c^2 (b2 + 4 c) + 6 b (2 b2 + 9 c) d + 81 d^2];
		f = (2 b^3 + 9 b c + 3 (9 d + e))^p;
		
		If[f == 0., def,
			(*first root is first eigenvalue*) 
			e1 = Abs[(1/6) (2 b + (2^(1 + p) (b2 + 3 c))/f + 4^p f)];
		
			(*QR decomp Using the Gram\Schmidt process to find first vector*)
			u1 = {u11, u12, u13} = {xx - e1, xy, xz};
			u2 = {xy, yy - e1, yz};
			
			u1n = u1 . u1;
			If[u1n <= 0., def,
				{u21, u22, u23} = u2 - ((u1 . u2)/u1n) u1;
				u3 = {-u13 u22 + u12 u23, u13 u21 - u11 u23, -u12 u21 + u11 u22};
				u3n = u3 . u3;
				If[u3n <= 0., def, Sign[u3[[1]]] (u3/Sqrt[u3n])]
			]
		]
	]		
], RuntimeAttributes->{Listable}, RuntimeOptions->"Speed"];


(* ::Subsection::Closed:: *)
(*FitTract*)


Options[FitTract] = {FittingOrder -> 4};

SyntaxInformation[FitTract] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};

FitTract[tract_, OptionsPattern[]] := ToPackedArray/@FitTractC[tract, Round@OptionValue[FittingOrder]]

FitTractC = Compile[{{trf, _Real, 2}, {ord, _Real, 0}}, Block[{mat = #^Range[0., ord] & /@ Range[Length[trf]]},
	mat . (PseudoInverse[mat] . trf)
], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsection:: *)
(*PlotTracts*)


(* ::Subsubsection::Closed:: *)
(*PlotTracts*)


Options[PlotTracts] = {
	MaxTracts -> 2000, 
	ImageSize -> 800, 
	Method->"line", 
	TractColoring->"Direction",
	ColorFunction->"SouthwestColors",
	Boxed->True
}

PlotTracts[tracts_, voxi_, opts : OptionsPattern[]] := PlotTracts[tracts, voxi, 0, opts]

PlotTracts[tracts_, voxi_, dimi_, OptionsPattern[]] := Block[{
	range, vox, size, select, opts, col, tube, line, plot, colOpt, met, set, max, colf
	},
	
	vox = Reverse@voxi;
	
	range = If[dimi === 0,
		Round[Reverse[MinMax /@ Transpose@Flatten[tracts, 1]]/vox],
		Reverse@Thread[{{0, 0, 0}, dimi}]
	];
	
	size = vox Flatten[Differences /@ range];
	
	max = OptionValue[MaxTracts];
	If[OptionValue[Method]==="tube", max =  Min[4000, max]];
	select = ToPackedArray /@ Map[Reverse, RandomSample[tracts, Min[max, Length[tracts]]], {-2}];
	
	opts = Sequence[{
		Method -> {"TubePoints" -> {6, 2}}, Lighting -> "Accent", 
		ImageSize -> OptionValue[ImageSize], SphericalRegion -> True, Boxed -> OptionValue[Boxed],
		Background -> Lighter@Gray, BoxRatios -> size, PlotRange -> range, 
		Axes -> OptionValue[Boxed], LabelStyle -> Directive[{Bold, 16, White}]
	}];
	
	colf = OptionValue[ColorFunction];
	colOpt = OptionValue[TractColoring];
	{met, set} = If[ListQ[colOpt], colOpt, {colOpt, Automatic}];
	
	col = Switch[met,
		"Length",
		MakeLengthColor[select, set, colf],
		"Angle",
		MakeAngleColor[select, set, colf],
		"Direction",
		MakeColor[select]
	];
	
	plot = Switch[OptionValue[Method],
		"tube", Scale[Tube[select, 0.5 Min[vox], VertexColors -> col], 1/size, {0, 0, 0}],
		"line", Scale[Line[select, VertexColors -> col], 1/vox, {0, 0, 0}],
		_, $Failed
	];
	
	Graphics3D[{CapForm["Square"], JoinForm["Miter"], plot}, opts]
]


(* ::Subsubsection::Closed:: *)
(*Make Color*)


MakeColor[tract : {{_?NumberQ, _?NumberQ, _?NumberQ} ..}] := Block[{dirs},
	dirs = Abs[Differences[tract]];
	ToPackedArray[Normalize[#] & /@ Mean[{Prepend[dirs, dirs[[1]]], Append[dirs, dirs[[-1]]]}]]
]

MakeColor[tracts : {_?ListQ ..}] := MakeColor /@ tracts

MakeColor[tract_, fun_, ran_] := Block[{prange, colFunc, cval, fore},
	prange = ran;
	colFunc = ColorData["ThermometerColors"];
	cval = fun @@ # & /@ tract;
	fore = Unitize[cval];
	(fore (colFunc /@ Rescale[cval, prange])) + Gray (1 - fore)
]


(* ::Subsubsection::Closed:: *)
(*MakeLengthColor*)


MakeLengthColor[tracts_, set_, colf_]:=Block[{len, sc},
	len = FiberLength[tracts];
	sc = If[set===Automatic, Quantile[len, {.05,0.95}], set];
	MapThread[ConstantArray[#2, Length@#1]&, {tracts, (ColorData[colf]/@Rescale[len, sc])}]
];


(* ::Subsubsection::Closed:: *)
(*MakeAngColor*)


MakeAngleColor[tracts_, set_, colf_] := Block[{ang, sc},
	ang = VecAngC[tracts];
	sc = If[set === Automatic, {0, 90}, set];
	ang = Rescale[(Mean[{Prepend[#, #[[1]]], Append[#, #[[-1]]]}]) & /@ ang, sc];
	ToPackedArray[ColorData[colf] /@ #] & /@ ang
];


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

FindTensorPermutation[tensor_, vox_, opts : OptionsPattern[]] := FindTensorPermutation[tensor, vox, {0. First@tensor, {0., 1.}}, opts]

FindTensorPermutation[tensor_, vox_, {par_?ArrayQ, {min_?NumberQ, max_?NumberQ}}, opts : OptionsPattern[]] := FindTensorPermutation[tensor, vox, {{par, {min, max}}}, opts]

FindTensorPermutation[tens_, vox_, stop_, opts : OptionsPattern[]] := Block[{
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


(* ::Subsection:: *)
(*DensityMap*)


(* ::Subsubsection::Closed:: *)
(*ThreadedDiv*)


ThreadedDiv = Compile[{{coor, _Real, 1}, {vox, _Real, 1}}, 
	Ceiling[coor/vox],
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


ThreadedDivD = Compile[{{coor, _Real, 2}, {vox, _Real, 1}}, 
	Ceiling[#/vox] & /@ Mean[{coor[[2 ;; -1]], coor[[1 ;; -2]]}], 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


ThreadedDivL = Compile[{{coor, _Real, 2}, {vox, _Real, 1}}, 
	Ceiling[#/vox] & /@ coor, 
RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"]


(* ::Subsubsection::Closed:: *)
(*MakeList*)


MakeList[len_] := Block[{ran},
	ran = Range[0, len, 50000];
	{1, 0} + # & /@ Partition[If[Last@ran =!= len, Append[ran, len], ran], 2, 1]
]


(* ::Subsubsection::Closed:: *)
(*SeedDensityMap*)


SyntaxInformation[SeedDensityMap] = {"ArgumentsPattern" -> {_, _, _}};

SeedDensityMap[seeds_, vox_, dim_] := Normal@SparseArray[Normal@Counts@ThreadedDiv[seeds, vox], dim];


(* ::Subsubsection::Closed:: *)
(*TractDensityMap*)


Options[TractDensityMap] = {NormalizeDensity -> True}

SyntaxInformation[TractDensityMap] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

TractDensityMap[tracts_, vox_, dim_,  OptionsPattern[]] := Block[{dens},
	dens = ConstantArray[0., dim];
	Table[dens += Normal@SparseArray[Select[
		Normal[Counts[Flatten[DeleteDuplicates /@ ThreadedDivL[tracts[[i[[1]] ;; i[[2]]]], vox], 1]]],
		(1 <= #[[1, 1]] <= dim[[1]] && 1 <= #[[1, 2]] <= dim[[2]] && 1 <= #[[1, 3]] <= dim[[3]]) &], dim],
	{i, MakeList[Length[tracts]]}];
	
	If[OptionValue[NormalizeDensity], dens/MedianNoZero[Flatten[dens]], dens]
]


(* ::Subsubsection::Closed:: *)
(*TractLengthMap*)


SyntaxInformation[TractLengthMap] = {"ArgumentsPattern" -> {_, _, _}};

TractLengthMap[tracts_, vox_, dim_] := Block[{tr, vals, lengs},
	vals = GatherBy[Flatten[Table[
		tr = tracts[[i[[1]] ;; i[[2]]]];
		lengs = GatherBy[Flatten[Thread /@ Thread[(DeleteDuplicates /@ ThreadedDivL[tr, vox]) -> FiberLength@tr], 1], First];
		Select[Thread[lengs[[All, 1, 1]] -> lengs[[All, All, 2]]], (1 <= #[[1, 1]] <= dim[[1]] && 1 <= #[[1, 2]] <= dim[[2]] && 1 <= #[[1, 3]] <= dim[[3]]) &]
    , {i, MakeList[Length[tracts]]}]], First];
	
	vals = Thread[vals[[All, 1, 1]] -> vals[[All, All, 2]]];
	Normal@SparseArray[Thread[vals[[All, 1]] -> (Median[Flatten[#]] & /@ vals[[All, 2]])], dim]
]


(* ::Subsubsection::Closed:: *)
(*TractAngleMap*)


SyntaxInformation[TractAngleMap] = {"ArgumentsPattern" -> {_, _, _}};

TractAngleMap[tracts_, vox_, dim_] := Block[{tr, angsCoor, vals},
	vals = GatherBy[Flatten[Table[
		tr = tracts[[i[[1]] ;; i[[2]]]];
		angsCoor = GatherBy[Flatten[Thread /@ Thread[ThreadedDivD[tr, vox] -> VecAngC[tr]], 1], First];
		vals = Select[Thread[angsCoor[[All, 1, 1]] -> angsCoor[[All, All, 2]]], (1 <= #[[1, 1]] <= dim[[1]] && 1 <= #[[1, 2]] <= dim[[2]] && 1 <= #[[1, 3]] <= dim[[3]]) &]
	, {i, MakeList[Length[tracts]]}]], First];
	vals = Thread[vals[[All, 1, 1]] -> vals[[All, All, 2]]];
	Normal@SparseArray[Thread[vals[[All, 1]] -> (Median[Flatten[#]] & /@ vals[[All, 2]])],dim]
]


VecAngC = Compile[{{tr, _Real, 2}},
	90 - ((180./Pi) ArcCos[{1, 0, 0}.Normalize[Sign[#[[1]]] #]] & /@ (tr[[2 ;; -1]] - tr[[1 ;; -2]]))
, RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*FiberLength*)


FiberLength[tracts_]:=FLengthC[tracts]

FLengthC = Compile[{{trc,_Real,2}},
	Total[Norm/@Differences[trc]]
,RuntimeAttributes->{Listable},RuntimeOptions->"Speed"]


(* ::Subsubsection::Closed:: *)
(*GetTractValues*)


GetTractValues[tracts_, val_, vox_, int_:1]:=MakeIntFunction[val,vox,int]@@@#&/@tracts;


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
