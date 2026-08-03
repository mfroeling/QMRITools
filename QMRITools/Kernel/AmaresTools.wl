(* ::Package:: *)

(* ::Title:: *)
(*QMRITools AmaresTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`AmaresTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`AmaresTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


EditAmaresMatrix::usage
"EditAmaresMatrix[matrix_Association] ..."

ToFullAmaresMatrix::usage = 
"ToFullAmaresMatrix[names_List] ..."

FromFullAmaresMatrix::usage = 
"FromFullAmaresMatrix[matrix] ..."

AmaresBasis::usage = 
"AmaresBasis[settings_, {metabolites_, split_}] {nsamp, bw, field, nuc, te, f0} = settings; ..."

AmaresInitialValues::usage = 
"AmaresInitialValues[dat_, time_, ppm_, gyro_, basis_, matrix_, pl_] ..."

AmaresFitFidI::usage = 
"AmaresFitFidI[dat_, time_, gyro_, basis_, matrix_?AssociationQ, cons_?MatrixQ]"

ShowParameterMatrix::usage = 
"ShowParameterMatrix[fit_, mets_]"

AmaresSignalCalc::usage = 
"AmaresSignalCalc[bs_, time_, pars_, t1b1_, b1Par_, gyro_] ..."

AmaresSignalJacobianCalc::usage = 
"AmaresSignalJacobianCalc[bs_, time_, pars_, t1b1_, b1Par_, gyro_] ..."

MakeT1Model::usage = "
MakeT1Model[names_List, matrix_]"

FitB1FromATP::usage = ""

SaturationProfileCorrection::usage = ""

(* ::Subsection:: *)
(*Options*)


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*Support functions*)


ShowParameterMatrix[fit_, mets_] := Block[{fiti=fit},
	fiti[[1]] = fiti[[1]]/fiti[[1, 8]];
	MatrixForm[Transpose@Prepend[Transpose[Prepend[Round[fiti, .01], mets]], 
	{"met", "amp", "sigma", "gamma", "eps", "phi", "B1"}[[;; Length[fiti] + 1]]]]
]


ComplexToVector = Join @@ Through[{Re, Im}[#]] &;


VectorToComplex = {1, I} . Partition[#, Length[#]/2] &;


(* ::Subsection::Closed:: *)
(*AmaresBasis*)


AmaresBasis[settings_, metabolites_?VectorQ] := AmaresBasis[settings, {metabolites, {}}]

AmaresBasis[settings_, {metabolites_?ListQ, split_?ListQ}] := Block[{
		nsamp, bw, field, nuc, te, f0, dw, gyro, time, ppm, names, fids, basisFunc, sp, spin
	},

	{nsamp, bw, field, nuc, te, f0} = settings;
	dw = 1./bw;
	gyro = GetGyro[field, nuc];

	{time, ppm} = GetTimePpmRange[nsamp, dw, gyro] + {te, 0};
	{names, fids, sp, spin} = GetSpectraBasisFunctions[metabolites, split,
		BasisSequence -> {"PulseAcquire", te},
		SpectraSamples -> nsamp, SpectraBandwidth -> bw,
		SpectraPpmShift -> f0, SpectraFieldStrength -> field,
		SpectraNucleus -> nuc];
	basisFunc = Association[Thread[names -> fids]];

	{basisFunc, {time, ppm, gyro}}
]


(* ::Subsection:: *)
(*AmaresFitMatrix*)


fitMatrixDefaults = <|"amp" -> -1, "sigma" -> 0, "gamma" -> 0, "shift" -> -1, "phase" -> 0, "B1" -> False, "init" -> False|>;
fitMatrixParams = Keys[fitMatrixDefaults];


(* ::Subsubsection::Closed:: *)
(*ToFullAmaresMatrix*)


ToFullAmaresMatrix[names_List] := Association[# -> fitMatrixDefaults & /@ names];

ToFullAmaresMatrix[matrix_Association] := Association[# -> Join[fitMatrixDefaults, matrix[#]] & /@ Keys[matrix]];

ToFullAmaresMatrix[matrix_Association, names_List] := Association[# -> Join[fitMatrixDefaults, Lookup[matrix, #, <||>]] & /@ names];


(* ::Subsubsection::Closed:: *)
(*FromFullAmaresMatrix*)


FromFullAmaresMatrix[matrix_Association] := Select[
	Association[# -> Association[
		KeyValueMap[Function[{k, v}, If[v === fitMatrixDefaults[k], Nothing, k -> v]], matrix[#]]
		] & /@ Keys[matrix]],
	# =!= <||> &]


(* ::Subsubsection::Closed:: *)
(*EditAmaresMatrix*)


EditAmaresMatrix[names_List] := EditAmaresMatrix[<||>, names];

EditAmaresMatrix[matrix_Association] := EditAmaresMatrix[matrix, Keys[matrix]];

EditAmaresMatrix[matrix_Association, names_List] := Block[{start, return, list, nn, np, cellVars, cellStart, fitPar},
	fitPar = fitMatrixParams;
	list = Join[{-1 -> "unique", 0 -> "same"}, # -> "group-" <> IntegerString[#, 10, 2] & /@ Range[1, Length[names]]];
	start = ToFullAmaresMatrix[matrix, names];
	{nn, np} = {Length[names], Length[fitPar]};
	cellStart = Table[start[names[[i]], fitPar[[j]]], {i, nn}, {j, np}];

	return = DialogInput[DynamicModule[{data = cellStart, datai = cellStart},
	Column[{Grid[Prepend[Table[Prepend[Table[With[{i = i, j = j, param = fitPar[[j]]},
			If[param === "init" || param === "B1",
			Checkbox[Dynamic[data[[i, j]]]],
			PopupMenu[Dynamic[data[[i, j]]], list]
			]], {j, np}], names[[i]]], {i, nn}], 
		Prepend[fitPar, "Metabolite"]]
		, Frame -> All, Alignment -> Center, Spacings -> {1, 1}],
		Row[{
		Button["Save", DialogReturn[data], ImageSize -> 80],
		Spacer[10],
		Button["Cancel", DialogReturn[datai], ImageSize -> 80]
		}]}]
	]];
	Association[Thread[names -> (Association[Thread[fitPar -> #]] & /@ return)]]
]


(* ::Subsection::Closed:: *)
(*AmaresInitialValues*)


AmaresInitialValues[dat_, time_, basis_, matrix_, {gyro_, dw_}] := AmaresInitialValues[dat, time, basis, matrix, {gyro, dw}, False]

AmaresInitialValues[dat_, time_, basis_, matrix_, {gyro_, dw_}, pl_] := Block[{
		fids, base, gsm, ep, sig, bspec, spec, corr, wgth, pos, eps, model, 
		a, g, t, s, p, e, fit, fun ,gm, sm, aa, ph, ppm,
		modelC, objC, decay, amps, init,pEps, pGam, pPh, pAmp, plots
	},

	(*get basis functions and base for calibration*)
	fids = Values[basis];
	base = Total[Pick[Values[basis], Table[matrix[name, "init"], {name, Keys[basis]}]]];
	(*fit contraints for linewdith and shift and blur kernel*)
	{gsm, ep, sig} = {50, 0.2, 5};

	(*find initial shift using autocorrelation*)
	bspec = Rescale@Abs@ShiftedFourier[base];
	spec = Rescale@Abs@ShiftedFourier[dat];
	corr = ListCorrelate[bspec, spec, Round[Length[base]/2], 0];

	(*find max peak closest to 0*)
	ppm = GetPpmRange[dat, dw, gyro];
	wgth = 1/(Exp[(ppm^2/(2 sig^2))] (Sqrt[2 Pi] sig));
	pos = Position[wgth corr, Max[wgth corr]][[1, 1]];
	eps = ppm[[pos]];

	(*fit the signal appodization*)
	model = a Exp[-Pi (g t + (s t)^2)];
	fit = FindFit[Transpose[{time, Abs[Exp[2 Pi eps gyro I time] dat]}],
		{model, 0 <= g <= gsm, 0 <= s <= gsm}, {a, g, s}, t];
	fun = Function[{t}, Evaluate[model /. fit]];
	{gm, sm} = {g, s} /. fit;

	(*fit the phase using the basis function that are shifted and apodized*)
	modelC[a_, p_, e_, g_, s_, t_, b_] := a Exp[-p I] (Exp[-Pi (g t + (s t)^2)] Exp[2 Pi e gyro I t]) b; 
	objC[dd_, a_?NumericQ, p_, e_, g_, s_, t_, b_] := Total[Abs[dd - modelC[a, p, e, g, s, t, b]]^2];
	fit = FindMinimum[{objC[dat, a, p, e, g, s, time, base],
		-2 Pi <= p <= 2 Pi, 10 <= g <= gsm, 10 <= s < gsm, 
		eps - ep < e < eps + ep},
	{a, {p, 0}, {e, eps}, {g, gm}, {s, sm}}, MaxIterations -> 100];
	{aa, sm, gm, eps, ph} = {a, s, g, e, p} /. Last[fit];

	(*guess the amps with least squares fitting*)
	decay = Exp[-ph I] (Exp[-Pi (gm time + (sm time)^2) + 2 Pi eps gyro I time]);
	amps = Abs[LeastSquares[Transpose[fids] decay, dat]];

	(*generate parameter matrix*)
	init = Join[{amps}, ConstantArray[#, Length@amps] & /@ {sm, gm, eps, ph, 1}];

	If[pl,
		(*make autocorrelation test plot*)
		pEps = ListLinePlot[
			Transpose[{ppm, #}] & /@ {spec, Rescale@Abs@ShiftedFourier[Exp[2 Pi eps gyro I time] base],	
				Rescale[wgth corr]}, 
			ScalingFunctions -> {"Reverse", None},
			ImageSize -> 200, Frame -> True, Axes -> True, 
			FrameTicks -> {{False, False}, {True, False}},
			PlotRange -> {Full, {0, 1.2}}, PlotLabel -> "Estimation of shift", 
			PlotHighlighting -> None];
		(*make the gamma fit test plot*)
		pGam = ListLinePlot[
			Transpose[{1000 time, Abs@#}] & /@ {dat, fun[time]},
			PlotRange -> Full, ImageSize -> 200, Frame -> True, Axes -> True, 
			FrameTicks -> {{False, False}, {True, False}},
			PlotLabel -> "Estimation of linewidth", PlotHighlighting -> None];
		(*make the phase test plot*)
		pPh = ListLinePlot[
			Transpose[{1000 time, #}] & /@ {Re@dat, Re[modelC[aa, ph, eps, gm, sm, time, base]]},
			PlotRange -> Full, ImageSize -> 200, Frame -> True, Axes -> True, 
			FrameTicks -> {{False, False}, {True, False}},
			PlotLabel -> "Estimation of phase", PlotHighlighting -> None];
		(*make the linear amplitude fit test plot*)
		pAmp = ListLinePlot[
			Transpose[{ppm, Rescale@Abs@ShiftedFourier[#]}] & /@ {dat, amps . Transpose[Transpose[fids] decay]}, 
			ScalingFunctions -> {"Reverse", None},
			ImageSize -> 200, Frame -> True, Axes -> True, 
			FrameTicks -> {{False, False}, {True, False}},
			PlotRange -> {Full, {0, 1.2}}, PlotLabel -> "Estimation of amplitudes", 
			PlotHighlighting -> None];
		
		(*make the plot row*)
		plots = Grid[{{pEps, pGam, pAmp, pPh}}, Spacings -> {2, 1}];
		
		(*output including plots*)
		{init, plots}
		,
		(*output without plots*)
		init
	]
]


(* ::Subsection:: *)
(*AmaresFitFidI*)


(* ::Subsubsection::Closed:: *)
(*AmaresFitFidI*)


AmaresFitFidI[dat_, time_, gyro_, basis_, matrix_?AssociationQ, cons_?MatrixQ] := AmaresFitFidI[dat, time, gyro, basis, matrix, {1, 1, 1}, {cons, 1}]

AmaresFitFidI[dat_, time_, gyro_, basis_, matrix_?AssociationQ, {cons_?MatrixQ, init_?MatrixQ}] := AmaresFitFidI[dat, time, gyro, basis, matrix, {1, 1, 1}, {cons, init}]

AmaresFitFidI[dat_, time_, gyro_, basis_, matrix_?AssociationQ, {tau_, tr_, fa_}, cons_?MatrixQ] := AmaresFitFidI[dat, time, gyro, basis, matrix, {tau, tr, fa}, {cons, 1}]

AmaresFitFidI[dat_, time_, gyro_, basis_, matrix_?AssociationQ, {tau_, tr_, fa_}, {cons_?MatrixQ, init_?MatrixQ}] := Block[{
		fids, names, nb, sc, datF, groups, parsV, parsM, con, int, parsFit,
		jacobian, residual, sigma, cov, se, t1Model, atp, conA
	},

	(*define fit parameters and constrians*)
	fids = Values[basis];
	names = Keys[basis];
	nb = Length@fids;

	(*make the t1model for b1 correction, if B1 is False wil be 1,0,.. and then tau tr and fa dont matter*)
	t1Model = MakeT1Model[names, matrix];
	t1Model = {
		GetSpinResonance /@ names,
		t1Values /@ names,
		Boole[matrix[#, "B1"] & /@ names]
	};

	(*determine data scaling for normalized data range*)
	sc = 10 / Max@Abs@dat; (*force all amplitudes between 0 and 10*)
	datF = ComplexToVector[dat sc];

	(*see how the parameters are grouped and make cons and inits accordingly*)
	groups = GetGroupIds[matrix, names];
	parsV = MakeParVec[groups];
	parsM = ParVecToMat[parsV, groups, 1.];(*always 6 in B1*)

	int = Thread[{parsV, ParMatToVec[init, groups]}];

	(*gamma atp always has to be biggest*)
	con = ParMatToVec[Transpose@ConstantArray[cons, nb], groups];
	atp = DeleteDuplicates[parsV[[groups[[1,Flatten[Position[names,#]&/@{"ATP-\[Gamma]","ATP-\[Alpha]","ATP-\[Beta]"}]]]]]];
	conA = If[Length[atp] > 1, #[[1]] >= #[[2]] & /@ Partition[atp, 2, 1], {}];
	con = Join[conA, Thread[con[[All, 1]] <= parsV <= con[[All, 2]]]];

	(*perform the minimization*)
	parsFit = parsV /. Last[Quiet@FindMinimum[
		{ObjectiveFunc[datF, fids, time, gyro, {fa, tr, tau}, t1Model, parsM], con}, int,
		Gradient -> GradientFunc[datF, fids, time, gyro, {fa, tr, tau}, t1Model, parsM, groups],
		MaxIterations -> 100, AccuracyGoal -> 10, PrecisionGoal -> 10]];
	parsFit = ParVecToMat[parsFit, groups, 1.];(*always 6 inc B1*)
	parsFit[[-2]] = Wrap[parsFit[[-2]]]; (*make phase as small as possible*)

	(*calculate the SE and noise sigma (CRLB)*)
	jacobian = JacobianFunc[datF, fids, time, gyro, {fa, tr, tau}, t1Model, parsFit, groups];
	residual = ResidualFunc[datF, fids, time, gyro, {fa, tr, tau}, t1Model, parsFit];
	sigma = Total[residual^2] / (Length[datF] - Length@parsV);
	cov = sigma PseudoInverse[jacobian . Transpose[jacobian]];
	se = ParVecToMat[Sqrt[Diagonal[cov]], groups, 0.];

	(*give the fit output but correct scaling*)
	parsFit[[1]] = parsFit[[1]] / sc;
	se[[1]] = se[[1]] / sc;
	{
		parsFit, 
		Round[se, .001], 
		Round[100 se/ReplacePart[Abs[parsFit], -2 -> 2 Pi], .1], 
		Round[Sqrt[sigma]/sc, .1]
	}
]


t1Values = <|
	"PE" -> 3900, "PC" -> 2300, "Piin" -> 4300, "Piex" -> 2300, "GPE" -> 4400, "GPC" -> 3900, "PTDC" -> 1000, 
	"PCr" -> 3800, "ATP-\[Gamma]" -> 560, "ATP-\[Alpha]" -> 570, "ATP-\[Beta]" -> 560, "NAD" -> 2000, "UDPG" -> 3300
|>;


GetSpinResonance[nn_?ListQ] := Association[# -> GetSpinResonance[#] & /@ nn]
GetSpinResonance[nn_?StringQ] := Block[{res, tags},
	{res, tags} = GetSpinSystem[nn /. Thread[{"ATP-\[Gamma]", "ATP-\[Alpha]", "ATP-\[Beta]"} -> "ATP"]][[{-4, -3}]];
	If[StringContainsQ[nn, "ATP-"],
		Mean[res[[First@Position[tags, Last[StringSplit[nn, "-"]]]]]],
		Mean[res]
	]
]


MakeT1Model[matrix_Association]:=MakeT1Model[Keys[matrix], matrix]

MakeT1Model[basis_Association, matrix_]:= MakeT1Model[Keys[basis], matrix]

MakeT1Model[names_List, matrix_]:={
	GetSpinResonance /@ names,
	t1Values /@ names,
	Boole[matrix[#, "B1"] & /@ names]
};


(* ::Subsubsection::Closed:: *)
(*AmaresSignalCalc*)


AmaresSignalCalc[basis_, matrix_, time_, pars_, gyro_?NumberQ]:=VectorToComplex[
	SignalCalcC[Values[basis], time, pars, MakeT1Model[basis, matrix], {1., 1. ,1.}, gyro]
]

AmaresSignalCalc[basis_, matrix_, time_, pars_, b1Par_?VectorQ, gyro_?NumberQ]:=VectorToComplex[
	SignalCalcC[Values[basis], time, pars, MakeT1Model[basis, matrix], b1Par, gyro]
]

SignalCalcC = Compile[{{bs, _Complex, 2}, {t, _Real, 1}, {pars, _Real, 2}, 
	{t1b1, _Real, 2}, {b1Par, _Real, 1}, {gyro, _Real}}, Block[{
		bRe, bIm, bt, as, sm, gm ,ep, ph, b1, epR, t1, b1Mask, fab1,
		fa, tr, tau, trt1, decay, theta, cosT, sinT, epTau, b1Amp,
		abRe, abIm
	},

	(*split the complex signal of the basis vectors*)
	bRe = Re[bs];
	bIm = Im[bs];

	(*make time vector to matrix*)
	bt = Table[t, {Length@bs}];

	(*Split parameter vector and b1 model parameters*)
	{as, sm, gm, ep, ph, b1} = pars;
	{epR, t1, b1Mask} = t1b1; (*reference ppm, t1 and b1 correction mask*)
	{fa, tr, tau} = b1Par;
	fa = fa Degree;
	fab1 = fa b1;
	tau = tau / 1000;
	trt1 = Exp[tr / t1];

	(*generate the parts of the signal equation.*)
	(*the line widths or the decay function*)
	decay = Exp[-Pi (gm bt + (sm bt)^2)];
	(*the phase rotation is complex rotation can be written as Sin Cos*)
	theta = -ph + 2 Pi gyro ep bt;
	cosT = decay Cos[theta];
	sinT = decay Sin[theta];
	(*the amplitudes with or without b1 correction*)
	epTau = (epR - ep) gyro Pi tau;
	b1Amp = b1Mask as ((-trt1 + Cos[fab1]) Csc[fab1] Sin[fab1 Sinc[epTau]]
		) / (-trt1 + Cos[fab1 Sinc[epTau]]) + as (1 - b1Mask);

	(*complex splitted basis signals with decay and phase*)
	abRe = b1Amp . (bRe cosT - bIm sinT);
	abIm = b1Amp . (bIm cosT + bRe sinT);

	(*The final signal*)
	Join[abRe, abIm]

], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*AmaresSignalJacobianCalc*)


AmaresSignalJacobianCalc[bs_, time_, pars_, t1b1_, b1Par_, gyro_]:=VectorToComplex@SignalJacCalcC[bs, time, pars, t1b1, b1Par, gyro]

SignalJacCalcC = Compile[{{bs, _Complex, 2}, {t, _Real, 1}, {pars, _Real, 2}, 
	{t1b1, _Real, 2}, {b1Par, _Real, 1}, {gyro, _Real, 0}}, Block[{
		bRe, bIm, bt, as, sm, gm ,ep, ph, b1, epR, t1, b1Mask, fa, tr, tau, trt1, decay, theta, 
		cosT, sinT, epTau, b1Amp, sigRe, sigIm, abRe, abIm, b1Dep, b1Db1, daRe, daIm, fab1,
		dsmRe, dsmIm, dgmRe, dgmIm, depRe, depIm, dphRe, dphIm, db1Re, db1Im, dRe, dIm, 
		sig, jac
	},
	(*split the complex signal of the basis vectors*)
	bRe = Re[bs];
	bIm = Im[bs];

	(*make time vector to matrix*)
	bt = Table[t, {Length@bs}];

	(*Split parameter vector and b1 model parameters*)
	{as, sm, gm, ep, ph, b1} = pars;
	{epR, t1, b1Mask} = t1b1; (*reference ppm, t1 and b1 correction mask*)
	{fa, tr, tau} = b1Par;
	fa = fa Degree;
	tau = tau / 1000;
	trt1 = Exp[tr / t1];

	(*generate the parts of the signal equation.*)
	(*the line widths or the decay function*)
	decay = Exp[-Pi (gm bt + (sm bt)^2)];
	(*the phase rotation is complex rotation can be written as Sin Cos*)
	theta = -ph + 2 Pi gyro ep bt;
	cosT = decay Cos[theta];
	sinT = decay Sin[theta];
	(*the amplitudes with or without b1 correction*)
	epTau = (epR - ep) gyro Pi tau;
	fab1 = fa b1;
	b1Amp = b1Mask as ((-trt1 + Cos[fab1]) Csc[fab1] Sin[fab1 Sinc[epTau]]
		) / (-trt1 + Cos[fab1 Sinc[epTau]]) + as (1 - b1Mask);

	(*complex splitted basis signals with decay and phase*)
	sigRe = bRe cosT - bIm sinT;
	sigIm = bIm cosT + bRe sinT;
	(*split signal modulated with amplitudes and b1*)
	abRe = b1Amp sigRe;
	abIm = b1Amp sigIm;

	(*start derivatives for jacobian*)
	(*derivatives of the b1 dependancy to ep and b1*)
	b1Dep = b1Mask (fab1 (trt1 - Cos[fab1]) (-1 + trt1 Cos[fab1 Sinc[epTau]]) Csc[fab1] (
		(ep - epR) gyro Pi tau Cos[epTau] + Sin[epTau])) / (
			(ep - epR)^2 gyro Pi tau (trt1 - Cos[fab1 Sinc[epTau]])^2);
	b1Db1 = b1Mask (fa Csc[fab1] ((-1 + trt1 Cos[fab1]) (-trt1 + Cos[fab1 Sinc[epTau]]
		) Csc[fab1] Sin[fab1 Sinc[epTau]] + (trt1 - Cos[fab1]) (-1 + trt1 Cos[fab1 Sinc[epTau]]) Sinc[epTau])
			) / (trt1 - Cos[fab1 Sinc[epTau]])^2;
	(* derivatives to amplitudes*)
	daRe = sigRe;
	daIm = sigIm;
	(*real derivatives parameters for gausiang and laurentzian decay parameters*)
	dsmRe = -2 Pi sm bt^2 abRe ;
	dsmIm = -2 Pi sm bt^2 abIm ;
	dgmRe = -Pi bt abRe;
	dgmIm = -Pi bt abIm;
	(*complex derivatives for the shift including b1part*)
	depRe = -2 Pi gyro bt abIm + as b1Dep daRe;
	depIm = 2 Pi gyro bt abRe + as b1Dep daIm;
	(*complex derivatives for the phase*)
	dphRe = abIm;
	dphIm = -abRe;
	(*derivative to b1*)
	db1Re = as b1Db1 daRe;
	db1Im = as b1Db1 daIm;

	(*combining all derivatives into the full jacobian*)
	dRe = Join[daRe, dsmRe, dgmRe, depRe, dphRe, db1Re];
	dIm = Join[daIm, dsmIm, dgmIm, depIm, dphIm, db1Im];

	(*output the signal and the jacobian*)
	sig = {Total[Join[abRe, abIm, 2]]};
	jac = Join[dRe, dIm, 2];

	Join[sig, jac]

], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*GetGroupIds*)


GetGroupIds[matrix_] := GetGroupIds[matrix] = GetGroupIds[matrix, Keys[matrix]]

GetGroupIds[matrix_, names_] := Table[
	If[lab === "B1",
		If[AnyTrue[Values[matrix[[All, lab]]], # === True &], 
			ConstantArray[1, Length@names], 
			Nothing
		],
		ResolveGroupIds[matrix[#, lab] & /@ names]
	]
, {lab, Most[fitMatrixParams]}]


ResolveGroupIds[codes_List] := Module[{result, nextId, pos},
	result = ConstantArray[0, Length[codes]];
	nextId = 1;
	
	(*step 1: all default-group (code 0) entries share one id*)
	pos = Flatten[Position[codes, 0]];
	If[Length[pos] > 0, result[[pos]] = nextId; nextId += 1];
	
	(*step 2: each distinct explicit positive code becomes its own group, in order of first appearance*)
	Do[
		pos = Flatten[Position[codes, g]];
		result[[pos]] = nextId; nextId += 1;
	, {g, DeleteDuplicates[Select[codes, # > 0 &]]}];
	
	(*step 3: each unique (-1) entry gets its own group*)
	pos = Flatten[Position[codes, -1]];
	Do[result[[p]] = nextId; nextId += 1, {p, pos}];

	(*give the resolved group id*)
	result
];


(* ::Subsubsection::Closed:: *)
(*MakeParVec*)


MakeParVec[groupsId_] := Flatten[Table[Unique[#[[1]]], {#[[2]]}] & /@ Thread[
	{Most[fitMatrixParams][[;; Length@groupsId]], Max /@ groupsId}]]


(* ::Subsubsection::Closed:: *)
(*ParVecToMat*)


ParVecToMat[pvec_, groupsId_] := #[[1]][[#[[2]]]] & /@ Thread[{TakeList[pvec, Max /@ groupsId], groupsId}]


ParVecToMat[pvec_, groupsId_, fixedB1_] := Block[{expanded, nRows},
	expanded = ParVecToMat[pvec, groupsId];(*however many rows groupsId actually has*)
	nRows = Length[groupsId];
	(*full matrix always contains 6*)
	If[nRows < 6, Append[expanded, ConstantArray[fixedB1, Length[First[expanded]]]], expanded]
];


(* ::Subsubsection::Closed:: *)
(*ParMatToVec*)


ParMatToVec[mat_, groupsId_] := Block[{np, pp},
	{np, pp} = Dimensions[groupsId];
	Flatten[Table[(
		Mean[mat[[i, Flatten[Position[groupsId[[i]], #]]]]]) & /@ Range[Max[groupsId[[i]]]]
	, {i, 1, np}], 1]
]


(* ::Subsubsection::Closed:: *)
(*ObjectiveFunc*)


ObjectiveFunc[data_, fids_, time_, gyro_, b1par_, b1Model_, p_?(MatrixQ[#, NumericQ] &)] := 
	Total[(ResidualFunc[data, fids, time, gyro, b1par, b1Model, p])^2]


(* ::Subsubsection::Closed:: *)
(*ResidualFunc*)


ResidualFunc[data_, fids_, time_, gyro_, b1par_, b1Model_, p_?(MatrixQ[#, NumericQ] &)] := 
	data - SignalCalcC[fids, time, p, b1Model, b1par, gyro]


(* ::Subsubsection::Closed:: *)
(*JacobianFunc*)


JacobianFunc[data_, fids_, time_, gyro_, b1par_, b1Model_, p_?(MatrixQ[#, NumericQ] &), grp_] := 
	CollapseJacMat[Rest@SignalJacCalcC[fids, time, p, b1Model, b1par, gyro], grp]


(* ::Subsubsection::Closed:: *)
(*GradientFunc*)


GradientFunc[data_, fids_, time_, gyro_, b1par_, b1Model_, p_?(MatrixQ[#, NumericQ] &), grp_] := Block[{sigjac},
	sigjac = SignalJacCalcC[fids, time, p, b1Model, b1par, gyro];
	-2 CollapseJacMat[Rest@sigjac, grp] . (data - First@sigjac)
]


CollapseJacMat[jac_, groupsId_] := Block[{np, pp, jacMat},
	{np, pp} = Dimensions[groupsId];
	jacMat = Partition[jac, pp];
	Flatten[Table[
		(Total[jacMat[[i, Flatten[Position[groupsId[[i]], #]]]]]) & /@ Range[Max[groupsId[[i]]]]
	, {i, 1, np}], 1]
]


(* ::Subsection::Closed:: *)
(*FitB1FromATP*)


FitB1FromATP[fit_, matrix_, {tau_, fa_, tr_, gyro_}] := FitB1FromATP[fit, matrix, {tau, fa, tr, gyro}, False]

FitB1FromATP[fit_, matrix_, {tau_, fa_, tr_, gyro_}, plot_] := Block[{
		tp, pos, sig, amp, b1f, b1, par, cor, b1fit, res, t1s, p1, p2, pl, atp, met
	},
	met = Keys[matrix];
	{res, t1s} = MakeT1Model[met, matrix][[1 ;; 2]];
	atp = {"ATP-\[Gamma]", "ATP-\[Alpha]", "ATP-\[Beta]"};
	pos = Flatten[Position[met, #] & /@ atp];
	sig = Thread[{
		fit[[1, pos]],
		fit[[-3, pos]] + res[[pos]],
		t1s[[pos]]}
		];

	b1fit = 
	Last@FindMinimum[{errorB1[sig, {tau, fa, tr, gyro}, {amp, b1f}], 
		0.5 < b1f < 2.5}, {{amp, sig[[1, 1]]}, {b1f, 1}}, MaxIterations -> 100];

	b1 = b1f /. b1fit;
	par = Thread[{met, fit[[-2]] + res, t1s}];
	cor = 
	Association[#[[1]] -> 
		CalcSaturationProfile[tau, {#[[2]], gyro}, {fa, {tr, #[[3]]}},
			b1] & /@ par];

	If[plot,
	p1 = ListPlot[Thread[{
		Join[sig[[All, 2]], -sig[[All, 2]]],
		Join[sig[[All, 1]], sig[[All, 1]]]}],
		PlotRange -> {{-25, 25}, {0, 1.2 Max@sig[[All, 1]]}}, 
		PlotHighlighting -> None, ScalingFunctions -> {"Reverse", None}, Ticks -> None, 
		PlotStyle -> Directive[Red, PointSize[Large]]
		];
	p2 = 
		Plot[Evaluate[(amp CalcSaturationProfile[
				tau, {pp, gyro}, {fa, {tr, #}}, b1f] & /@ 
			sig[[All, 3]]) /. b1fit], {pp, -25, 25},
		ScalingFunctions -> {"Reverse", None}, PlotHighlighting -> None];
	pl = Show[p1, p2, PlotLabel -> Round[b1f /. b1fit, .01]];
	{b1, cor, pl},
	{b1, cor}
	]
	];

errorB1[sig_, {tau_, fa_, tr_, gyro_}, {a_?NumericQ, b1_?NumericQ}] := Block[{amp, eps, t1, err},
	err = 
		({amp, eps, t1} = #; amp - a CalcSaturationProfile[tau, {eps, gyro}, {fa, {tr, t1}}, b1]
	) & /@ sig;
	Total[err^2]
]

CalcSaturationProfile[tau_, {eps_, gyro_}, {fa_, {tr_, t1_}}, b1_] := Block[{e1, s0, faAct, loc, sLoc},
	faAct = b1*fa;
	e1 = Exp[-tr/t1];
	loc = Sinc[Pi (eps gyro) (tau/1000.)] ;
	sLoc = ((1 - e1) Sin[loc faAct Degree]/(1 - e1 Cos[loc faAct Degree]));
	s0 = (1 - e1) Sin[faAct Degree]/(1 - e1 Cos[faAct Degree]);
	sLoc/s0
];

SaturationProfileCorrection[fa_, tr_, t1_?MatrixQ, b1_] := SaturationProfileCorrection[fa, tr, t1[[2]], b1]
SaturationProfileCorrection[fa_, tr_, t1_?VectorQ, b1_] := Block[{e1, s0Act, s0, faAct, fa2, s},
	If[b1 === 0., ConstantArray[0., Length@t1],
		faAct = b1*fa;
		e1 = Exp[-tr/t1];
		s0 = (1 - e1) Sin[fa Degree]/(1 - e1 Cos[fa Degree]);
		s0Act = (1 - e1) Sin[faAct Degree]/(1 - e1 Cos[faAct Degree]);
		s0/s0Act
	]
];




(* ::Section:: *)
(*End Package*)


End[](* End Private Context *)

EndPackage[]
