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


EditAmaresMatrix::usage = 
"EditAmaresMatrix[matrix] is used to edit the Amares fitting matrix which can be compressed or extended using
ToFullAmaresMatrix or FromFullAmaresMatrix. 
The default properties are \"amp\": -1, \"sigma\": 0, \"gamma\": 0, \"shift\": -1, \"phase\": 0, \"B1\": False, \"init\": False.

EditAmaresMatrix[matrix, names] will do the same but will only do it for the values in the list names. If values are not in matrix
they will be added, if values are names but not in matrix they will be removed from matrix."

ToFullAmaresMatrix::usage = 
"ToFullAmaresMatrix[names_List] completes the Amares fitting matrix to its full form by using defaults for the values not defined.
To edit the matrix the function EditAmaresMatrix can be used and to compress it FromFullAmaresMatrix can be used."

FromFullAmaresMatrix::usage = 
"FromFullAmaresMatrix[matrix] compresses the Amares fitting matrix to its full form by removing all defaults from the defined.
To edit the matrix the function EditAmaresMatrix can be used and to extend it ToFullAmaresMatrix can be used."

AmaresBasis::usage = 
"AmaresBasis[settings, {metabolites, split}] generates basis function needed for Amares fitting.
The settings are defined as {nSamp, bw, field, nuc, te, f0} and the metabolites are a member of
\"PE\", \"PC\", \"Piin\", \"Piex\", \"GPE\", \"GPC\", \"PTDC\", \"PCr\", \"ATP\", \"NAD\" and \"UDPG\".
When \"ATP\" is defined in split the 3 ATP peaks will be separate spectra to fit. 
The output is {basisFunc, {time, ppm, gyro}} where basisFunc is an Association."

AmaresInitialValues::usage = 
"AmaresInitialValues[fid, time, basis, matrix, {gyro, dw}] makes a guess for the initial fit values to initialize AmaresFitFidI.
The time, basis and gyro can be obtained from AmaresBasis and the matrix can be obtained form EditMaresMatrix."

AmaresFitFidI::usage = 
"AmaresFitFidI[fid, time, gyro, basis, matrix, cons] fits the fid with using the basis spectra.
AmaresFitFidI[fid, time, gyro, basis, matrix, {cons, init}] fits the fid using initial values init which can be obtained from a 
previous fit or from AmaresInitialValues.
AmaresFitFidI[fid, time, gyro, basis, matrix, b1pars, cons] the b1pars are needed if in the matrix B1 correction is set True.
The b1pars are defined as {tau, fa, TR}  with tau and TR in ms and fa in degrees.
AmaresFitFidI[fid, time, gyro, basis, matrix, b1pars, {cons, init}] fits using initial values using b1pars for B1 correction.

The time, basis and gyro can be obtained from AmaresBasis and the matrix can be obtained form EditMaresMatrix.
The cons have to be defined for {\"amp\", \"sigma\", \"gamma\", \"shift\", \"phase\", \"B1\"} as {min, max} for each.
Output is the {fitted parameters, the SE, the SE in %, the standard deviation of the residuals.}."

ShowParameterMatrix::usage = 
"ShowParameterMatrix[fit, metabolites] shows the fitted or initial fit parameter matrix from AmaresFitFidI and AmaresInitialValues using met as labels.
ShowParameterMatrix[fit, true, metabolites] shows the relative difference between the true and fitted parameters."

AmaresSignalCalc::usage = 
"AmaresSignalCalc[basis, matrix, pars, gyro] generates a fid from the basis functions according to matrix and pars.
AmaresSignalCalc[basis, matrix, pars, b1Pars, gyro] also generates a fid but also using the B1 effect if set True in matrix.
The values of b1Pars are {tau, fa, TR}  with tau and TR in ms and fa in degrees."

AmaresSignalJacobianCalc::usage = 
"AmaresSignalJacobianCalc[basis, matrix, pars, gyro] generates a fid and the Jacobian from the basis functions according to matrix and pars.
AmaresSignalJacobianCalc[basis, matrix, pars, b1Pars, gyro] also generates a fid but also using the B1 effect if set True in matrix.
The values of b1Pars are {tau, fa, TR}  with tau and TR in ms and fa in degrees."


$AmaresB1Regularization::usage = 
"$AmaresB1Regularization is de regularization factor used for b1 regularization during Amares fitting."

$AmaresT1Values::usage = 
"$AmaresT1Values is an association that holds the T1 values of the metabolites used in Amares fitting."

$AmaresMatrixDefaults::usage = 
"$AmaresMatrixDefaults is an association that holds the default values for the Amares fit matrix."

(* ::Subsection:: *)
(*Options*)


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


$AmaresT1Values = <|
	"PE" -> 3900, "PC" -> 2300, "Piin" -> 4300, "Piex" -> 2300, "GPE" -> 4400, "GPC" -> 3900, "PTDC" -> 1000, 
	"PCr" -> 3800, "ATP-\[Gamma]" -> 560, "ATP-\[Alpha]" -> 570, "ATP-\[Beta]" -> 560, "NAD" -> 2000, "UDPG" -> 3300
|>;


$AmaresB1Regularization = 0.000;


$AmaresMatrixDefaults = <|"amp" -> -1, "sigma" -> 0, "gamma" -> -1, "shift" -> -1, "phase" -> 0, "B1" -> False, "init" -> False|>;


ComplexToVector = Join @@ Through[{Re, Im}[#]] &;


VectorToComplex = {1, I} . Partition[#, Length[#]/2] &;


(* ::Subsection::Closed:: *)
(*ShowParameterMatrix*)


ShowParameterMatrix[fit_, mets_] := Block[{fiti=fit},
	fiti[[1]] = Round[100 fiti[[1]] / fiti[[1, 8]], .1];
	ShowMat[fiti, mets]
]


ShowParameterMatrix[fit_, ref_, mets_] := Block[{fiti, refi},
	fiti= fit; refi = ref;
	fiti[[1]] = Round[100 fit[[1]] / fit[[1, 8]], .1];
	refi[[1]] = Round[100 ref[[1]] / ref[[1, 8]], .1];
	fiti = fiti - refi;
	ShowMat[fiti, mets]
]


ShowMat[fiti_, mets_]:=MatrixForm[
	Transpose@Prepend[Transpose[Prepend[Round[fiti, .01], mets]], 
	{"met", "amp", "sigma", "gamma", "eps", "phi", "B1"}[[;; Length[fiti] + 1]]]
]


(* ::Subsection::Closed:: *)
(*AmaresBasis*)


AmaresBasis[settings_, metabolites_?VectorQ] := AmaresBasis[settings, {metabolites, {}}]

AmaresBasis[settings_, {metabolites_?ListQ, split_?ListQ}] := Block[{
		nSamp, bw, field, nuc, te, f0, dw, gyro, time, ppm, names, fids, basisFunc, sp, spin
	},

	{nSamp, bw, field, nuc, te, f0} = settings;
	dw = 1./bw;
	gyro = GetGyro[field, nuc];

	{time, ppm} = GetTimePpmRange[nSamp, dw, gyro] + {te, 0};
	{names, fids, sp, spin} = GetSpectraBasisFunctions[metabolites, split,
		BasisSequence -> {"PulseAcquire", te},
		SpectraSamples -> nSamp, SpectraBandwidth -> bw,
		SpectraPpmShift -> f0, SpectraFieldStrength -> field,
		SpectraNucleus -> nuc];
	basisFunc = Association[Thread[names -> fids]];

	{basisFunc, {time, ppm, gyro}}
]


(* ::Subsection::Closed:: *)
(*AmaresFitMatrix*)


(* ::Subsubsection::Closed:: *)
(*ToFullAmaresMatrix*)


ToFullAmaresMatrix[names_List] := Association[# -> $AmaresMatrixDefaults & /@ names];

ToFullAmaresMatrix[matrix_Association] := Association[# -> Join[$AmaresMatrixDefaults, matrix[#]] & /@ Keys[matrix]];

ToFullAmaresMatrix[matrix_Association, names_List] := Association[# -> Join[$AmaresMatrixDefaults, Lookup[matrix, #, <||>]] & /@ names];


(* ::Subsubsection::Closed:: *)
(*FromFullAmaresMatrix*)


FromFullAmaresMatrix[matrix_Association] := Select[
	Association[# -> Association[
		KeyValueMap[Function[{k, v}, If[v === $AmaresMatrixDefaults[k], Nothing, k -> v]], matrix[#]]
		] & /@ Keys[matrix]],
	# =!= <||> &]


(* ::Subsubsection::Closed:: *)
(*EditAmaresMatrix*)


EditAmaresMatrix[names_List] := EditAmaresMatrix[<||>, names];

EditAmaresMatrix[matrix_Association] := EditAmaresMatrix[matrix, Keys[matrix]];

EditAmaresMatrix[matrix_Association, names_List] := Block[{start, return, list, nn, np, cellVars, cellStart, fitPar},
	fitPar = Keys[$AmaresMatrixDefaults];
	list = Join[{-1 -> "unique", 0 -> "same"}, # -> "group-" <> IntegerString[#, 10, 2] & /@ Range[1, Length[names]]];
	start = ToFullAmaresMatrix[matrix, names];
	{nn, np} = {Length[names], Length[fitPar]};
	cellStart = Table[start[names[[i]], fitPar[[j]]], {i, nn}, {j, np}];

	return = DialogInput[DynamicModule[{data = cellStart, datai = cellStart},
		Column[{
			Grid[Prepend[Table[Prepend[Table[
				With[{i = i, j = j, param = fitPar[[j]]},
					If[param === "init" || param === "B1",
						Checkbox[Dynamic[data[[i, j]]]],
						PopupMenu[Dynamic[data[[i, j]]], list]
					]
				]
				, {j, np}], names[[i]]], {i, nn}], Prepend[fitPar, "Metabolite"]]
			, Frame -> All, Alignment -> Center, Spacings -> {1, 1}]
			,
			Row[{
				Button["Save", DialogReturn[data], ImageSize -> 80],
				Spacer[10],
				Button["Cancel", DialogReturn[datai], ImageSize -> 80]
			}]
		}]
	]];
	Association[Thread[names -> (Association[Thread[fitPar -> #]] & /@ return)]]
]


(* ::Subsection:: *)
(*Amares fit groups*)


(* ::Subsubsection::Closed:: *)
(*GetGroupIds*)


GetGroupIds[matrix_] := GetGroupIds[matrix] = GetGroupIds[matrix, Keys[matrix]]

GetGroupIds[matrix_, names_] := Table[
	(*only make B1 group IDs if any B1 fitting is true*)
	If[lab === "B1",
		If[AnyTrue[Values[matrix[[All, lab]]], # === True &], 
			ConstantArray[1, Length@names], 
			Nothing
		],
		ResolveGroupIds[matrix[#, lab] & /@ names]
	]
, {lab, Most[Keys[$AmaresMatrixDefaults]]}]


ResolveGroupIds[codes_List] := Module[{result, nextId, pos},
	(*step 0: initialize group matrix and id*)
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
	{Most[Keys[$AmaresMatrixDefaults]][[;; Length@groupsId]], Max /@ groupsId}]]


(* ::Subsubsection::Closed:: *)
(*ParVecToMat*)


ParVecToMat[pVec_, groupsId_] := #[[1]][[#[[2]]]] & /@ Thread[{TakeList[pVec, Max /@ groupsId], groupsId}]


ParVecToMat[pVec_, groupsId_, fixedB1_] := Block[{expanded, nRows},
	expanded = ParVecToMat[pVec, groupsId];(*however many rows groupsId actually has*)
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


(* ::Subsection::Closed:: *)
(*MakeT1Model*)


(* ::Subsubsection::Closed:: *)
(*MakeT1Model*)


MakeT1Model[matrix_Association]:=MakeT1Model[Keys[matrix], matrix]

MakeT1Model[basis_Association, matrix_]:= MakeT1Model[Keys[basis], matrix]

MakeT1Model[names_List, matrix_]:={
	GetSpinResonance /@ names,
	$AmaresT1Values /@ names,
	Boole[matrix[#, "B1"] & /@ names]
};


(* ::Subsubsection::Closed:: *)
(*GetSpinResonance*)


GetSpinResonance[nn_?ListQ] := Association[# -> GetSpinResonance[#] & /@ nn]
GetSpinResonance[nn_?StringQ] := Block[{res, tags},
	{res, tags} = GetSpinSystem[nn /. Thread[{"ATP-\[Gamma]", "ATP-\[Alpha]", "ATP-\[Beta]"} -> "ATP"]][[{-4, -3}]];
	If[StringContainsQ[nn, "ATP-"],
		Mean[res[[First@Position[tags, Last[StringSplit[nn, "-"]]]]]],
		Mean[res]
	]
]


(* ::Subsection::Closed:: *)
(*AmaresInitialValues*)


AmaresInitialValues[dat_, time_, basis_, matrix_, {gyro_, dw_}] := AmaresInitialValues[dat, time, basis, matrix, {gyro, dw}, False]

AmaresInitialValues[dat_, time_, basis_, matrix_, {gyro_, dw_}, pl_] := Block[{
		sc, datf, model, rf, a, p, e, g, s, am, gm, sm, eps, ph, epi, ami, gmi,
		fids, base, sig, corr, ppm, wght, fit, datA, obj,
		decay, amp, opts, pEps, pGam, pPh, pAmp, plots, init
	},
	
	(*Normalize the signal*)
	sc = 1. / Max[Abs@dat];
	datf = sc dat;
	
	(*the signal model for linewidth*)
	model[a_?NumberQ, g_] := a Exp[-Pi (g time + 2 Pi g^2 time^2)];
	model[a_?NumberQ, p_, e_, g_, s_] := a Exp[-Pi (g time + 2 Pi s^2 time^2) + (-p + 2 Pi e gyro time) I];
	model[a_?VectorQ, p_, e_, g_, s_] := model[#, p, e, g, s] & /@ a;
	rf = Rescale[Abs[ShiftedFourier[#]]]&;

	(*get basis functions and base for calibration*)
	fids = Values[basis];
	base = Total[Pick[fids, Table[matrix[name, "init"], {name, Keys[basis]}]]];

	(*Step 1. find initial shift using autocorrelation*)
	corr = ListCorrelate[rf[base], rf[datf], Round[Length[base]/2], 0];
	(*find max peak closest to 0 after smoothing auto correlation*)
	ppm = GetPpmRange[datf, dw, gyro];
	sig = 5 (*ppm*);
	wght = 1 / (Exp[(ppm^2 / (2 sig^2))] (Sqrt[2 Pi] sig));
	epi = ppm[[Position[wght corr, Max[wght corr]][[1, 1]]]];
	
	(*Step 2. find initial linewidth using the absolute signal*)
	obj[a_?NumberQ, g_?NumericQ] := Total[(Abs[datf] - model[a, g])^2];
	fit = Last[Quiet@FindMinimum[{obj[a, g], g > 1}, 
		{{a, 1}, {g, 10}}, MaxIterations -> 50]];
	{ami, gmi} = {a, g} /. fit;

	(*Step 3. guess the amplitudes*)
	am = Abs[LeastSquares[Transpose[fids] model[1, 0, epi, gmi, gmi], datf]];

	(*Step 4. find the phase and fine-tune the shift, lw*)
	obj[a_?VectorQ, p_?NumericQ, e_?NumericQ, g_?NumericQ, s_?NumericQ] := Total[Abs[datf - Total[fids model[a, p, e, g, s]]]^2];
	fit = Last[Quiet@FindMinimum[{obj[am, p, epi, gmi, gmi], 
			-2 Pi < p < 2 Pi, Max[{1, gmi - 5}] <= g <= gmi + 5, epi - 0.2 < e < epi + 0.2}, 
		{{p, 0}, {e, epi}, {g, gmi}}], MaxIterations -> 50];
	{ph, eps, gm, sm} = {p, e, g, g} /. fit;

	(*Step 5. fine-tune the amplitudes*)
	decay = model[1, ph, eps, gm, sm];
	amp = Abs[LeastSquares[Transpose[fids] decay, datf]];

	(*Finalize by generating parameter matrix*)
	init = Join[{amp / sc}, ConstantArray[#, Length@amp] & /@ {sm, gm, eps, ph, 1}];

	(*Make the plots if needed*)
	If[pl,
		opts = Sequence[ImageSize -> 200, Frame -> True, Axes -> True, PlotHighlighting -> None,
			FrameTicks -> {{False, False}, {True, False}}];
		(*Step 1. make autocorrelation test plot*)
		pEps = ListLinePlot[
			Transpose[{ppm, #}] & /@ {rf[datf], rf[base model[1, 0, epi, 0, 0]], Rescale[wght corr]}, 
			ScalingFunctions -> {"Reverse", None}, PlotRange -> {Full, {0, 1.2}}, PlotLabel -> "1. Estimation of shift", opts];
		(*Step 2. make the Linewidth test plot*)
		pGam = ListLinePlot[
			Transpose[{1000 time, Abs[#]}] & /@ {datf, model[ami, gmi]},
			PlotRange -> Full, PlotLabel -> "2. Estimation of linewidth", opts];
		(*Step 3. make the phase test plot*)
		pPh = ListLinePlot[
			Transpose[{1000 time, Re[#]}] & /@ {datf, decay (am . fids)},
			PlotRange -> Full, PlotLabel -> "3. Estimation of phase", opts];
		(*Step 4. make the linear amplitude fit test plot*)
		pAmp = ListLinePlot[
			Transpose[{ppm, rf[#]}] & /@ {datf, decay (amp . fids)}, 
			ScalingFunctions -> {"Reverse", None}, PlotRange -> {Full, {0, 1.2}}, PlotLabel -> "4. Fine-tune estimations", opts];
		
		(*output including plots*)
		{init, Grid[{{pEps, pGam, pPh, pAmp}}, Spacings -> {2, 1}]}
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

AmaresFitFidI[dat_, time_, gyro_, basis_, matrix_?AssociationQ, {tau_, fa_, tr_}, cons_?MatrixQ] := AmaresFitFidI[dat, time, gyro, basis, matrix, {tau, fa, tr}, {cons, 1}]

AmaresFitFidI[dat_, time_, gyro_, basis_, matrix_?AssociationQ, {tau_, fa_, tr_}, {cons_?MatrixQ, init_?MatrixQ}] := Block[{
		fids, names, nb, sc, datF, groups, parsV, parsM, con, int, parsFit,
		jacobian, residual, sigma, cov, se, t1Model, atp, conA
	},

	(*define fit parameters and constrains*)
	fids = Values[basis];
	names = Keys[basis];
	nb = Length@fids;

	(*make the t1model for b1 correction, if B1 is False wil be 1,0,.. and then tau tr and fa does not matter*)
	t1Model = MakeT1Model[names, matrix];

	(*determine data scaling for normalized data range*)
	sc = 1. / Max[Abs@dat]; (*force all amplitudes between 0 and 1*)
	datF = ComplexToVector[sc dat];

	(*see how the parameters are grouped and make cons and inits accordingly*)
	groups = GetGroupIds[matrix, names];
	parsV = MakeParVec[groups];
	parsM = ParVecToMat[parsV, groups, 1.];(*always 6 inc B1*)

	(*create the initial fit values vector with scaled amplitudes*)
	int = init;
	int[[1]] = sc int[[1]];
	int = Thread[{parsV, ParMatToVec[int, groups]}];

	(*create the fit constraints*)
	con = ParMatToVec[Transpose@ConstantArray[cons, nb], groups];
	con = Thread[con[[All, 1]] <= parsV <= con[[All, 2]]];

	(*perform the minimization*)
	parsFit = parsV /. Last[Quiet@FindMinimum[
		{ObjectiveFunc[datF, fids, time, gyro, {tau, fa, tr}, t1Model, parsM], con}, int,
		Gradient -> GradientFunc[datF, fids, time, gyro, {tau, fa, tr}, t1Model, parsM, groups]],
		MaxIterations -> 100];
	parsFit = ParVecToMat[parsFit, groups, 1.];(*always 6 inc B1*)
	parsFit[[-2]] = Wrap[parsFit[[-2]]]; (*wrap the phase to [-Pi, PI]*)

	(*calculate the SE and noise sigma (CRLB)*)
	residual = ResidualFunc[datF, fids, time, gyro, {tau, fa, tr}, t1Model, parsFit];
	sigma = Total[residual^2] / (Length[datF] - Length@parsV);
	jacobian = JacobianFunc[datF, fids, time, gyro, {tau, fa, tr}, t1Model, parsFit, groups];
	cov = sigma PseudoInverse[jacobian . Transpose[jacobian]];
	se = ParVecToMat[Sqrt[Diagonal[cov]], groups, 0.];

	(*give the fit output with correct scaling*)
	parsFit[[1]] = parsFit[[1]] / sc;
	se[[1]] = se[[1]] / sc;
	{
		parsFit, 
		Round[se, .001], 
		Round[100 se / ReplacePart[Abs[parsFit], -2 -> 2 Pi], .1], 
		Round[Sqrt[sigma] / sc, .1]
	}
]


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
		bRe, bIm, bt, as, sm, gm ,ep, ph, b1, epR, t1, b1Mask, fab1, fSyncE,
		fa, tr, tau, trt1, decay, theta, cosT, sinT, epTau, b1Amp, s0, s0Loc,
		abRe, abIm, mAmp
	},

	(*split the complex signal of the basis vectors*)
	bRe = Re[bs];
	bIm = Im[bs];

	(*make time vector to matrix*)
	bt = Table[t, {Length@bs}];

	(*Split parameter vector and b1 model parameters*)
	{as, sm, gm, ep, ph, b1} = pars;
	{epR, t1, b1Mask} = t1b1; (*reference ppm, t1 and b1 correction mask*)
	{tau, fa, tr} = b1Par;
	fa = fa Degree;
	tau = tau / 1000;

	(*generate the parts of the signal equation.*)
	
	(*the line widths or the decay function*)
	decay = Exp[-Pi (gm bt + 2 Pi sm^2 bt^2)];
	(*the phase rotation is complex rotation can be written as Sin[] Cos[]*)
	theta = -ph + 2 Pi gyro ep bt;
	cosT = decay Cos[theta];
	sinT = decay Sin[theta];

	(*The ernst angle b1 model*)
	epTau = (epR - ep) gyro Pi tau;
	trt1 = Exp[-tr / t1];
	fab1 = fa b1;
	fSyncE = fab1 Sinc[epTau];
	s0 = (1 - trt1) Sin[fa] / (1 - trt1 Cos[fa]);(*amp normalization*)
	s0Loc = (1 - trt1) Sin[fSyncE] / (1 - trt1 Cos[fSyncE]); (*t1 b1 correction at each location*)
	
	(*the amplitudes with or without b1 correction*)
	mAmp = (1 - b1Mask) + b1Mask (s0Loc / s0);
	b1Amp = as mAmp;

	(*complex split basis signals with decay and phase*)
	abRe = b1Amp . (bRe cosT - bIm sinT);
	abIm = b1Amp . (bIm cosT + bRe sinT);

	(*The final signal*)
	Join[abRe, abIm]

], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*AmaresSignalJacobianCalc*)


AmaresSignalJacobianCalc[basis_, matrix_, time_, pars_, gyro_?NumberQ]:= SignalJacCalcC[
	Values[basis], time, pars, MakeT1Model[basis, matrix], {1., 1. ,1.}, gyro]

AmaresSignalJacobianCalc[basis_, matrix_, time_, pars_, b1Par_?VectorQ, gyro_?NumberQ]:= SignalJacCalcC[
	Values[basis], time, pars, MakeT1Model[basis, matrix], b1Par, gyro]


SignalJacCalcC = Compile[{{bs, _Complex, 2}, {t, _Real, 1}, {pars, _Real, 2}, 
	{t1b1, _Real, 2}, {b1Par, _Real, 1}, {gyro, _Real, 0}}, Block[{
		bRe, bIm, bt, as, sm, gm ,ep, ph, b1, epR, t1, b1Mask, fa, tr, tau, trt1, decay, theta, s0, s0Loc,
		cosT, sinT, epTau, b1Amp, mAmp, sigRe, sigIm, abRe, abIm, b1Dep, b1Db1, daRe, daIm, fab1,
		dsmRe, dsmIm, dgmRe, dgmIm, depRe, depIm, dphRe, dphIm, db1Re, db1Im, dRe, dIm, fSyncE,
		sig, jac, dsLocDep, dsLocDb1
	},
	(*split the complex signal of the basis vectors*)
	bRe = Re[bs];
	bIm = Im[bs];

	(*make time vector to matrix*)
	bt = Table[t, {Length@bs}];

	(*Split parameter vector and b1 model parameters*)
	{as, sm, gm, ep, ph, b1} = pars;
	{epR, t1, b1Mask} = t1b1; (*reference ppm, t1 and b1 correction mask*)
	{tau, fa, tr} = b1Par;
	fa = fa Degree;
	tau = tau / 1000;

	(*generate the parts of the signal equation.*)
	
	(*the line widths or the decay function*)
	decay = Exp[-Pi (gm bt + 2 Pi sm^2 bt^2)];
	(*the phase rotation is complex rotation can be written as Sin Cos*)
	theta = -ph + 2 Pi gyro ep bt;
	cosT = decay Cos[theta];
	sinT = decay Sin[theta];
	
	(*The ernst angle b1 model*)
	epTau = (epR - ep) gyro Pi tau;
	trt1 = Exp[-tr / t1];
	fab1 = fa b1;
	fSyncE = fab1 Sinc[epTau];
	s0 = (1 - trt1) Sin[fa] / (1 - trt1 Cos[fa]);(*amp normalization*)
	s0Loc = (1 - trt1) Sin[fSyncE] / (1 - trt1 Cos[fSyncE]); (*t1 b1 correction at each location*)

	(*the amplitudes with or without b1 correction*)
	mAmp = (1 - b1Mask) + b1Mask (s0Loc / s0);
	b1Amp = as mAmp;

	(*complex split basis signals with decay and phase*)
	sigRe = bRe cosT - bIm sinT;
	sigIm = bIm cosT + bRe sinT;
	(*split signal modulated with amplitudes and b1*)
	abRe = b1Amp sigRe;
	abIm = b1Amp sigIm;

	(*start derivatives for jacobian*)
	
	(*derivatives of the b1 dependency to ep and b1*)
	dsLocDep = (fab1 (-1 + trt1) (trt1 - Cos[fSyncE]) (
			(ep - epR) gyro Pi tau Cos[epTau] + Sin[epTau])
				) / ((ep - epR)^2 gyro Pi tau (-1 + trt1 Cos[fSyncE])^2);
	b1Dep = as b1Mask (dsLocDep / s0);
	dsLocDb1 = (fa (-1 + trt1) (trt1 - Cos[fSyncE]) Sinc[epTau]
			) / (-1 + trt1 Cos[fSyncE])^2;
	b1Db1 = as b1Mask (dsLocDb1 / s0);
	(* derivatives to amplitudes*)
	daRe = mAmp sigRe;
	daIm = mAmp sigIm;
	(*real derivatives parameters for Gaussian and Lorentzian decay parameters*)
	dsmRe = -4 Pi^2 sm bt^2 abRe;
	dsmIm = -4 Pi^2 sm bt^2 abIm;
	dgmRe = -Pi bt abRe;
	dgmIm = -Pi bt abIm;
	(*complex derivatives for the shift including b1 part*)
	depRe = -2 Pi gyro bt abIm + b1Dep sigRe;
	depIm = 2 Pi gyro bt abRe + b1Dep sigIm;
	(*complex derivatives for the phase*)
	dphRe = abIm;
	dphIm = -abRe;
	(*derivative to b1*)
	db1Re = b1Db1 sigRe;
	db1Im = b1Db1 sigIm;

	(*combining all derivatives into the full jacobian*)
	dRe = Join[daRe, dsmRe, dgmRe, depRe, dphRe, db1Re];
	dIm = Join[daIm, dsmIm, dgmIm, depIm, dphIm, db1Im];

	(*output the signal and the jacobian*)
	sig = {Total[Join[abRe, abIm, 2]]};
	jac = Join[dRe, dIm, 2];
	(*First column is signal rest is jacobian*)
	Join[sig, jac]

], RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];


(* ::Subsubsection::Closed:: *)
(*ObjectiveFunc*)


ObjectiveFunc[data_, fids_, time_, gyro_, b1par_, b1Model_, p_?(MatrixQ[#, NumericQ] &)] := 
	Total[(ResidualFunc[data, fids, time, gyro, b1par, b1Model, p])^2] + 
		$AmaresB1Regularization (p[[-1, 1]] - 1.)^2


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


GradientFunc[data_, fids_, time_, gyro_, b1par_, b1Model_, p_?(MatrixQ[#, NumericQ] &), grp_] := Block[{sigjac, gradVec},
	sigjac = SignalJacCalcC[fids, time, p, b1Model, b1par, gyro];
	gradVec = -2 CollapseJacMat[Rest@sigjac, grp] . (data - First@sigjac);
	If[Length[grp] == 6, gradVec[[-1]] += 2 $AmaresB1Regularization (p[[-1, 1]] - 1.)];
	gradVec
]


(* ::Subsubsection::Closed:: *)
(*CollapseJacMat*)


CollapseJacMat[jac_, groupsId_] := Block[{np, pp, jacMat},
	{np, pp} = Dimensions[groupsId];
	jacMat = Partition[jac, pp];
	Flatten[Table[
		(Total[jacMat[[i, Flatten[Position[groupsId[[i]], #]]]]]) & /@ Range[Max[groupsId[[i]]]]
	, {i, 1, np}], 1]
]


(* ::Section:: *)
(*End Package*)


End[](* End Private Context *)

EndPackage[]
