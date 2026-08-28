(* ::Package:: *)

(* ::Title:: *)
(*QMRITools IVIMTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["QMRITools`IVIMTools`", Join[{"Developer`"}, Complement[QMRITools`$Contexts, {"QMRITools`IVIMTools`"}]]];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection::Closed:: *)
(*Functions*)


IVIMCalc::usage =
"IVIMCalc[data, bInp, init] calculates the IVIM fit.

data should be 1D ,2D, 3D or 4D. 
bInp should be full bmatrix which can be calculated from the bvecs en bvals using Bmatrix with the bvalues in s/mm^2. 
init should are the initialization parameters for 2 comp this is {s0, f, D, Dp} for 3 comp this is {s0, f1, f2, D, Dp1, Dp2}.

The fraction is defined between 0 and 1, the D, Dp, Dp1 and Dp2 is in mm^2/s.

output is {s0, f1, D, pD1} or {s0, f1, f2, D, pD1, pD2}."

IVIMFunction::usage =
"IVIMFunction[] gives the IVIM function with 2 comps.
IVIMFunction[comp] gives the IVIM function.
IVIMFunction[comp, type] gives the IVIM function. 

type can be \"Normal\" or \"Exp\".
comp can be 2 or 3.

output is the function with b, s0, f1, f2, D, pD1, pD2 as parameters. The fraction is defined between 0 and 1, the D, Dp, Dp1 and Dp2 is in mm^2/s."


IVIMCorrectData::usage = 
"IVIMCorrectData[data, {s0, f, pdc}, bval] removes the ivim signal from the data.

data is the original data.
{s0, f, pdc} are the solution to a 2 compartment IVIM fit using IVIMCalc or BayesianIVIMFit2.
bval are the bvalues.

The fraction is defined between 0 and 1, the pdc is in mm^2/s.

output is the corrected data."


IVIMResiduals::usage = 
"IVIMResiduals[data, bInp, pars] calculates the root mean square residuals of an IVIM fit using IVIMCalc, BayesianIVIMFit2 or BayesianIVIMFit3."

MeanBvalueSignal::usage = 
"MeanBvalueSignal[data, bval] calculates the geometric mean of the data for each unique bval. 
output is the mean data and the unique bvalues."


(* ::Subsection::Closed:: *)
(*General Options*)


IVIMConstrained::usage = 
"IVIMConstrained is an option for IVIMCalc. When set True the fit wil be constrained to the values given in IVIMConstrains."

IVIMTensFit::usage =
"IVIMTensFit is an option for IVIMCalc. When set True the tissue diffusion component wil be calculated as a tensor."

IVIMComponents::usage =
"IVIMComponents is an option for IVIMCalc. Default value is 2, the tissue and the blood component. can also be set to 3."

IVIMConstrains::usage = 
"IVIMConstrains is an option for IVIMCalc.
Default values are: {{0.8, 1.2}, {0, 1}, {0.0005, 0.0035}, {0.005, 0.5}, {0.002, 0.015}}.
Where {{s0 in percentage},{fractions},{tissue diffusion},{blood compartment Dp},{third compartment}}."

IVIMFixed::usage =
"IVIMFixed is an option for IVIMCalc and the default value is False. 
When set True the pseudo diffusion wil be fixed to the parameter given as init.
When set to \"One\" only the fast component of a 3 compartment fit is fixed."


FilterMaps::usage = 
"FilterMaps is an option for IVIMCorrectData. If True the IVIM parameter maps are filtered before signal correction."

FilterType::usage = 
"FilterType is an option for IVIMCorrectData. If FilterMaps is True it tells which filter to use. can be \"Median\" of \"Gaussian\"."

FilterSize::usage = 
"FilterSize is an option for IVIMCorrectData. If FilterMaps is True it gives the kernel size."


(* ::Subsection::Closed:: *)
(*Error Messages*)


IVIMCalc::init = "Number of initialization values is `1` and should be `2`."

IVIMCalc::time = "Number of comp is `1` but the number of relaxations times is `2`."

IVIMCalc::comp = "Number of comp should be 2 or 3 not `1`."

IVIMCalc::bvec = "The length of the data (`1`) should be the same as the length of the bmatrix (`2`)."


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*IVIMCalc*)


Options[IVIMCalc] = {
	Method -> Automatic, 
	Parallelize->True, 
	IVIMFixed -> True, 
	IVIMConstrained -> True, 
	IVIMComponents -> 2,
	IVIMConstrains -> {{0.1, 20}, {0, 1}, {0.0005, 0.0035}, {0.001, 0.5}, {0.001, 0.5}}
};

SyntaxInformation[IVIMCalc] = {"ArgumentsPattern" -> {_, _, _, _., OptionsPattern[]}};

IVIMCalc[data_, bInp_, init_, opts:OptionsPattern[]]:=IVIMCalc[data, bInp, init, False, opts]

IVIMCalc[data_, bInp_, init_, coil_, OptionsPattern[]] := Block[{
		comp, fix, con, method, s0Min, s0Max, fMin, fMax, dcMin, dcMax, pdc1Min, pdc1Max, pdc2Min, pdc2Max,
		depthD, mask, dat, coor, scale, coilFit, vox, start, dint, gradField, bvec, grad, bvecV, datM, bvecM,
		bound, invBound, allVars, selVars, startVals, s0e, f1e, f2e, dce, pdc1e, pdc2e,
		s0In, f1In, f2In, dcIn, pdc1In, pdc2In, func, vars, outVars, mapFunc, ivim, noVars,
		bm, s0, f1, f2, dc, pdc1, pdc2
	},

	(*get the options*)
	{comp, fix, con, method} = OptionValue[{IVIMComponents, IVIMFixed, IVIMConstrained, Method}];
	comp = Clip[comp, {1, 3}];
	{{s0Min, s0Max}, {fMin, fMax}, {dcMin, dcMax}, {pdc1Min, pdc1Max}, {pdc2Min, pdc2Max}} = OptionValue[IVIMConstrains];
	(*data checks*)
	depthD = ArrayDepth[data];
	If[depthD > 4, Return[Message[IVIMCalc::data, ArrayDepth[data]]]];

	(*convert data to vector if data is 2D or 3D or make data vector for 1D*)
	mask = Unitize@Mean@Ramp@Round[If[depthD==4, Transpose@data, data], .00001];
	Which[
		depthD >= 3, {dat, coor} = DataToVector[data, mask],
		depthD == 1, dat = {data},
		True, dat = data
	];

	(*normalize the data to the mean of the first value*)
	scale = Mean@dat[[All, 1]];
	dat = dat / scale;

	(*calculate the bmatrix using coil tensor if needed*)
	coilFit = If[coil =!= False && depthD >= 3,
		(*get the coil tensor*)
		{vox, start, dint} = coil;
		gradField = GradientCoilTensor[mask, vox, start, dint];
		{bvec, grad} = bInp;
		bvecV = BVector[bvec, grad, gradField];
		True
		,
		bvec = bInp;
		False
	];

	(*generate the signal vector and bvector per b-value*)
	{datM, bvecM} = MeanBvalueSignal[dat, bvec];
	bvecM = If[coilFit, First[MeanBvalueSignal[bvecV, bvec]], ConstantArray[bvecM, Length@dat]];
	If[Dimensions[datM]=!=Dimensions[bvecM], Return[Message[IVIMCalc::data, Dimensions@datM, Dimensions@bvecM]]];
	dat = Transpose[{bvecM, datM}, {3, 1, 2}];

	(*clear fit parameters*)
	Clear[bm, s0, f1, f2, dc, pdc1, pdc2];
	(*define bounds for the constraints*)
	bound[min_, max_, r_] := min + (max - min) LogisticSigmoid[r];
	invBound[min_, max_, x_] := With[{y = Clip[(x - min) / (max - min), {10.^-2, 1 - 10.^-2}]},	Log[y / (1 - y)]];

	(*the fit parameters*)
	allVars = {s0, f1, f2, dc, pdc1, pdc2};
	selVars = Join[
		{s0},
		If[comp > 1, {f1}, {}],
		If[comp > 2, {f2}, {}],
		{dc},
		If[comp > 1 && fix =!= True, {pdc1}, {}],
		If[comp > 2 && fix =!= True && fix =!= "One", {pdc2}, {}]
	];

	(*Initial fit values*)
	Switch[comp,
		1, {s0In, dcIn} = init;,
		2, {s0In, f1In, dcIn, pdc1In} = init;,
		3, {s0In, f1In, f2In, dcIn, pdc1In, pdc2In} = init;
	];
	s0In = If[s0In < 1.5, s0In, s0In / scale];

	startVals = selVars /. Thread[allVars -> {
		If[con, invBound[s0Min, s0Max, s0In], s0In],
		If[con, invBound[0, 1, f1In], f1In],
		If[con, invBound[0, 1, f2In / (1 - f1In)], f2In],
		If[con, invBound[dcMin, dcMax, dcIn], dcIn],
		If[con, invBound[pdc1Min, pdc1Max, pdc1In], pdc1In],
		If[con, invBound[pdc2Min, pdc2Max, pdc2In], pdc2In]
	}];

	(*physical params as functions of raw (unconstrained) fit variables*)
	s0e = If[con, bound[s0Min, s0Max, s0], s0];
	f1e = If[comp > 1, If[con, LogisticSigmoid[f1], f1], 0];
	f2e = If[comp > 2, If[con, (1 - f1e)*LogisticSigmoid[f2], f2], 0];
	dce = If[con, bound[dcMin, dcMax, dc], dc];
	pdc1e = If[comp > 1, If[fix === True, pdc1In, If[con, bound[pdc1Min, pdc1Max, pdc1], pdc1]], 0];
	pdc2e = If[comp > 2, If[fix === True || fix === "One", pdc2In, If[con, bound[pdc2Min, pdc2Max, pdc2], pdc2]], 0];

	(*generalized (constrained)IVIM model, expressed in raw vars*)
	func = s0e ((1 - f1e - f2e) Exp[-bm dce] + f1e Exp[-bm pdc1e] + f2e*Exp[-bm pdc2e]);
	vars = Thread[{selVars, startVals}];
	outVars = Join[{s0e, f1e, f2e}[[1 ;; comp]], {dce},
		If[comp == 1, {}, {pdc1e, pdc2e}[[1 ;; comp - 1]]]];
	noVars = ConstantArray[0., Length@outVars];

	(*perform the fitting using parallel if needed*)
	mapFunc = If[!(OptionValue[Parallelize] && depthD > 2), Map,
		DistributeDefinitions[noVars, outVars, func, vars];	ParallelMap];
	ivim = mapFunc[If[#[[1, 2]] < 0.1, noVars, Quiet[outVars /. FindFit[#, {func}, vars, bm]]] &, dat];

	(*rescale the signal and clip the values to logical range*)
	ivim[[All, 1]] = ivim[[All, 1]] scale;
	ivim[[All, 2 ;; comp]] = Clip[ivim[[All, 2 ;; comp]], {0., 1.}];
	ivim[[All, comp + 1]] = Clip[ivim[[All, comp + 1]], {0., 4./1000}];
	
	(*reconstruct the data*)
	If[depthD >= 3, ivim = VectorToData[ivim, coor]];
	Which[depthD==4 || depthD ===2, 
		Transpose@ivim, depthD == 1, 
		First@ivim, 
		True, ivim
	]
]


(* ::Subsection::Closed:: *)
(*IVIMFunction*)


SyntaxInformation[IVIMFunction] = {"ArgumentsPattern" -> {_, _}};

IVIMFunction[]:=IVIMFunction[2, "Normal"]

IVIMFunction[pars_]:=IVIMFunction[pars, "Normal"]

IVIMFunction[pars_, fun_] := Block[{
		func, Global`s0, Global`f1, Global`bm, Global`dc, Global`pdc1, 
		Global`f2, Global`pdc2, ff1, ff2, fdc, fpdc1, fpdc2
	},
	Switch[fun,
		"Normal", ff1 = Global`f1; ff2 = Global`f2; fdc = Global`dc; fpdc1 = Global`pdc1; fpdc2 = Global`pdc2;,
		"Exp", ff1 = Exp[Global`f1]/(1 + Exp[Global`f1]); ff2 = Exp[Global`f2]/(1 + Exp[Global`f2]); 
		fdc = Exp[Global`dc]; fpdc1 = Exp[Global`pdc1]; fpdc2 = Exp[Global`pdc2];
	];

	func = Switch[pars,
		2, Simplify[Global`s0*(((1 - ff1)*Exp[-Global`bm fdc]) + (ff1*Exp[-Global`bm fpdc1]))],
		3, Simplify[
		Global`s0*( (1 - ff1 - ff2)*Exp[-Global`bm*fdc] + ff1*Exp[-Global`bm*fpdc1] + 
			ff2*Exp[-Global`bm*fpdc2] )]
	]
]


(* ::Subsection::Closed:: *)
(*IVIMCorrectData*)


Options[IVIMCorrectData] = {FilterMaps -> True, FilterType -> "Median", FilterSize -> 1};

SyntaxInformation[IVIMCorrectData] = {"ArgumentsPattern" -> {_, {_, _,_} , _ , OptionsPattern[]}};

IVIMCorrectData[data_, {s0_, f_, pdc_}, bval_, OptionsPattern[]] := Module[{ff, pdcf, filt, dataSyn, dataCor},
	{ff, pdcf} = If[OptionValue[FilterMaps],
		filt = Switch[OptionValue[FilterType],
			"Median", MedianFilter,
			"Laplacian", LapFilt, 
			_, GaussianFilter
		];
		{filt[f, OptionValue[FilterSize]], filt[pdc, OptionValue[FilterSize]]}
		,
		{f, pdc}
	];

	dataSyn = Round@Clip[SynDataI[s0, ff, pdcf, bval], {0, Infinity}];
	dataCor = Clip[data - dataSyn, {0, 1.1 Max[data]}];

	{dataCor, dataSyn}
]

SynDataI = Compile[{{s0, _Real, 3}, {f, _Real, 3}, {pdc, _Real, 3}, {bval, _Real, 1}}, 
	Transpose[Map[(f s0 Exp[-# pdc]) &, bval]]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];

LapFilt[data_, fil_:0.8] := Clip[Chop[ImageData[TotalVariationFilter[Image3D[N@data, "Real"], fil, 
	Method -> "Laplacian", MaxIterations -> 15]]], MinMax[data]]


(* ::Subsection::Closed:: *)
(*IVIMResiduals*)


SyntaxInformation[MeanBvalueSignal] = {"ArgumentsPattern" -> {_, _}};

MeanBvalueSignal[data_, val_] := Block[{valU, pos, mean},
	{valU, pos} = UniqueBvalPosition[val];
	(*mean = Transpose[GeometricMean[Transpose[data[[All, #]]]] & /@ pos];*)
	mean = Transpose[Median[Transpose[data[[All, #]]]] & /@ pos];
	{mean, valU}
]


(* ::Subsection::Closed:: *)
(*IVIMResiduals*)

SyntaxInformation[IVIMResiduals] = {"ArgumentsPattern" -> {_, _, _}};

IVIMResiduals[data_, bInp_, pars_] := Block[{depthD,depthP,dat,par,res},
	(*data checks*)
	depthD = ArrayDepth[data];
	depthP = ArrayDepth[pars];

	dat = N[If[depthD == 4, Transpose[data], data]];
	dat = If[depthD > 1, RotateDimensionsLeft[dat], dat];
	par = If[depthP > 1, RotateDimensionsLeft[pars], pars];

	res = IVIMResCalcC[dat, bInp, par];

	res = If[depthD > 1, RotateDimensionsRight[res], res];

	Sqrt[Mean[Drop[res, 1]^2]] // N
];


IVIMResCalcC = Compile[{{dat, _Real, 1}, {bInp, _Real, 1}, {pars, _Real, 1}}, Block[{l, s0, f1, f2, dc, pdc1, pdc2, out}, 
	out = 0. dat;
	l = Length[pars];
	Which[
		l == 2,	{s0, dc} = pars;
		out = (s0*(((Exp[-bInp dc])))),
		l == 4,	{s0, f1, dc, pdc1} = pars;
		out = (s0*((((1 - f1) Exp[-bInp dc]) + (f1 Exp[-bInp pdc1])))),
		l == 6,	{s0, f1, f2, dc, pdc1, pdc2} = pars;
		out = (s0*((((1 - f1 - f2) Exp[-bInp dc]) + (f1 Exp[-bInp pdc1]) + (f2*Exp[-bInp pdc2]))))
	];
	dat - out
], RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
