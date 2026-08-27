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
"IVIMCalc[data, binp, init] calculates the IVIM fit.

data should be 1D ,2D, 3D or 4D. 
binp should be full bmatrix which can be calculated from the bvecs en bvals using Bmatrix with the bvalues in s/mm^2. 
init should are the initialization parameters for 2 components this is {s0, f, D, Dp} for 3 components this is {s0, f1, f2, D, Dp1, Dp2}.

The fraction is defined between 0 and 1, the D, Dp, Dp1 and Dp2 is in mm^2/s.

output is {s0, f1, D, pD1} or {s0, f1, f2, D, pD1, pD2}."

IVIMFunction::usage =
"IVIMFunction[] gives the IVIM function with 2 comps.
IVIMFunction[components] gives the IVIM function.
IVIMFunction[components, type] gives the IVIM function. 

type can be \"Normal\" or \"Exp\".
components can be 2 or 3.

output is the function with b, s0, f1, f2, D, pD1, pD2 as parameters. The fraction is defined between 0 and 1, the D, Dp, Dp1 and Dp2 is in mm^2/s."


IVIMCorrectData::usage = 
"IVIMCorrectData[data, {s0, f, pdc}, bval] removes the ivim signal from the data.

data is the original data.
{s0, f, pdc} are the solution to a 2 compartment IVIM fit using IVIMCalc or BayesianIVIMFit2.
bval are the bvalues.

The fraction is defined between 0 and 1, the pdc is in mm^2/s.

output is the corrected data."


IVIMResiduals::usage = 
"IVIMResiduals[data, binp, pars] calculates the root mean square residuals of an IVIM fit using IVIMCalc, BayesianIVIMFit2 or BayesianIVIMFit3."

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

IVIMCalc::time = "Number of components is `1` but the number of relaxations times is `2`."

IVIMCalc::comp = "Number of components should be 2 or 3 not `1`."

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
	IVIMConstrained -> False, 
	IVIMComponents -> 2,
	IVIMConstrains -> {{0.8, 1.2}, {0, 1}, {0.0005, 0.0035}, {0.001, 0.5}, {0.001, 0.5}}
};

SyntaxInformation[IVIMCalc] = {"ArgumentsPattern" -> {_, _, _, _., OptionsPattern[]}};

IVIMCalc[data_, binp_, init_, opts:OptionsPattern[]]:=IVIMCalc[data, binp, init, False, opts]

IVIMCalc[data_, binp_, init_, coil_, OptionsPattern[]] := Block[{
		depthD, dirD, mask, dat, coor, mdat, datn, cfit, vox, start, dint, gradField, 
		bvec, grad, bvecV, bvecf, components, fixed, constrained, method,
		s0min, s0max, fmin, fmax, dcmin, dcmax, pdc1min, pdc1max, pdc2min, pdc2max,
		fVars, pdcVars, fixValues, fixRule, pdcmin, pdcmax, activePdc, func, cons,
		vars, funcf, fpars, out, mapfun, ivim, s0s, frin1, frin2, dcin, pdcin1, pdcin2
	},

	(*data checks*)
	depthD = ArrayDepth[data];
	If[depthD > 4, Return[Message[IVIMCalc::data, ArrayDepth[data]]]];

	(*convert data to vector if data is 2D or 3D or make data vector for 1D*)
	mask = Unitize@Mean@If[depthD==4, Transpose@data, data];
	Which[
		depthD >= 3, {dat, coor} = DataToVector[data, mask],
		depthD == 1, dat = {data},
		True, dat = data
	];
	mdat = First[Mean@dat];
	datn = dat/mdat;

	(*calculate the bmatrix using coil tensor if needed*)
	cfit = If[coil =!= False && depthD >= 3,
		(*get the coil tensor*)
		{vox, start, dint} = coil;
		gradField = GradientCoilTensor[mask, vox, start, dint];
		{bvec, grad} = binp;
		bvecV = BVector[bvec, grad, gradField];
		True
		,
		bvec = binp;
		False
	];

	(*generate the signal vector and bvector per b-value*)
	{datn, bvecf} = MeanBvalueSignal[datn, bvec];
	bvecf = If[cfit,
		First[MeanBvalueSignal[bvecV, bvec]],
		ConstantArray[bvecf, Length@datn]
	];
	If[Dimensions[datn]=!=Dimensions[bvecf], Return[Message[IVIMCalc::data, Dimensions@datn, Dimensions@bvecf]]];

	(*get the options*)
	{components, fixed, constrained, method} = OptionValue[{IVIMComponents, IVIMFixed, IVIMConstrained, Method}];
	components = Clip[components, {1, 3}];
	{{s0min, s0max}, {fmin, fmax}, {dcmin, dcmax}, {pdc1min, pdc1max}, {pdc2min, pdc2max}} = OptionValue[IVIMConstrains];

	Clear[bm, f1, f2, s0, dc, pdc1, pdc2];

	(*Indices for pdc and f based on components*)
	fVars = If[components > 1, {f1, f2}[[;; components - 1]], {}];
	pdcVars = If[components > 1, {pdc1, pdc2}[[;; components - 1]], {}];

	Switch[components,
		1, f1 = f2 = 0; {s0s, dcin} = init;,
		2, f2 = 0; {s0s, frin1, dcin, pdcin1} = init;,
		3, {s0s, frin1, frin2, dcin, pdcin1, pdcin2} = init;
	];

	(*Fixed values for pdc*)
	fixValues = {pdc1 -> pdcin1, pdc2 -> pdcin2};
	fixRule = Switch[fixed,
		True, Take[fixValues, components - 1],
		"One", {pdc2 -> pdcin2},
		_, {}];
	{pdcmin, pdcmax} = Switch[fixed,
		False, {Take[{pdc1min, pdc2min}, components - 1], Take[{pdc1max, pdc2max}, components - 1]},
		"One", {pdc1min, pdc1max},
		_, {{}, {}}];
	activePdc = Complement[pdcVars, fixRule[[All, 1]]];

	(*Generalized IVIM model*)
	func = s0*((1 - Total[fVars])*Exp[-bm*dc] + Total[fVars*Exp[-bm*pdcVars]]) // Simplify // ReplaceAll[fixRule];

	(*fitting contrains if needed*)
	cons = If[constrained,
	Flatten[{
		{s0 > 0, dcmin < dc < dcmax},
		If[components > 1, Total[fVars] <= 1, {}],
		If[components > 1, Thread[fmin < fVars < fmax], {}],
		{Less @@ Flatten[{dc, pdcVars /. fixRule}]},
		If[fixed === True, {}, Thread[pdcmin < activePdc < pdcmax]]
	}], {}];

	(*fitting variables and function*)
	vars = Flatten[{s0, fVars, dc, activePdc}];
	vars = Thread[{vars, vars /. Thread[{s0, f1, f2, dc, pdc1, pdc2} -> {s0s, frin1, frin2, dcin, pdcin1, pdcin2}]}];
	funcf = If[constrained, {func, cons}, func];
	fpars = {bm};

	(*define output*)
	out = Join[
		{s0, f1, f2}[[1 ;; components]], {dc}, 
		If[components == 1, {}, {pdc1, pdc2}[[1 ;; components - 1]]]
	] /. fixRule;

	mapfun = If[OptionValue[Parallelize] && depthD>1,
		DistributeDefinitions[out, funcf, vars, fpars];
		ParallelMap,
		Map];

	ivim = mapfun[Quiet[
		out /. FindFit[Thread[#], funcf, vars, fpars, MaxIterations -> 150]
	] &,Thread[{bvecf, datn}]];

	ivim[[All, 1]] = ivim[[All, 1]] mdat;
	ivim[[All, 2 ;; components]] = Clip[ivim[[All, 2 ;; components]], {0., 1.}];
	ivim[[All, components + 1]] = Clip[ivim[[All, components + 1]], {0., 4./1000}];
	
	If[depthD >= 3, ivim = VectorToData[ivim, coor]];
	If[depthD==4 || depthD ===2, Transpose@ivim, ivim]
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
		{filt[f, OptionValue[FilterSize]],filt[pdc, OptionValue[FilterSize]]}
		,
		{f, pdc}
	];

	dataSyn = Round@Clip[SynDatai[s0, ff, pdcf, bval], {0, Infinity}];
	dataCor = Clip[data - dataSyn, {0, 1.1 Max[data]}];

	{dataCor, dataSyn}
]

SynDatai = Compile[{{s0, _Real, 3}, {f, _Real, 3}, {pdc, _Real, 3}, {bval, _Real, 1}}, 
	Transpose[Map[(f s0 Exp[-# pdc]) &, bval]]
, RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];

LapFilt[data_, fil_:0.8] := Clip[Chop[ImageData[TotalVariationFilter[Image3D[N@data, "Real"], fil, 
	Method -> "Laplacian", MaxIterations -> 15]]], MinMax[data]]


(* ::Subsection::Closed:: *)
(*IVIMResiduals*)


SyntaxInformation[MeanBvalueSignal] = {"ArgumentsPattern" -> {_, _}};

MeanBvalueSignal[data_, val_] := Block[{valU, pos, mean},
	{valU, pos} = UniqueBvalPosition[val];
	mean = Transpose[GeometricMean[Transpose[data[[All, #]]]] & /@ pos];
	{mean, valU}
]


(* ::Subsection::Closed:: *)
(*IVIMResiduals*)

SyntaxInformation[IVIMResiduals] = {"ArgumentsPattern" -> {_, _, _}};

IVIMResiduals[data_, binp_, pars_] := Block[{depthD,depthP,dat,par,res},
	(*data checks*)
	depthD = ArrayDepth[data];
	depthP = ArrayDepth[pars];

	dat = N[If[depthD == 4, Transpose[data], data]];
	dat = If[depthD > 1, RotateDimensionsLeft[dat], dat];
	par = If[depthP > 1, RotateDimensionsLeft[pars], pars];

	res = IVIMResCalcC[dat, binp, par];

	res = If[depthD > 1, RotateDimensionsRight[res], res];

	Sqrt[Mean[Drop[res, 1]^2]] // N
];


IVIMResCalcC = Compile[{{dat, _Real, 1}, {binp, _Real, 1}, {pars, _Real, 1}}, Block[{s0, f1, f2, dc, pdc1, pdc2, out}, 
	out = 0. dat;
	out = Switch[Length[pars],
		2,
		{s0, dc} = pars;
		(s0*(((Exp[-binp dc])))),
		4,
		{s0, f1, dc, pdc1} = pars;
		(s0*((((1 - f1) Exp[-binp dc]) + (f1 Exp[-binp pdc1])))),
		6,
		{s0, f1, f2, dc, pdc1, pdc2} = pars;
		(s0*((((1 - f1 - f2) Exp[-binp dc]) + (f1 Exp[-binp pdc1]) + (f2*Exp[-binp pdc2]))))
	];
	dat - out
], RuntimeAttributes -> {Listable}, RuntimeOptions -> {"Speed", "WarningMessages" -> False}];


(* ::Section:: *)
(*End Package*)


End[]

EndPackage[]
