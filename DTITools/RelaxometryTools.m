(* ::Package:: *)

(* ::Title:: *)
(*DTITools RelaxometryTools*)


(* ::Subtitle:: *)
(*Written by: Martijn Froeling, PhD*)
(*m.froeling@gmail.com*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["DTITools`RelaxometryTools`", {"Developer`"}];
$ContextPath=Union[$ContextPath,System`$DTIToolsContextPaths];

Unprotect @@ Names["DTITools`RelaxometryTools`*"];
ClearAll @@ Names["DTITools`RelaxometryTools`*"];


(* ::Section:: *)
(*Usage Notes*)


(* ::Subsection:: *)
(*Functions*)


T1rhoFit::usage = 
"T1rhoFit[data, EchoTimes]"

T2Fit::usage = 
"T2Fit[data, EchoTimes]"

TriExponentialT2Fit::usage = 
"TriExponentialT2Fit[data, EchoTimes] fits the T2 based on Azzabou N et.al. Validation of a generic approach to muscle water 
T2 determination at 3T in fat-infiltrated skeletal muscle. J. Magn. Reson. 2015."

EPGSignal::usage = 
"EPGSignal[Necho, echoSpace, T1, T2, angle, B1] generates a EPG T2 curve with stimulated echos. T1, T2 and echoSpace are in ms, angel is in degree, B1 is between 0 and 1."

EPGT2Fit::usage = 
"EPGT2Fit[data, EchoTimes, angle, relax] fits the T2 based on Marty B et.al. Simultaneous muscle water T2 and fat fraction mapping using transverse relaxometry with stimulated echo compensation.
angle is the refocussing angle in degree."

CreateT2Dictionary::usage = 
"CreateT2Dictionary[{T1m, T1f, T2f}, {Necho, echoSpace, angle}] Creates a EPG signal dictionary used for EPGT2fit."

DictionaryMinSearch::usage = 
"DictionaryMinSearch[dictionary, y] performs dictionary minimization of data y. dictionary is generated with CreateT2Dictionary."

NonLinearEPGFit::usage = 
"NonLinearEPGFit[{vals, T2cons}, y] performs dictionary minimization of data y. vals = {{T1muscle, T1fat, T2fat}, {Necho, echoSpace, angle}}."

CalibrateEPGT2Fit::usage = 
"CalibrateEPGT2Fit[datan, times, angle] calculates the Fat T2 ralaxation that will be used in the EPGT2fit."



(* ::Subsection:: *)
(*General Options*)


DictT2Range::usage = 
"DictT2Range is an option for CreateT2Dictionary and EPGT2Fit. is specifies the range and step of the T2 values in the dictionary {min, max, step} in ms."

DictB1Range::usage = 
"DictB1Range is an option for CreateT2Dictionary and EPGT2Fit. It specifies the range and step of the B1 values in the dictionary {min, max, step}."

EPGRelaxPars::usage = 
"EPGRelaxPars is and option for EPGT2Fit. Needs to be {T1muscl, T1Fat, T2Fat} in ms, defaul is {1400,365,137}."

MonitorEPGFit::usage = 
"MonitorEPGFit show waitbar during EPGT2Fit."

EPGFitPoints::usage = 
"EPGFitPoints is a option for CalibrateEPGT2Fit and EPGT2Fit. Number of points is 200 by default."

EPGCalibrate::usage = 
"EPGCalibrate is an option for EPGT2Fit. If set to True it does autmatic callibration of the T2 fat relaxation time."

OutputCalibration::usage = 
"OuputCalibration is an option for EPGT2Fit and TriExponentialT2Fit. If true it outputs the calibartion values."


(* ::Subsection:: *)
(*Error Messages*)


(* ::Section:: *)
(*Functions*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*T1rhoFit*)


Options[T1rhoFit] = {Method -> "Linear"};

SyntaxInformation[T1rhoFit]= {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

T1rhoFit[datan_, times_, OptionsPattern[]] := 
 Switch[OptionValue[Method],
  "Linear", LinFit[datan, times],
  _, LogFit[datan, times]
  ]


(* ::Subsection::Closed:: *)
(*T2Fit*)


Options[T2Fit] = {Method -> "Linear"};

SyntaxInformation[T2Fit]= {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

T2Fit[datan_, times_, OptionsPattern[]] := Switch[OptionValue[Method],
  "Linear", LinFit[N[datan], times],
  _, LogFit[N[datan], times]
  ]


(* ::Subsection:: *)
(*linFit and LogFit*)


(* ::Subsubsection::Closed:: *)
(*LinFit*)


LinFit[datan_, times_] := 
 Module[{datal, result, fdat, offset, T1r, t, ad,off,t1rho},
  ad = ArrayDepth[datan];
  datal = Log[datan] /. {Indeterminate -> 0, -Infinity -> 0} // N;
  datal = Switch[ad,
    3, Transpose[datal, {3, 1, 2}],
    4, Transpose[datal, {1, 4, 2, 3}]
    ];
  
  result = ParallelMap[If[Total[#]==0.,
  		{0.,0.},
  		fdat = Transpose[{times, #}];
  		(*Quiet[LinearModelFit[fdat, t, t]["BestFitParameters"]]*)
  		{off, t1rho} /. Quiet[FindFit[fdat, off + t1rho t, {off, t1rho}, t]]
  	]&, datal, {ad - 1}];
  
  {offset, T1r} = Switch[ad,
    3, Transpose[result, {2, 3, 1}],
    4, Transpose[result, {2, 3, 4, 1}]
    ];
  
  {offset, T1r} = {Exp[offset], -1/(T1r /. 0. -> Infinity)};
  T1r = Clip[T1r, {0, 500}, {0, 500}];
  {offset, T1r}
  ]


(* ::Subsubsection::Closed:: *)
(*LogFit*)


LogFit[datan_, times_] := 
 Module[{result, fdat, offset, T1r, off, t1rho, t, ad, datal},
  ad = ArrayDepth[datan];
  datal = Switch[ad,
    3, Transpose[datan, {3, 1, 2}],
    4, Transpose[datan, {1, 4, 2, 3}]
    ];
  
  result = ParallelMap[
    (
      fdat = Transpose[{times, #}];
      {off, t1rho} /. Quiet[NonlinearModelFit[
          fdat,
          {off Exp[-t/t1rho], 10 < t1rho < 301},
          {{off, Last@Last@fdat}, {t1rho, 100}},
          t]["BestFitParameters"]]
      ) &, datal, {ad - 1}];
  
  {offset, T1r} = Switch[ad,
    3, Transpose[result, {2, 3, 1}],
    4, Transpose[result, {2, 3, 4, 1}]
    ];
  
  T1r = Clip[T1r, {0, 500}, {0, 500}];
  {offset, T1r}
  ]


(* ::Subsection:: *)
(*EPGSignal*)


(* ::Subsubsection::Closed:: *)
(*EPGSignal*)


SyntaxInformation[EPGSignal]= {"ArgumentsPattern" -> {_, _, _, _, _, _, OptionsPattern[]}}

EPGSignal[Nechoi_, echoSpace_, T1_, T2_, angle_, B1_] := EPGSignali[Nechoi, echoSpace, T1, T2, angle, B1]
  
EPGSignali[Nechoi_, echoSpace_, T1_, T2_, angle_, B1_] := Block[
  {tau, T0, R0, alpha, Smat, Tmat, Rmat, Pmat, Emat, xvec, 
   out, t2r, t1r,Necho},
  (*define internal paramters*)
  Necho=Round[Nechoi];
  alpha = N[B1 angle Degree];
  tau = echoSpace/2.;
  t2r = Exp[-tau/T2];
  t1r = Exp[-tau/T1];
  
  (*Selection matrix to move all traverse states up one coherence level*)
  Smat = MixMatrix[Necho];
  (*RF mixing matrix*)
  T0 = RotMatrixT[alpha];
  Tmat = MakeDiagMat[T0, Necho, 1];
  (* Relaxation matrix*)
  R0 = DiagonalMatrix[N[{t2r, t2r, t1r}]];
  Rmat = MakeDiagMat[R0, Necho, t2r];
  (*Precession and relaxation matrix*)
  Pmat = Rmat.Smat;
  (*Matrix representing the inter-echo duration*)
  Emat = Pmat.Tmat.Pmat;
  (*Recursively apply matrix to get echo amplitudes*)
  xvec = ConstantArray[0, 3*Necho + 1]; xvec[[1]] = 1;
  (*create outpu*)
  out = Table[xvec = Emat.xvec; First[xvec], {i, Necho}];
  out/out[[1]]
  ]


(* ::Subsubsection::Closed:: *)
(*MakeDiagMat*)


MakeDiagMat[mat_, Necho_, first_] := Block[{out},
   out = ArrayPad[Flatten[Flatten[(DiagonalMatrix[ConstantArray[1, Necho]]*ConstantArray[mat, {Necho, Necho}]), {2, 4}], {2, 3}], {{1, 0}, {1, 0}}];
   out[[1, 1]] = first;
   out
   ];


(* ::Subsubsection::Closed:: *)
(*MixMatrix*)


(*if run once with Necho definition is stored*)
MixMatrix[Necho_] := MixMatrix[Necho] = Block[{len, Smat, off1, off2},
   len = 3*Necho + 1;
   (*mixing matirx*)
   Smat = ConstantArray[0, {len, len}];
   Smat[[1, 3]] = Smat[[2, 1]] = Smat[[3, 6]] = Smat[[4, 4]] = 1;
   Table[
    off1 = ((o - 1) - 1)*3 + 2;
    If[off1 <= len, Smat[[3*o - 1, off1]] = 1];
    off2 = ((o + 1) - 1)*3 + 3;
    If[off2 <= len, Smat[[3*o, off2]] = 1];
    Smat[[3*o + 1, 3*o + 1]] = 1;
    , {o, 2, Necho}];
   (*output mixing matrix*)
   Smat
   ]


(* ::Subsubsection::Closed:: *)
(*RotMatrixT*)


(*if run once with alpha definition is stored*)
RotMatrixT[alpha_] := RotMatrixT[alpha] = {
    {Cos[alpha/2]^2, Sin[alpha/2]^2, Sin[alpha]},
    {Sin[alpha/2]^2, Cos[alpha/2]^2, -Sin[alpha]},
    {-0.5 Sin[alpha], 0.5 Sin[alpha], Cos[alpha]}
    };


(* ::Subsection:: *)
(*EPGT2Fit*)


(* ::Subsubsection::Closed:: *)
(*EPGT2Fit*)


Options[EPGT2Fit]= {DictT2Range -> {20., 80., 0.3}, DictB1Range -> {0.4, 1., 0.02}, 
	EPGRelaxPars->{1400., 365., 137.},Method->"dictionary",MonitorEPGFit->True, 
	EPGCalibrate->True, EPGFitPoints -> 200,OutputCalibration->False}

SyntaxInformation[EPGT2Fit]= {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}}

EPGT2Fit[datan_, times_, angle_, OptionsPattern[]]:=Block[{Necho,echoSpace,T1m, T1f, T2f, ad, datal, sol,
	wat, fat, fatMap, T2map, B1Map, clip, cons, B1i, T2i, T2s, B1s, soli, error, 
	dictf, valsf,ydat,fwf,residualError,T2mc,B1c,cal},
	
	{T1m, T1f, T2f} = OptionValue[EPGRelaxPars];
	clip = OptionValue[DictT2Range][[1;;2]];
	
	If[OptionValue[EPGCalibrate],
		Print["Callibrating EPG fat T2."];
		cal = {T2mc, T2f, B1c} = CalibrateEPGT2Fit[datan, times, angle, 
			EPGRelaxPars -> {clip, {50, 300}, {T1m, T1f}},EPGFitPoints->OptionValue[EPGFitPoints]];
		T2f=N@Round@T2f;
		Print["EPG fat callibration:  ", T2f, " ms"];
	];
  
  	ad = ArrayDepth[datan];
	datal = N@Switch[ad,
		3, Transpose[datan, {3, 1, 2}],
		4, Transpose[datan, {1, 4, 2, 3}]
	];
	
	Necho = Length[times];
	echoSpace = First[times];

	(*monitor calculation*)
	i=j=0;
	SetSharedVariable[i];ParallelEvaluate[j=0];
	If[OptionValue[MonitorEPGFit], PrintTemporary[ProgressIndicator[Dynamic[i],{ 0,Times@@Dimensions[datan[[All,1]]] }]] ];
	
	sol = Switch[OptionValue[Method],
		"dictionary",
	  	(*create the dictionary*)
		{dictf, valsf} = CreateT2Dictionary[{T1m, T1f, T2f}, {Necho, echoSpace, angle}, 
			DictB1Range->OptionValue[DictB1Range], DictT2Range->OptionValue[DictT2Range]];
		(*monitor calculation*) 
		PrintTemporary["starting dictionary min search: ",DateString[]];
		(*perform the fit using parallel kernels*)
		DistributeDefinitions[dictf, valsf,LeastSquaresC,ErrorC,DictionaryMinSearchi];
		ParallelMap[(
			(*monitor calculation*)
			j++;If[j>100,i+=j;j=1;];
			DictionaryMinSearchi[{dictf,valsf},#]
		) &, datal, {ad - 1}]
		
		,"NLLS",
	  	(*define fit values*)
	  	valsf = {{T1m, T1f, T2f}, {Necho, echoSpace, angle}};
		(*monitor calculation*)
		PrintTemporary["starting NLLS fitting: ",DateString[]];
		(*perform the fit using parallel kernels*)
		DistributeDefinitions[ErrorFunc, LeastSquaresC, ErrorC, EPGSignali, NonLinearEPGFiti, valsf, clip];
		ParallelMap[(
			(*monitor calculation*)
			j++;If[j>100,i+=j;j=1;];
			NonLinearEPGFiti[{valsf,clip},#]
		) &, datal, {ad - 1}]
	];
	
	(*restructure fit solution*)
	sol = TransData[sol,"r"];
	
	Print[Dimensions[sol]];
		
	{wat, fat} = Clip[{sol[[3]], sol[[4]]}, {0, Max[sol[[3;;4]]]}];
  	fatMap = Clip[fat / ((wat + fat) /. 0. -> Infinity), {0, 1}];
  	T2map = sol[[1]];
  	B1Map = MedianFilter[sol[[2]], 1];
  	error = sol[[5]];
  		
  	(*if neede also output callibaration*)
  	If[OptionValue[OutputCalibration],
  		{{{T2map,B1Map},{wat, fat, fatMap}},cal},
  		{{T2map,B1Map},{wat, fat, fatMap}}
  	]
]


(* ::Subsubsection::Closed:: *)
(*NonLinearEPGFit*)


SyntaxInformation[NonLinearEPGFit]= {"ArgumentsPattern" -> {_, _}}

NonLinearEPGFit[{valsf_,cons_}, yi_] := NonLinearEPGFiti[{valsf, cons}, N[yi]]

NonLinearEPGFiti[{_,_}, {0. ..}] = {0., 0., 0., 0., 0.};

NonLinearEPGFiti[{valsf_,cons_}, ydat_] := Block[
	{fwf, residualError, soli, T1m, T1f, T2f, Necho, echoSpace, angle, T2s, B1s},
	(*perform the fit*)
	{residualError, soli} = Quiet@FindMinimum[{
			ErrorFunc[ydat, T2i, B1i, valsf],  
			{0.4 <= B1i, B1i <= 1, cons[[1]] <= T2i, T2i <= cons[[2]]}
		},{{T2i,35}, {B1i, 0.8}}, 	
		AccuracyGoal -> 5, PrecisionGoal -> 5, MaxIterations -> 25];
	{T2s, B1s} = {T2i, B1i} /. soli;
	(*get corresponding fat fractions*)
	{{T1m, T1f, T2f}, {Necho, echoSpace, angle}} = valsf;	
	fwf = LeastSquaresC[Transpose[{
		EPGSignali[Necho, echoSpace, T1m, T2s, angle, B1s],
		EPGSignali[Necho, echoSpace, T1f, T2f, angle, B1s]
		}], ydat];
	(*export paramters*)
	Flatten[{{T2s, B1s}, fwf, residualError/Total[fwf]^2}]
   ]


(* ::Subsubsection::Closed:: *)
(*DictionaryMinSearch*)


SyntaxInformation[DictionaryMinSearch]= {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

DictionaryMinSearch[{dictf_, valsf_}, yi_] := DictionaryMinSearchi[{dictf, valsf}, N[yi]]

DictionaryMinSearchi[{_, _}, {0. ..}] = {0., 0., 0., 0., 0.};

DictionaryMinSearchi[{dictf_, valsf_}, ydat_] := Block[{fwf, residualError},
	(*calcualte dictionary error*)
  	fwf = LeastSquaresC[dictf, ydat];
  	residualError = ErrorC[ydat, fwf, dictf];
  	(*find Min value*)
  	Flatten[{valsf, fwf, residualError}[[All, First@Ordering[residualError, 1]]]]
   ]


(* ::Subsubsection::Closed:: *)
(*ErrorFunc*)


ErrorFunc[y_, T2m_Real, B1_Real, vals_] := Quiet@Block[
	{sig, f,T1m,T1f,T2f,Necho,echoSpace,angle},
		{{T1m,T1f,T2f},{Necho,echoSpace,angle}}=vals;
		sig = Transpose[{
			EPGSignali[Necho, echoSpace, T1m, T2m, angle, B1],
			EPGSignali[Necho, echoSpace, T1f, T2f, angle, B1]
			}];
		LeastSquaresErrorC[sig, y]
   ]

ErrorFunc[y_, T2m_Real, T2f_Real, B1_Real, vals_] := Quiet@Block[
	{T1m,T1f,Necho,echoSpace,angle},
		{{T1m,T1f},{Necho,echoSpace,angle}}=vals;
		LeastSquaresErrorC[Transpose[{
			EPGSignali[Necho, echoSpace, T1m, T2m, angle, B1],
			EPGSignali[Necho, echoSpace, T1f, T2f, angle, B1]
			}], y]
   ]


(* ::Subsubsection::Closed:: *)
(*LeastSquaresC*)


LeastSquaresC = Compile[{{A, _Real, 2}, {y, _Real, 1}}, Block[{T = Transpose[A]}, 
	(Inverse[T.A].T).y], 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", CompilationTarget -> System`$DTIToolsCompiler];


(* ::Subsubsection::Closed:: *)
(*LeastSquaresErrorC*)


LeastSquaresErrorC = Compile[{{A, _Real, 2}, {y, _Real, 1}}, Block[{T = Transpose[A]},
    Total[(y - A.(Inverse[T.A].T).y)^2]], 
    RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", CompilationTarget -> System`$DTIToolsCompiler];


(* ::Subsubsection::Closed:: *)
(*ErrorC*)


ErrorC = Compile[{{y, _Real, 1}, {f, _Real, 1}, {A, _Real, 2}}, 
	Total[(y - A.f)^2], 
	RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed", CompilationTarget -> System`$DTIToolsCompiler];


(* ::Subsection::Closed:: *)
(*CalibrateEPGT2Fit*)


Options[CalibrateEPGT2Fit] = {EPGRelaxPars -> {{20, 80}, {50, 300}, {1400., 365.}}, EPGFitPoints -> 200};

CalibrateEPGT2Fit[datan_, times_, angle_, OptionsPattern[]] := 
 Block[{ad,Necho,echoSpace,maskT2,dataT2,fmask,fitData,step,
 	T2mmin, T2mmax, T2fmin, T2fmax, T1m, T1f, valsf,cons,
 	soli,fits,residualError},
  
  ad = ArrayDepth[datan];
  
  Necho = Length[times];
  echoSpace = First[times];
  
  Switch[ad,
   3,
   (*single slice*)
   (*make mask an normalize data to first echo*)
   maskT2 = SmoothMask[{Mask[Mean[datan], {2}]}, MaskComponents -> 2, MaskClosing -> 1][[1]];
   dataT2 = maskT2 # & /@ datan;
   dataT2 = dataT2/MeanNoZero[Flatten[dataT2]];
   (*create mask selecting fat*)
   fmask = Mask[dataT2[[-1]], {0.5}];
   fmask = ImageData[SelectComponents[Image[fmask], "Count", -2]];
   (*data for calibration fit*)
   fitData = Transpose[Flatten[GetMaskData[#, fmask]] & /@ dataT2]
   ,
   4,
   (*mulit slice*)
   (*make mask an normalize data to first echo*)
   maskT2 = SmoothMask[Mask[Mean[Transpose[datan]], {2}], MaskComponents -> 2, MaskClosing -> 1];
   dataT2 = NormalizeData[MaskDTIdata[datan, maskT2]];
   (*create mask selecting fat*)
   fmask = Mask[dataT2[[All, -1]], {0.5}];
   fmask = ImageData[SelectComponents[Image3D[fmask], "Count", -2]];
   (*data for calibration fit*)
   fitData = Transpose[Flatten[GetMaskData[#, fmask]] & /@ Transpose@dataT2]
   ];
  
  step = Ceiling[Length[fitData]/OptionValue[EPGFitPoints]];
  {{T2mmin, T2mmax}, {T2fmin, T2fmax}, {T1m, T1f}} = OptionValue[EPGRelaxPars];
  
  (*define fit values*)
  valsf = {{T1m, T1f}, {Necho, echoSpace, angle}};
  cons = {0.4 <= B1i, B1i <= 1, T2mmin <= T2mi, T2mi <= T2mmax, T2fmin <= T2fi, T2fi <= T2fmax};
  
  (*perform the fit using parallel kernels*)
  DistributeDefinitions[ErrorFunc, LeastSquaresC, ErrorC, EPGSignali, valsf, cons];
  fits = ParallelMap[(
      (*calcualte NLLS error*)
      {residualError, soli} = Quiet@FindMinimum[{ErrorFunc[#, T2mi, T2fi, B1i, valsf], cons}, {{T2mi, 35.}, {T2fi, 137.}, {B1i, 0.8}},
       	AccuracyGoal -> 5, PrecisionGoal -> 5, MaxIterations -> 25];
      soli[[All, 2]]
      ) &, fitData[[1 ;; ;; step]]];
  
  Mean[fits]
  ]


(* ::Subsection::Closed:: *)
(*CreateT2Dictionary*)


Options[CreateT2Dictionary] = {DictT2Range -> {20., 80., 0.3}, DictB1Range -> {0.4, 1., 0.02}};

SyntaxInformation[CreateT2Dictionary]= {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

CreateT2Dictionary[relax_, ang_, OptionsPattern[]] := CreateT2Dictionaryi[N[relax], N[ang], N[OptionValue[DictT2Range]],N[OptionValue[DictB1Range]]]

(*save each unique dictionary*)
CreateT2Dictionaryi[relax_, ang_, T2range_, B1range_] := CreateT2Dictionaryi[relax, ang, T2range, B1range] = Block[
   {dict, dictf, valsf, T1m, T1f, T2f, Necho, echoSpace, angle, t2s, t2e, t2i, b1s, b1e, b1i},
   (*set parameters*)
   {T1m, T1f, T2f} = relax;
   {Necho, echoSpace, angle} = ang;
   (*get dictionary values*)
   {t2s, t2e, t2i} = T2range;
   {b1s, b1e, b1i} = B1range;
   (*Print["Creating dictionary for new values :", 
   "\n		{T1muscle, T1fat, T2fat} =",	{T1m, T1f, T2f}, 
   "\n		{Necho, echoSpace, ange} =", {Necho, echoSpace, angle}, 
   "\n		{T2min, T2max, T2step} =", {t2s, t2e, t2i},
   "\n		{B1min, B1max, B1step} =", {b1s, b1e, b1i}];*)
   (*create dictionary*)
   dict = Table[
     {{EPGSignal[Necho, echoSpace, T1m, T2m, angle, B1], EPGSignal[Necho, echoSpace, T1f, T2f, angle, B1]}, {T2m, B1}}, 
     	{T2m, t2s, t2e, t2i}, {B1, b1s, b1e, b1i}];
   (*flatten dictionary*)  	
   valsf = Flatten[dict[[All, All, 2]], 1];
   dictf = Transpose /@ Flatten[dict[[All, All, 1]], 1];
   Print["The dictionary contains "<>ToString[Length[dictf]]<>" values."];
   {dictf, valsf}
   ]


(* ::Subsection::Closed:: *)
(*TriExponentialT2Fit*)


Options[TriExponentialT2Fit]={OutputCalibration->False}

SyntaxInformation[TriExponentialT2Fit]= {"ArgumentsPattern" -> {_, _, OptionsPattern[]}}

TriExponentialT2Fit[datan_, times_,OptionsPattern[]] := 
 Block[{result, fdat, offset, T1r, off, t1rho, t, ad, datal,
   model, cs, cf, T2s, T2f, Af, Am, T2m, x, Aff, Amm, T2mf, S0, ffr, 
   mfr, T2, model2,cal,
   maskT2, dataT2, fmask, fitData,
   Afi, Ami, csi, T2mi, T2fi, T2si, sci
   },
  
  ad = ArrayDepth[datan];
  
  Switch[ad,
   
   3,(*single slice*)
   (*make mask an normalize data to first echo*)
   maskT2 = SmoothMask[{Mask[Mean[datan], {2}]}, MaskComponents -> 2, MaskClosing -> 1][[1]];
   dataT2 = maskT2 # & /@ datan;
   dataT2 = dataT2/MeanNoZero[Flatten[dataT2]];
   (*create mask selecting fat*)
   fmask = Mask[dataT2[[-1]], {0.4}];
   fmask = ImageData[SelectComponents[Image[fmask], "Count", -2]];
   (*data for calibration fit*)
   fitData = Transpose[{times, Mean[Flatten[GetMaskData[#, fmask]]] & /@ dataT2}];
   datal = Transpose[dataT2, {3, 1, 2}],
   
   4,(*mulit slice*)
   (*make mask an normalize data to first echo*)
   maskT2 = SmoothMask[Mask[Mean[Transpose[datan]], {2}], MaskComponents -> 2, MaskClosing -> 1]; 
   dataT2 = NormalizeData[MaskDTIdata[datan, maskT2]];
   (*create mask selecting fat*)
   fmask = Mask[dataT2[[All, -1]], {0.4}];
   fmask = ImageData[SelectComponents[Image3D[fmask], "Count", -2]];
   (*data for calibration fit*)
   fitData = Transpose[{times, Mean[Flatten[GetMaskData[#, fmask]]] & /@ Transpose[dataT2]}];
   datal = Transpose[dataT2, {1, 4, 2, 3}]
   ];
 
  model = Af *(cs*Exp[-x/T2s] + (1 - cs)*Exp[-x/T2f]) + Am * Exp[-x/T2m];
  
  (*perform callibration fit*)
 {Afi, Ami, csi, T2mi, T2fi, T2si} = {Af, Am, cs, T2m, T2f, T2s} /. 
   FindFit[fitData, {model, 
   	{0.1 <= cs <= 0.9, 30 < T2m < 50, T2m < T2f, T2f < T2s, 0 <= Am, Am < Af}},
   	{{Af, 1 fitData[[1, 2]]}, {Am, 0.25 fitData[[1, 2]]},{cs, 0.33}, {T2m, 35}, {T2f, 81}, {T2s, 250}}, 
   	x, Method -> "NMinimize"];
  (*normalize signal fractions*)  	
  sci = Afi + Ami;
  {Afi, Ami} = {Afi, Ami}/sci;
  cal = Round[{sci, {100 Ami, T2mi}, {100 csi, T2si, T2fi}}, .1];
  
  (*Print the callibration results*)
  Print[Row[{
    Column[Row /@ {
       {"signal: ", Round[sci, .1]},
       {"f mus: ", Round[100 Ami, .1] , " - T2 mus: ", 
        Round[T2mi, .1] },
       {""},
       {"f fat-slow: ", Round[100 csi, .1]}, 
       {" T2 slow: ", Round[T2si, .1], " - T2 fast: ", Round[T2fi, .1]}
       }, Alignment -> Center], 
    Show[
    	ListLinePlot[fitData, PlotStyle -> Directive[Thick, Red], PlotRange -> {{0, 150}, {0, sci}}, ImageSize -> 150], 
    	Plot[sci model /. Thread[{Af, Am, cs, T2m, T2f, T2s} -> {Afi, Ami, csi, T2mi, T2fi, T2si}], {x, 0, 200}, PlotStyle -> Directive[{Black, Dashed}]]
    	]
    	}, "   "]];
  
  (*define the fat values in the model*)
  model2 = Af *(csi*Exp[-x/T2si] + (1 - csi)*Exp[-x/T2fi]) + Am * Exp[-x/T2m];
  DistributeDefinitions[model2];
  
  (*perform the per voxel fit*)
  result = ParallelMap[(
      If[N@Total[#] === 0.,
       {0., 0., 0., 0.},
       
       (*fit the muscle T2*)
       fdat = Transpose[{times, #}];
       {Aff, Amm, T2mf} = {Af, Am, T2m} /. Quiet[FindFit[fdat, model2, {{Af, 0.125 fdat[[1, 2]]}, {Am, 1.125 fdat[[1, 2]]}, {T2m, 35}}, x]];
       
       S0 = Aff + Amm;
       {Aff, Amm} = {Aff, Amm}/S0;
       {S0, Aff, Amm, T2mf}
       
       ]) &, datal, {ad - 1}];
  
  (*generate output*)
  {S0, ffr, mfr, T2} = Switch[ad, 
  	3, Transpose[result, {2, 3, 1}], 
  	4, Transpose[result, {2, 3, 4, 1}]];
  	
  T2 = Clip[T2, {0, 500}, {0, 500}];
  S0 = Clip[S0, {Min[datan], 1.5 Max[datan]}, {0, 0}];
  mfr = Clip[mfr, {0, 1}, {0, 1}];
  ffr = Clip[ffr, {0, 1}, {0, 1}];
  
  (*in needed also output the calibration values*)
  If[OptionValue[OutputCalibration],{N@{S0, ffr, mfr, T2},cal},N@{S0, ffr, mfr, T2}]
  ]


(* ::Section:: *)
(*End Package*)


End[]

If[DTITools`verbose,Print[Names["DTITools`RelaxometryTools`*"]]];
SetAttributes[#,{Protected, ReadProtected}]&/@Names["DTITools`RelaxometryTools`*"];

EndPackage[]
